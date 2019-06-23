/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong & L3 Bioinformatics Limited
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* contact: Dinghua Li <dhli@cs.hku.hk> */

#include "cx1_read2sdbg.h"

#include <omp.h>
#include <string.h>
#include <zlib.h>
#include <algorithm>
#include <mutex>

#include "sequence/kmer.h"
#include "sequence/packed_reads.h"
#include "sequence/lib_io.h"
#include "sequence/readers/kseq.h"
#include "sequence/sequence_package.h"
#include "utils/safe_open.h"
#include "utils/utils.h"

#include "sorting.h"
// helping functions

namespace cx1_read2sdbg {

/**
 * @brief encode read_id and its offset in one int64_t
 */
inline int64_t EncodeOffset(int64_t read_id, int offset, int strand, SeqPackage &p) {
  return ((p.StartPos(read_id) + offset) << 1) | strand;
}

// helper: see whether two lv2 items have the same (k-1)-mer
inline bool IsDiffKMinusOneMer(uint32_t *item1, uint32_t *item2, int64_t spacing, int kmer_k) {
  // mask extra bits
  int chars_in_last_word = (kmer_k - 1) % kCharsPerEdgeWord;
  int num_full_words = (kmer_k - 1) / kCharsPerEdgeWord;

  if (chars_in_last_word > 0) {
    uint32_t w1 = item1[num_full_words * spacing];
    uint32_t w2 = item2[num_full_words * spacing];

    if ((w1 >> (kCharsPerEdgeWord - chars_in_last_word) * kBitsPerEdgeChar) !=
        (w2 >> (kCharsPerEdgeWord - chars_in_last_word) * kBitsPerEdgeChar)) {
      return true;
    }
  }

  for (int i = num_full_words - 1; i >= 0; --i) {
    if (item1[i * spacing] != item2[i * spacing]) {
      return true;
    }
  }

  return false;
}

inline uint8_t ExtractHeadTail(const uint32_t *item, int64_t spacing, int words_per_substring) {
  return *(item + spacing * (words_per_substring - 1)) & ((1 << 2 * kBWTCharNumBits) - 1);
}

inline uint8_t ExtractPrevNext(int64_t readinfo) { return readinfo & ((1 << 2 * kBWTCharNumBits) - 1); }

// cx1 core functions

int64_t CX1Read2SdbgS1::encode_lv1_diff_base_func_(int64_t read_id, int &g) {
  return EncodeOffset(read_id, 0, 0, seq_pkg_->package);
}

void CX1Read2SdbgS1::prepare_func_(int &globals) {
  bool is_reverse = true;
  int64_t num_bases, num_reads;
  GetBinaryLibSize(opt.read_lib_file, num_bases, num_reads);

  if (opt.assist_seq_file != "") {
    FILE *assist_seq_info = xfopen((opt.assist_seq_file + ".info").c_str(), "r");
    long long num_ass_bases, num_ass_seq;
    if (fscanf(assist_seq_info, "%lld%lld", &num_ass_seq, &num_ass_bases) != 2) {
      xfatal("Invalid format\n");
    }
    fclose(assist_seq_info);

    num_bases += num_ass_bases;
    num_reads += num_ass_seq;
  }

  seq_pkg_->package.ReserveSequences(num_reads);
  seq_pkg_->package.ReserveBases(num_bases);

  std::vector<lib_info_t> lib_info;
  ReadBinaryLibs(opt.read_lib_file, seq_pkg_->package, lib_info, is_reverse);
  // set up these figures before reading assist seq
  int max_read_len = seq_pkg_->package.MaxSequenceLength();

  if (opt.assist_seq_file != "") {
    FastxReader reader({opt.assist_seq_file});
    reader.ReadAll(&seq_pkg_->package, is_reverse);
  }

  seq_pkg_->package.BuildIndex();

  xinfo("%lu reads, %d max read length, %lld total bases\n", seq_pkg_->package.Size(), max_read_len,
        seq_pkg_->package.BaseCount());

  int64_t mem_packed_reads = 0;
  if (opt.kmer_freq_threshold == 1) {
    // do not need to count solid kmers
    mem_packed_reads = seq_pkg_->package.SizeInByte();
  } else {
    mem_packed_reads = DivCeiling(seq_pkg_->package.BaseCount(), 8) + seq_pkg_->package.SizeInByte();
  }

  int64_t mem_low_bound = mem_packed_reads + kNumBuckets * sizeof(int64_t) * (opt.num_cpu_threads * 3 + 1) +
      (kMaxMul + 1) * (opt.num_cpu_threads + 1) * sizeof(int64_t);
  mem_low_bound *= 1.05;

  if (mem_low_bound > opt.host_mem) {
    xfatal("%lld bytes is not enough for CX1 sorting, please set -m parameter to at least %lld\n", opt.host_mem,
           mem_low_bound);
  }

  // set cx1 param
  SetNumCpuThreads(opt.num_cpu_threads);;
  SetNumItems(seq_pkg_->package.Size());
}

void CX1Read2SdbgS1::lv0_calc_bucket_size_func_(ReadPartition *_data) {
  auto &rp = *_data;
  int &globals = *(rp.globals);
  auto &bucket_sizes = rp.rp_bucket_sizes;
  std::fill(bucket_sizes.begin(), bucket_sizes.end(), 0);
  GenericKmer k_minus1_mer, rev_k_minus1_mer;  // (k-1)-mer and its rc

  for (int64_t read_id = rp.rp_start_id; read_id < rp.rp_end_id; ++read_id) {
    int read_length = seq_pkg_->package.SequenceLength(read_id);

    if (read_length < opt.kmer_k + 1) {
      continue;
    }

    auto ptr_and_offset = seq_pkg_->package.WordPtrAndOffset(read_id);
    int64_t offset = ptr_and_offset.second;
    const uint32_t *read_p = ptr_and_offset.first;

    k_minus1_mer.InitFromPtr(read_p, offset, opt.kmer_k - 1);
    rev_k_minus1_mer = k_minus1_mer;
    rev_k_minus1_mer.ReverseComplement(opt.kmer_k - 1);

    // the first one special handling
    bucket_sizes[k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
    bucket_sizes[rev_k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;

    int last_char_offset = opt.kmer_k - 1;
    int c = seq_pkg_->package.GetBase(read_id, last_char_offset);
    k_minus1_mer.ShiftAppend(c, opt.kmer_k - 1);
    rev_k_minus1_mer.ShiftPreappend(3 - c, opt.kmer_k - 1);

    while (last_char_offset < read_length - 1) {
      int cmp = k_minus1_mer.cmp(rev_k_minus1_mer, opt.kmer_k - 1);

      if (cmp > 0) {
        bucket_sizes[rev_k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
      } else {
        bucket_sizes[k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
      }

      int c = seq_pkg_->package.GetBase(read_id, ++last_char_offset);
      k_minus1_mer.ShiftAppend(c, opt.kmer_k - 1);
      rev_k_minus1_mer.ShiftPreappend(3 - c, opt.kmer_k - 1);
    }

    // last one special handling
    bucket_sizes[k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
    bucket_sizes[rev_k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
  }
}

void CX1Read2SdbgS1::init_global_and_set_cx1_func_(int &globals) {
  int64_t max_bucket_size = *std::max_element(GetBucketSizes().begin(), GetBucketSizes().end());
  int64_t tot_bucket_size = 0;
  int num_non_empty = 0;

  for (int i = 0; i < kNumBuckets; ++i) {
    tot_bucket_size += GetBucketSizes()[i];
    num_non_empty += (GetBucketSizes()[i] > 0);
  }

  // to count (k+1)-mers, sort by the internal (k-1)-mer
  // (k+1)-mer = abS[0..k-2]cd
  // is solid: number of bSc >= threshold
  // bS has in coming: for some a, num of abS >= threshold
  // Sc has outgoing: for some a, num of Scd >= threshold
  words_per_substring =
      DivCeiling((opt.kmer_k - 1) * kBitsPerEdgeChar + 2 * kBWTCharNumBits, kBitsPerEdgeWord);
  xinfo("%d words per substring\n", words_per_substring);

  // lv2 bytes: substring, readinfo
  int64_t lv2_bytes_per_item = (words_per_substring) * sizeof(uint32_t) + sizeof(uint64_t);
  num_non_empty = std::max(1, num_non_empty);

  for (int i = 0; i < kNumBuckets; ++i) {
    if (GetBucketSizes()[i] > 2 * tot_bucket_size / num_non_empty) {
      // xinfo("Bucket %d size = %lld > %lld = 2 * avg\n", i, (long
      // long)GetBucketSizes()[i], (long long)2 * tot_bucket_size / num_non_empty);
    }
  }

  int64_t max_sorting_items =
      std::max(3 * tot_bucket_size / num_non_empty * opt.num_cpu_threads, max_bucket_size);

  int64_t mem_packed_reads = 0;
  if (opt.kmer_freq_threshold == 1) {
    // do not need to count solid kmers
    mem_packed_reads = seq_pkg_->package.SizeInByte();
  } else {
    seq_pkg_->is_solid.reset(seq_pkg_->package.BaseCount());
    mem_packed_reads = DivCeiling(seq_pkg_->package.BaseCount(), 8) + seq_pkg_->package.SizeInByte();
  }

  int64_t mem_remained = opt.host_mem - mem_packed_reads -
      opt.num_cpu_threads * 65536 * sizeof(uint64_t)  // radix sort buckets
      - kNumBuckets * sizeof(int64_t) * (opt.num_cpu_threads * 3 + 1) -
      (kMaxMul + 1) * (opt.num_cpu_threads + 1) * sizeof(int64_t);
  int64_t min_lv1_items = tot_bucket_size / (kMaxLv1ScanTime - 0.5);
  int64_t max_lv1_items = 0;

  if (opt.mem_flag == 1) {
    // auto set memory
    max_lv1_items = int64_t(tot_bucket_size / (kDefaultLv1ScanTime - 0.5));
    max_lv1_items = std::max(max_lv1_items, max_bucket_size);
    int64_t mem_needed =
        max_lv1_items * kLv1BytePerItem + max_sorting_items * lv2_bytes_per_item;

    if (mem_needed > mem_remained) {
      adjust_mem(mem_remained, lv2_bytes_per_item, min_lv1_items, max_bucket_size,
                              max_sorting_items, max_lv1_items, max_sorting_items);
    }

  } else if (opt.mem_flag == 0) {
    // min memory
    max_lv1_items = int64_t(tot_bucket_size / (kMaxLv1ScanTime - 0.5));
    max_lv1_items = std::max(max_lv1_items, max_bucket_size);
    int64_t mem_needed =
        max_lv1_items * kLv1BytePerItem + max_sorting_items * lv2_bytes_per_item;

    if (mem_needed > mem_remained) {
      adjust_mem(mem_remained, lv2_bytes_per_item, min_lv1_items, max_bucket_size,
                              max_sorting_items, max_lv1_items, max_sorting_items);
    } else {
      adjust_mem(mem_needed, lv2_bytes_per_item, min_lv1_items, max_bucket_size,
                              max_sorting_items, max_lv1_items, max_sorting_items);
    }

  } else {
    // use all
    adjust_mem(mem_remained, lv2_bytes_per_item, min_lv1_items, max_bucket_size,
                            max_sorting_items, max_lv1_items, max_sorting_items);
  }

  if (max_lv1_items < min_lv1_items) {
    xfatal("No enough memory to process.");
  }


  SetMaxLv1Lv2Items(max_lv1_items, max_sorting_items, lv2_bytes_per_item / sizeof(uint32_t), 2);
  xinfo("Memory for reads: %lld\n", mem_packed_reads);
  xinfo("max # lv.1 items = %lld\n", max_lv1_items);

  // --- initialize output mercy files ---
  seq_pkg_->n_mercy_files = 1;

  while (seq_pkg_->n_mercy_files * 10485760LL < seq_pkg_->package.Size() && seq_pkg_->n_mercy_files < 64) {
    seq_pkg_->n_mercy_files <<= 1;
  }

  xinfo("Number of files for mercy candidate reads: %d\n", seq_pkg_->n_mercy_files);

  for (int i = 0; i < seq_pkg_->n_mercy_files; ++i) {
    mercy_files.push_back(xfopen(FormatString("%s.mercy_cand.%d", opt.output_prefix.c_str(), i), "wb"));
  }

  // --- initialize stat ---
  thread_edge_counting.resize(opt.num_cpu_threads);
  for (auto &c : thread_edge_counting) {
    c.resize(kMaxMul + 1);
    std::fill(c.begin(), c.end(), 0);
  }
}

void CX1Read2SdbgS1::lv1_fill_offset_func_(ReadPartition *_data) {
  auto &rp = *_data;
  std::array<int64_t, kNumBuckets> prev_full_offsets{};  // temporary array for computing differentials
  for (auto b = GetLv1StartBucket(); b < GetLv1EndBucket(); ++b)
    prev_full_offsets[b] = rp.rp_lv1_differential_base;

  // this loop is VERY similar to that in PreprocessScanToFillBucketSizesThread
  GenericKmer k_minus1_mer, rev_k_minus1_mer;  // (k+1)-mer and its rc

  int key;

  for (int64_t read_id = rp.rp_start_id; read_id < rp.rp_end_id; ++read_id) {
    int read_length = seq_pkg_->package.SequenceLength(read_id);

    if (read_length < opt.kmer_k + 1) {
      continue;
    }

    auto ptr_and_offset = seq_pkg_->package.WordPtrAndOffset(read_id);
    int64_t offset = ptr_and_offset.second;
    const uint32_t *read_p = ptr_and_offset.first;

    k_minus1_mer.InitFromPtr(read_p, offset, opt.kmer_k - 1);
    rev_k_minus1_mer = k_minus1_mer;
    rev_k_minus1_mer.ReverseComplement(opt.kmer_k - 1);

    // ===== this is a macro to save some copy&paste ================
#define CHECK_AND_SAVE_OFFSET(offset, strand)                                                         \
  do {                                                                                                \
    if (HandlingBucket(key)) {                                                           \
      int key_ = GetBucketRank(key);                                                     \
      int64_t full_offset = EncodeOffset(read_id, offset, strand, seq_pkg_->package);                   \
      int64_t differential = full_offset - prev_full_offsets[key_];                                   \
      int64_t index = rp.rp_bucket_offsets[key_]++;                                                   \
      WriteOffset(index, differential, full_offset);                                     \
      assert(rp.rp_bucket_offsets[key_] <= GetLv1NumItems());                            \
      prev_full_offsets[key_] = full_offset;                                                          \
    }                                                                                                 \
  } while (0)
    // ^^^^^ why is the macro surrounded by a do-while? please ask Google
    // =========== end macro ==========================

    // the first one special handling
    key = k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
    CHECK_AND_SAVE_OFFSET(0, 0);
    key = rev_k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
    CHECK_AND_SAVE_OFFSET(0, 1);

    int last_char_offset = opt.kmer_k - 1;
    int c = seq_pkg_->package.GetBase(read_id, last_char_offset);
    k_minus1_mer.ShiftAppend(c, opt.kmer_k - 1);
    rev_k_minus1_mer.ShiftPreappend(3 - c, opt.kmer_k - 1);

    // shift the key char by char
    while (last_char_offset < read_length - 1) {
      int cmp = k_minus1_mer.cmp(rev_k_minus1_mer, opt.kmer_k - 1);

      if (cmp > 0) {
        key = rev_k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
        CHECK_AND_SAVE_OFFSET(last_char_offset - opt.kmer_k + 2, 1);
      } else if (cmp < 0) {
        key = k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
        CHECK_AND_SAVE_OFFSET(last_char_offset - opt.kmer_k + 2, 0);
      } else {
        // a not-that-math-correct solution if the edge is palindrome, but works well enough
        int prev = seq_pkg_->package.GetBase(read_id, last_char_offset - (opt.kmer_k - 1));
        int next = seq_pkg_->package.GetBase(read_id, last_char_offset + 1);

        if (prev <= 3 - next) {
          key = rev_k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
          CHECK_AND_SAVE_OFFSET(last_char_offset - opt.kmer_k + 2, 0);
        } else {
          key = rev_k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
          CHECK_AND_SAVE_OFFSET(last_char_offset - opt.kmer_k + 2, 1);
        }
      }

      int c = seq_pkg_->package.GetBase(read_id, ++last_char_offset);
      k_minus1_mer.ShiftAppend(c, opt.kmer_k - 1);
      rev_k_minus1_mer.ShiftPreappend(3 - c, opt.kmer_k - 1);
    }

    // the last one special handling
    key = k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
    CHECK_AND_SAVE_OFFSET(last_char_offset - opt.kmer_k + 2, 0);
    key = rev_k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
    CHECK_AND_SAVE_OFFSET(last_char_offset - opt.kmer_k + 2, 1);
  }

#undef CHECK_AND_SAVE_OFFSET
}

void CX1Read2SdbgS1::lv2_extract_substr_(unsigned bp_from, unsigned bp_to, int &globals, uint32_t *substr) {
  auto lv1_p = GetLv1Iterator(bp_from);

  for (auto b = bp_from; b < bp_to; ++b) {
    for (int t = 0; t < opt.num_cpu_threads; ++t) {
      int64_t full_offset = GetReadPartition(t).rp_lv1_differential_base;
      int64_t num = GetReadPartition(t).rp_bucket_sizes[b];

      for (int64_t i = 0; i < num; ++i) {
        if (*lv1_p >= 0) {
          full_offset += *(lv1_p++);
        } else {
          full_offset = GetSpecialOffset(-1 - *(lv1_p++));
        }

        int64_t read_id = seq_pkg_->package.GetSeqID(full_offset >> 1);
        int strand = full_offset & 1;
        int offset = (full_offset >> 1) - seq_pkg_->package.StartPos(read_id);
        int read_length = seq_pkg_->package.SequenceLength(read_id);
        int num_chars_to_copy = opt.kmer_k - 1;
        unsigned char prev, next, head, tail;  // (k+1)=abScd, prev=a, head=b, tail=c, next=d

        assert(offset < read_length);

        if (offset > 1) {
          head = seq_pkg_->package.GetBase(read_id, offset - 1);
          prev = seq_pkg_->package.GetBase(read_id, offset - 2);
        } else {
          prev = kSentinelValue;

          if (offset > 0) {
            head = seq_pkg_->package.GetBase(read_id, offset - 1);
          } else {
            head = kSentinelValue;
          }
        }

        if (offset + opt.kmer_k < read_length) {
          tail = seq_pkg_->package.GetBase(read_id, offset + opt.kmer_k - 1);
          next = seq_pkg_->package.GetBase(read_id, offset + opt.kmer_k);
        } else {
          next = kSentinelValue;

          if (offset + opt.kmer_k - 1 < read_length) {
            tail = seq_pkg_->package.GetBase(read_id, offset + opt.kmer_k - 1);
          } else {
            tail = kSentinelValue;
          }
        }

        auto ptr_and_offset = seq_pkg_->package.WordPtrAndOffset(read_id);
        int start_offset = ptr_and_offset.second;
        const uint32_t *read_p = ptr_and_offset.first;
        int words_this_read = DivCeiling(start_offset + read_length, 16);
        auto readinfo_ptr = reinterpret_cast<uint64_t *>(substr + words_per_substring);

        if (strand == 0) {
          CopySubstring(substr, read_p, offset + start_offset, num_chars_to_copy, 1, words_this_read,
                        words_per_substring);
          uint32_t *last_word = substr + int64_t(words_per_substring - 1) * 1;
          *last_word |= (head << kBWTCharNumBits) | tail;
          *readinfo_ptr = (full_offset << 6) | (prev << 3) | next;
        } else {
          CopySubstringRC(substr, read_p, offset + start_offset, num_chars_to_copy, 1, words_this_read,
                          words_per_substring);
          uint32_t *last_word = substr + int64_t(words_per_substring - 1) * 1;
          *last_word |= ((tail == kSentinelValue ? kSentinelValue : 3 - tail) << kBWTCharNumBits) |
              (head == kSentinelValue ? kSentinelValue : 3 - head);
          *readinfo_ptr = (full_offset << 6) | ((next == kSentinelValue ? kSentinelValue : (3 - next)) << 3) |
              (prev == kSentinelValue ? kSentinelValue : (3 - prev));
        }

        substr += words_per_substring + 2;
      }
    }
  }
}

void CX1Read2SdbgS1::output_(int64_t from, int64_t to, int tid, int &globals, uint32_t *substr) {
  int64_t end_idx;
  int64_t count_prev_head[5][5];
  int64_t count_tail_next[5][5];
  int64_t count_head_tail[(1 << 2 * kBWTCharNumBits) - 1];
  auto &thread_edge_count = thread_edge_counting[tid];

  for (int64_t i = from; i < to; i = end_idx) {
    end_idx = i + 1;
    uint32_t *first_item = substr + i * (2 + words_per_substring);

    memset(count_prev_head, 0, sizeof(count_prev_head));
    memset(count_tail_next, 0, sizeof(count_tail_next));
    memset(count_head_tail, 0, sizeof(count_head_tail));

    {
      auto readinfo_ptr = reinterpret_cast<uint64_t *>(first_item + words_per_substring);
      uint8_t prev_and_next = ExtractPrevNext(*readinfo_ptr);
      uint8_t head_and_tail = ExtractHeadTail(first_item, 1, words_per_substring);
      count_prev_head[prev_and_next >> 3][head_and_tail >> 3]++;
      count_tail_next[head_and_tail & 7][prev_and_next & 7]++;
      count_head_tail[head_and_tail]++;
    }

    while (end_idx < to) {
      if (IsDiffKMinusOneMer(first_item, substr + end_idx * (2 + words_per_substring), 1, opt.kmer_k)) {
        break;
      }

      auto readinfo_ptr = reinterpret_cast<uint64_t *>(substr + end_idx * (2 + words_per_substring) +
          words_per_substring);
      uint8_t prev_and_next = ExtractPrevNext(*readinfo_ptr);
      uint8_t head_and_tail =
          ExtractHeadTail(substr + end_idx * (words_per_substring + 2), 1, words_per_substring);
      count_prev_head[prev_and_next >> 3][head_and_tail >> 3]++;
      count_tail_next[head_and_tail & 7][prev_and_next & 7]++;
      count_head_tail[head_and_tail]++;

      ++end_idx;
    }

    int has_in = 0, has_out = 0;

    for (int j = 0; j < 4; ++j) {
      for (int x = 0; x < 4; ++x) {
        if (count_prev_head[x][j] >= opt.kmer_freq_threshold) {
          has_in |= 1 << j;
          break;
        }
      }

      for (int x = 0; x < 4; ++x) {
        if (count_tail_next[j][x] >= opt.kmer_freq_threshold) {
          has_out |= 1 << j;
          break;
        }
      }
    }

    int l_has_out = 0, r_has_in = 0;

    for (int j = 0; j < 4; ++j) {
      for (int x = 0; x < 4; ++x) {
        if (count_head_tail[(j << kBWTCharNumBits) | x] >= opt.kmer_freq_threshold) {
          l_has_out |= 1 << j;
          r_has_in |= 1 << x;
        }
      }
    }

    while (i < end_idx) {
      uint8_t head_and_tail =
          ExtractHeadTail(substr + i * (2 + words_per_substring), 1, words_per_substring);
      uint8_t head = head_and_tail >> 3;
      uint8_t tail = head_and_tail & 7;

      if (head != kSentinelValue && tail != kSentinelValue) {
        ++thread_edge_count[std::min(int64_t(kMaxMul), count_head_tail[head_and_tail])];
      }
      if (head != kSentinelValue && tail != kSentinelValue &&
          count_head_tail[head_and_tail] >= opt.kmer_freq_threshold) {
        for (int64_t j = 0; j < count_head_tail[head_and_tail]; ++j, ++i) {
          auto readinfo_ptr = reinterpret_cast<uint64_t *>(substr + i * (2 + words_per_substring) +
              words_per_substring);
          int64_t read_info = *readinfo_ptr >> 6;
          int strand = read_info & 1;
          int64_t read_id = seq_pkg_->package.GetSeqID(read_info >> 1);
          int offset = (read_info >> 1) - seq_pkg_->package.StartPos(read_id) - 1;
          int l_offset = strand == 0 ? offset : offset + 1;
          int r_offset = strand == 0 ? offset + 1 : offset;

          // mark this is a solid edge
          seq_pkg_->is_solid.set((read_info >> 1) - 1);

          if (!(has_in & (1 << head))) {
            // no in
            int64_t packed_mercy_cand = ((seq_pkg_->package.StartPos(read_id) + l_offset) << 2) | (1 + strand);
            fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1,
                   mercy_files[read_id & (seq_pkg_->n_mercy_files - 1)]);
          }

          if (!(has_out & (1 << tail))) {
            // no out
            int64_t packed_mercy_cand = ((seq_pkg_->package.StartPos(read_id) + r_offset) << 2) | (2 - strand);
            fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1,
                   mercy_files[read_id & (seq_pkg_->n_mercy_files - 1)]);
          }
        }
      } else {
        // not solid, but we still need to tell whether its left/right kmer is solid
        for (int64_t j = 0; j < count_head_tail[head_and_tail]; ++j, ++i) {
          uint64_t *readinfo_ptr = reinterpret_cast<uint64_t *>(substr + i * (2 + words_per_substring) +
              words_per_substring);
          int64_t read_info = *readinfo_ptr >> 6;
          int strand = read_info & 1;
          int64_t read_id = seq_pkg_->package.GetSeqID(read_info >> 1);
          int offset = (read_info >> 1) - seq_pkg_->package.StartPos(read_id) - 1;
          int l_offset = strand == 0 ? offset : offset + 1;
          int r_offset = strand == 0 ? offset + 1 : offset;

          if (l_has_out & (1 << head)) {
            if (has_in & (1 << head)) {
              // has both in & out
              int64_t packed_mercy_cand = ((seq_pkg_->package.StartPos(read_id) + l_offset) << 2) | 0;
              fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1,
                     mercy_files[read_id & (seq_pkg_->n_mercy_files - 1)]);
            } else {
              // has out but no in
              int64_t packed_mercy_cand = ((seq_pkg_->package.StartPos(read_id) + l_offset) << 2) | (1 + strand);
              fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1,
                     mercy_files[read_id & (seq_pkg_->n_mercy_files - 1)]);
            }
          } else {
            if (has_in & (1 << head)) {
              // has in but no out
              int64_t packed_mercy_cand = ((seq_pkg_->package.StartPos(read_id) + l_offset) << 2) | (2 - strand);
              fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1,
                     mercy_files[read_id & (seq_pkg_->n_mercy_files - 1)]);
            }
          }

          if (r_has_in & (1 << tail)) {
            if (has_out & (1 << tail)) {
              // has both in & out
              int64_t packed_mercy_cand = ((seq_pkg_->package.StartPos(read_id) + r_offset) << 2) | 0;
              fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1,
                     mercy_files[read_id & (seq_pkg_->n_mercy_files - 1)]);
            } else {
              // has in but no out
              int64_t packed_mercy_cand = ((seq_pkg_->package.StartPos(read_id) + r_offset) << 2) | (2 - strand);
              fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1,
                     mercy_files[read_id & (seq_pkg_->n_mercy_files - 1)]);
            }
          } else {
            if (has_out & (1 << tail)) {
              // has out but no in
              int64_t packed_mercy_cand = ((seq_pkg_->package.StartPos(read_id) + r_offset) << 2) | (1 + strand);
              fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1,
                     mercy_files[read_id & (seq_pkg_->n_mercy_files - 1)]);
            }
          }
        }
      }
    }
  }
}

void CX1Read2SdbgS1::post_proc_func_(int &globals) {
  std::vector<int64_t> edge_counting(kMaxMul + 1, 0);
  for (int t = 0; t < opt.num_cpu_threads; ++t) {
    for (int i = 1; i <= kMaxMul; ++i) {
      edge_counting[i] += thread_edge_counting[t][i];
    }
  }

  // --- stat ---
  int64_t num_solid_edges = 0;

  for (int i = opt.kmer_freq_threshold; i <= kMaxMul; ++i) {
    num_solid_edges += edge_counting[i];
  }

  xinfo("Total number of solid edges: %llu\n", num_solid_edges);

  FILE *counting_file = xfopen((std::string(opt.output_prefix) + ".counting").c_str(), "w");

  for (int64_t i = 1, acc = 0; i <= kMaxMul; ++i) {
    acc += edge_counting[i];
    fprintf(counting_file, "%lld %lld\n", (long long) i, (long long) acc);
  }

  fclose(counting_file);

  // --- cleaning ---
  thread_edge_counting = std::vector<std::vector<int64_t>>();
  for (int i = 0; i < seq_pkg_->n_mercy_files; ++i) {
    fclose(mercy_files[i]);
  }
}

}  // namespace cx1_read2sdbg