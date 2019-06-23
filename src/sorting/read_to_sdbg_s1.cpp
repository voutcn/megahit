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

#include "read_to_sdbg.h"

#include <omp.h>
#include <algorithm>
#include <mutex>

#include "sequence/kmer.h"
#include "sequence/packed_reads.h"
#include "sequence/lib_io.h"
#include "sequence/readers/kseq.h"
#include "sequence/sequence_package.h"
#include "utils/safe_open.h"
#include "utils/utils.h"

namespace { // helper functions

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

inline uint64_t ComposeUint64(uint32_t *read_info_ptr) {
  return (static_cast<uint64_t>(*read_info_ptr) << 32u) | *(read_info_ptr + 1);
}

inline uint8_t ExtractHeadTail(const uint32_t *item, int64_t spacing, int words_per_substring) {
  return *(item + spacing * (words_per_substring - 1)) & ((1 << 2 * Read2SdbgS1::kBWTCharNumBits) - 1);
}

inline uint8_t ExtractPrevNext(int64_t readinfo) { return readinfo & ((1 << 2 * Read2SdbgS1::kBWTCharNumBits) - 1); }

} // namespace

// sorting core functions

int64_t Read2SdbgS1::Lv0EncodeDiffBase(int64_t read_id) {
  return EncodeOffset(read_id, 0, 0, seq_pkg_->package);
}

Read2SdbgS1::Meta Read2SdbgS1::Initialize() {
  bool is_reverse = true;
  int64_t num_bases, num_reads;
  GetBinaryLibSize(opt_.read_lib_file, num_bases, num_reads);

  seq_pkg_->package.ReserveSequences(num_reads);
  seq_pkg_->package.ReserveBases(num_bases);

  std::vector<lib_info_t> lib_info;
  ReadBinaryLibs(opt_.read_lib_file, seq_pkg_->package, lib_info, is_reverse);

  seq_pkg_->package.BuildIndex();

  xinfo("%lu reads, %d max read length, %lld total bases\n", seq_pkg_->package.SeqCount(),
        seq_pkg_->package.MaxSequenceLength(),
        seq_pkg_->package.BaseCount());

  words_per_substr_ =
      DivCeiling((opt_.k - 1) * kBitsPerEdgeChar + 2 * kBWTCharNumBits, kBitsPerEdgeWord);
  xinfo("%d words per substring\n", words_per_substr_);

  if (opt_.solid_threshold > 1) {
    seq_pkg_->is_solid.reset(seq_pkg_->package.BaseCount());
  }


  // --- initialize output mercy files ---
  seq_pkg_->n_mercy_files = 1;

  while (seq_pkg_->n_mercy_files * 10485760LL <
      static_cast<int64_t>(seq_pkg_->package.SeqCount()) && seq_pkg_->n_mercy_files < 64) {
    seq_pkg_->n_mercy_files <<= 1;
  }
  xinfo("Number of files for mercy candidate reads: %d\n", seq_pkg_->n_mercy_files);

  for (int i = 0; i < seq_pkg_->n_mercy_files; ++i) {
    mercy_files_.push_back(xfopen(FormatString("%s.mercy_cand.%d", opt_.output_prefix.c_str(), i), "wb"));
  }

  // --- initialize stat ---
  thread_edge_counting_.resize(opt_.n_threads);
  for (auto &c : thread_edge_counting_) {
    c.resize(kMaxMul + 1);
    std::fill(c.begin(), c.end(), 0);
  }

  int64_t memory_for_data = DivCeiling(seq_pkg_->is_solid.size(), 8)
      + seq_pkg_->package.SizeInByte()
      + (kMaxMul + 1) * (opt_.n_threads + 1) * sizeof(int64_t);

  return {
      num_reads,
      memory_for_data,
      words_per_substr_ + 2,
      2,
  };
}

void Read2SdbgS1::Lv0CalcBucketSize(SeqPartition *_data) {
  auto &rp = *_data;
  auto &bucket_sizes = rp.rp_bucket_sizes;
  std::fill(bucket_sizes.begin(), bucket_sizes.end(), 0);
  GenericKmer k_minus1_mer, rev_k_minus1_mer;  // (k-1)-mer and its rc

  for (int64_t read_id = rp.from; read_id < rp.to; ++read_id) {
    auto read_length = seq_pkg_->package.SequenceLength(read_id);

    if (read_length < opt_.k + 1) {
      continue;
    }

    auto ptr_and_offset = seq_pkg_->package.WordPtrAndOffset(read_id);
    int64_t offset = ptr_and_offset.second;
    const uint32_t *read_p = ptr_and_offset.first;

    k_minus1_mer.InitFromPtr(read_p, offset, opt_.k - 1);
    rev_k_minus1_mer = k_minus1_mer;
    rev_k_minus1_mer.ReverseComplement(opt_.k - 1);

    // the first one special handling
    bucket_sizes[k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
    bucket_sizes[rev_k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;

    unsigned last_char_offset = opt_.k - 1;
    int c = seq_pkg_->package.GetBase(read_id, last_char_offset);
    k_minus1_mer.ShiftAppend(c, opt_.k - 1);
    rev_k_minus1_mer.ShiftPreappend(3 - c, opt_.k - 1);

    while (last_char_offset < read_length - 1) {
      int cmp = k_minus1_mer.cmp(rev_k_minus1_mer, opt_.k - 1);

      if (cmp > 0) {
        bucket_sizes[rev_k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
      } else {
        bucket_sizes[k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
      }

      int c = seq_pkg_->package.GetBase(read_id, ++last_char_offset);
      k_minus1_mer.ShiftAppend(c, opt_.k - 1);
      rev_k_minus1_mer.ShiftPreappend(3 - c, opt_.k - 1);
    }

    // last one special handling
    bucket_sizes[k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
    bucket_sizes[rev_k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
  }
}

void Read2SdbgS1::Lv1FillOffsets(SeqPartition *_data) {
  auto &rp = *_data;
  std::array<int64_t, kNumBuckets> prev_full_offsets{};  // temporary array for computing differentials
  for (auto b = GetLv1StartBucket(); b < GetLv1EndBucket(); ++b)
    prev_full_offsets[b] = rp.rp_lv1_differential_base;

  // this loop is VERY similar to that in PreprocessScanToFillBucketSizesThread
  GenericKmer k_minus1_mer, rev_k_minus1_mer;  // (k+1)-mer and its rc

  int key;

  for (int64_t read_id = rp.from; read_id < rp.to; ++read_id) {
    auto read_length = seq_pkg_->package.SequenceLength(read_id);

    if (read_length < opt_.k + 1) {
      continue;
    }

    auto ptr_and_offset = seq_pkg_->package.WordPtrAndOffset(read_id);
    int64_t offset = ptr_and_offset.second;
    const uint32_t *read_p = ptr_and_offset.first;

    k_minus1_mer.InitFromPtr(read_p, offset, opt_.k - 1);
    rev_k_minus1_mer = k_minus1_mer;
    rev_k_minus1_mer.ReverseComplement(opt_.k - 1);

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

    unsigned last_char_offset = opt_.k - 1;
    int c = seq_pkg_->package.GetBase(read_id, last_char_offset);
    k_minus1_mer.ShiftAppend(c, opt_.k - 1);
    rev_k_minus1_mer.ShiftPreappend(3 - c, opt_.k - 1);

    // shift the key char by char
    while (last_char_offset < read_length - 1) {
      int cmp = k_minus1_mer.cmp(rev_k_minus1_mer, opt_.k - 1);

      if (cmp > 0) {
        key = rev_k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
        CHECK_AND_SAVE_OFFSET(last_char_offset - opt_.k + 2, 1);
      } else if (cmp < 0) {
        key = k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
        CHECK_AND_SAVE_OFFSET(last_char_offset - opt_.k + 2, 0);
      } else {
        // a not-that-math-correct solution if the edge is palindrome, but works well enough
        int prev = seq_pkg_->package.GetBase(read_id, last_char_offset - (opt_.k - 1));
        int next = seq_pkg_->package.GetBase(read_id, last_char_offset + 1);

        if (prev <= 3 - next) {
          key = rev_k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
          CHECK_AND_SAVE_OFFSET(last_char_offset - opt_.k + 2, 0);
        } else {
          key = rev_k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
          CHECK_AND_SAVE_OFFSET(last_char_offset - opt_.k + 2, 1);
        }
      }

      int c = seq_pkg_->package.GetBase(read_id, ++last_char_offset);
      k_minus1_mer.ShiftAppend(c, opt_.k - 1);
      rev_k_minus1_mer.ShiftPreappend(3 - c, opt_.k - 1);
    }

    // the last one special handling
    key = k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
    CHECK_AND_SAVE_OFFSET(last_char_offset - opt_.k + 2, 0);
    key = rev_k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
    CHECK_AND_SAVE_OFFSET(last_char_offset - opt_.k + 2, 1);
  }

#undef CHECK_AND_SAVE_OFFSET
}

void Read2SdbgS1::Lv2ExtractSubString(unsigned start_bucket, unsigned end_bucket, uint32_t *substr_ptr) {
  auto offset_iterator = GetOffsetIterator(start_bucket, end_bucket);

  while (offset_iterator.HasNext()) {
    int64_t full_offset = offset_iterator.Next();

    int64_t read_id = seq_pkg_->package.GetSeqID(full_offset >> 1);
    unsigned strand = full_offset & 1;
    unsigned offset = (full_offset >> 1) - seq_pkg_->package.StartPos(read_id);
    unsigned read_length = seq_pkg_->package.SequenceLength(read_id);
    unsigned num_chars_to_copy = opt_.k - 1;
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

    if (offset + opt_.k < read_length) {
      tail = seq_pkg_->package.GetBase(read_id, offset + opt_.k - 1);
      next = seq_pkg_->package.GetBase(read_id, offset + opt_.k);
    } else {
      next = kSentinelValue;

      if (offset + opt_.k - 1 < read_length) {
        tail = seq_pkg_->package.GetBase(read_id, offset + opt_.k - 1);
      } else {
        tail = kSentinelValue;
      }
    }

    auto ptr_and_offset = seq_pkg_->package.WordPtrAndOffset(read_id);
    int start_offset = ptr_and_offset.second;
    const uint32_t *read_p = ptr_and_offset.first;
    int words_this_read = DivCeiling(start_offset + read_length, 16);
    uint64_t read_info;

    if (strand == 0) {
      CopySubstring(substr_ptr, read_p, offset + start_offset, num_chars_to_copy, 1, words_this_read,
                    words_per_substr_);
      uint32_t *last_word = substr_ptr + int64_t(words_per_substr_ - 1) * 1;
      *last_word |= (head << kBWTCharNumBits) | tail;
      read_info = (full_offset << 6) | (prev << 3) | next;
    } else {
      CopySubstringRC(substr_ptr, read_p, offset + start_offset, num_chars_to_copy, 1, words_this_read,
                      words_per_substr_);
      uint32_t *last_word = substr_ptr + int64_t(words_per_substr_ - 1) * 1;
      *last_word |= ((tail == kSentinelValue ? kSentinelValue : 3 - tail) << kBWTCharNumBits) |
          (head == kSentinelValue ? kSentinelValue : 3 - head);
      read_info = (full_offset << 6) | ((next == kSentinelValue ? kSentinelValue : (3 - next)) << 3) |
          (prev == kSentinelValue ? kSentinelValue : (3 - prev));
    }
    
    DecomposeUint64(substr_ptr + words_per_substr_, read_info);
    substr_ptr += words_per_substr_ + 2;
  }
}

void Read2SdbgS1::Lv2Postprocess(int64_t from, int64_t to, int tid, uint32_t *substr) {
  int64_t end_idx;
  int64_t count_prev_head[5][5];
  int64_t count_tail_next[5][5];
  int64_t count_head_tail[(1 << 2 * kBWTCharNumBits) - 1];
  auto &thread_edge_count = thread_edge_counting_[tid];

  for (int64_t i = from; i < to; i = end_idx) {
    end_idx = i + 1;
    uint32_t *first_item = substr + i * (2 + words_per_substr_);

    memset(count_prev_head, 0, sizeof(count_prev_head));
    memset(count_tail_next, 0, sizeof(count_tail_next));
    memset(count_head_tail, 0, sizeof(count_head_tail));

    {
      auto read_info = ComposeUint64(first_item + words_per_substr_);
      uint8_t prev_and_next = ExtractPrevNext(read_info);
      uint8_t head_and_tail = ExtractHeadTail(first_item, 1, words_per_substr_);
      count_prev_head[prev_and_next >> 3][head_and_tail >> 3]++;
      count_tail_next[head_and_tail & 7][prev_and_next & 7]++;
      count_head_tail[head_and_tail]++;
    }

    while (end_idx < to) {
      if (IsDiffKMinusOneMer(first_item, substr + end_idx * (2 + words_per_substr_), 1, opt_.k)) {
        break;
      }

      auto read_info = ComposeUint64(first_item + words_per_substr_);
      uint8_t prev_and_next = ExtractPrevNext(read_info);
      uint8_t head_and_tail =
          ExtractHeadTail(substr + end_idx * (words_per_substr_ + 2), 1, words_per_substr_);
      count_prev_head[prev_and_next >> 3][head_and_tail >> 3]++;
      count_tail_next[head_and_tail & 7][prev_and_next & 7]++;
      count_head_tail[head_and_tail]++;

      ++end_idx;
    }

    int has_in = 0, has_out = 0;

    for (int j = 0; j < 4; ++j) {
      for (int x = 0; x < 4; ++x) {
        if (count_prev_head[x][j] >= opt_.solid_threshold) {
          has_in |= 1 << j;
          break;
        }
      }

      for (int x = 0; x < 4; ++x) {
        if (count_tail_next[j][x] >= opt_.solid_threshold) {
          has_out |= 1 << j;
          break;
        }
      }
    }

    int l_has_out = 0, r_has_in = 0;

    for (int j = 0; j < 4; ++j) {
      for (int x = 0; x < 4; ++x) {
        if (count_head_tail[(j << kBWTCharNumBits) | x] >= opt_.solid_threshold) {
          l_has_out |= 1 << j;
          r_has_in |= 1 << x;
        }
      }
    }

    while (i < end_idx) {
      uint8_t head_and_tail =
          ExtractHeadTail(substr + i * (2 + words_per_substr_), 1, words_per_substr_);
      uint8_t head = head_and_tail >> 3;
      uint8_t tail = head_and_tail & 7;

      if (head != kSentinelValue && tail != kSentinelValue) {
        ++thread_edge_count[std::min(int64_t(kMaxMul), count_head_tail[head_and_tail])];
      }
      if (head != kSentinelValue && tail != kSentinelValue &&
          count_head_tail[head_and_tail] >= opt_.solid_threshold) {
        for (int64_t j = 0; j < count_head_tail[head_and_tail]; ++j, ++i) {
          auto read_info_context = ComposeUint64(substr + i * (2 + words_per_substr_) + words_per_substr_);
          int64_t read_info = read_info_context >> 6;
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
                   mercy_files_[read_id & (seq_pkg_->n_mercy_files - 1)]);
          }

          if (!(has_out & (1 << tail))) {
            // no out
            int64_t packed_mercy_cand = ((seq_pkg_->package.StartPos(read_id) + r_offset) << 2) | (2 - strand);
            fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1,
                   mercy_files_[read_id & (seq_pkg_->n_mercy_files - 1)]);
          }
        }
      } else {
        // not solid, but we still need to tell whether its left/right kmer is solid
        for (int64_t j = 0; j < count_head_tail[head_and_tail]; ++j, ++i) {
          auto read_info_context = ComposeUint64(substr + i * (2 + words_per_substr_) + words_per_substr_);
          int64_t read_info = read_info_context >> 6;
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
                     mercy_files_[read_id & (seq_pkg_->n_mercy_files - 1)]);
            } else {
              // has out but no in
              int64_t packed_mercy_cand = ((seq_pkg_->package.StartPos(read_id) + l_offset) << 2) | (1 + strand);
              fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1,
                     mercy_files_[read_id & (seq_pkg_->n_mercy_files - 1)]);
            }
          } else {
            if (has_in & (1 << head)) {
              // has in but no out
              int64_t packed_mercy_cand = ((seq_pkg_->package.StartPos(read_id) + l_offset) << 2) | (2 - strand);
              fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1,
                     mercy_files_[read_id & (seq_pkg_->n_mercy_files - 1)]);
            }
          }

          if (r_has_in & (1 << tail)) {
            if (has_out & (1 << tail)) {
              // has both in & out
              int64_t packed_mercy_cand = ((seq_pkg_->package.StartPos(read_id) + r_offset) << 2) | 0;
              fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1,
                     mercy_files_[read_id & (seq_pkg_->n_mercy_files - 1)]);
            } else {
              // has in but no out
              int64_t packed_mercy_cand = ((seq_pkg_->package.StartPos(read_id) + r_offset) << 2) | (2 - strand);
              fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1,
                     mercy_files_[read_id & (seq_pkg_->n_mercy_files - 1)]);
            }
          } else {
            if (has_out & (1 << tail)) {
              // has out but no in
              int64_t packed_mercy_cand = ((seq_pkg_->package.StartPos(read_id) + r_offset) << 2) | (1 + strand);
              fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1,
                     mercy_files_[read_id & (seq_pkg_->n_mercy_files - 1)]);
            }
          }
        }
      }
    }
  }
}

void Read2SdbgS1::Lv0Postprocess() {
  std::vector<int64_t> edge_counting(kMaxMul + 1, 0);
  for (int t = 0; t < opt_.n_threads; ++t) {
    for (int i = 1; i <= kMaxMul; ++i) {
      edge_counting[i] += thread_edge_counting_[t][i];
    }
  }

  // --- stat ---
  int64_t num_solid_edges = 0;

  for (int i = opt_.solid_threshold; i <= kMaxMul; ++i) {
    num_solid_edges += edge_counting[i];
  }

  xinfo("Total number of solid edges: %llu\n", num_solid_edges);

  FILE *counting_file = xfopen((std::string(opt_.output_prefix) + ".counting").c_str(), "w");

  for (int64_t i = 1, acc = 0; i <= kMaxMul; ++i) {
    acc += edge_counting[i];
    fprintf(counting_file, "%lld %lld\n", (long long) i, (long long) acc);
  }

  fclose(counting_file);

  // --- cleaning ---
  thread_edge_counting_ = std::vector<std::vector<int64_t>>();
  for (int i = 0; i < seq_pkg_->n_mercy_files; ++i) {
    fclose(mercy_files_[i]);
  }
}
