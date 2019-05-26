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
#include <mutex>
#include <parallel/algorithm>
#include <string>
#include <vector>

#include "sequence/kmer.h"
#include "sequence/packed_reads.h"
#include "sequence/lib_io.h"
#include "utils/safe_open.h"
#include "utils/utils.h"

#include "sorting.h"

namespace cx1_read2sdbg {

namespace s2 {

typedef CX1<read2sdbg_global_t, kNumBuckets> cx1_t;
typedef CX1<read2sdbg_global_t, kNumBuckets>::ReadPartition readpartition_data_t;

// helper functions
inline int64_t EncodeOffset(int64_t read_id, int offset, int strand, SeqPackage &p, int edge_type) {
  // edge_type: 0 left $; 1 solid; 2 right $
  return ((p.StartPos(read_id) + offset) << 3) | (edge_type << 1) | strand;
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

// bS'a
inline int Extract_a(uint32_t *item, int num_words, int64_t spacing, int kmer_k) {
  int non_dollar = (item[(num_words - 1) * spacing] >> kBWTCharNumBits) & 1;

  if (non_dollar) {
    int which_word = (kmer_k - 1) / kCharsPerEdgeWord;
    int word_index = (kmer_k - 1) % kCharsPerEdgeWord;
    return (item[which_word * spacing] >> (kCharsPerEdgeWord - 1 - word_index) * kBitsPerEdgeChar) & kEdgeCharMask;
  } else {
    return kSentinelValue;
  }
}

inline int Extract_b(uint32_t *item, int num_words, int64_t spacing) {
  return item[(num_words - 1) * spacing] & ((1 << kBWTCharNumBits) - 1);
}

// cx1 core functions
int64_t s2_encode_lv1_diff_base(int64_t read_id, read2sdbg_global_t &globals) {
  return EncodeOffset(read_id, 0, 0, globals.package, 0);
}

void s2_read_mercy_prepare(read2sdbg_global_t &globals) {
  if (!globals.need_mercy || globals.kmer_freq_threshold == 1) return;

  SimpleTimer timer;

  timer.reset();
  timer.start();
  xinfo("Adding mercy edges...\n");

  std::vector<uint64_t> mercy_cand;
  uint64_t num_mercy = 0;
  AtomicBitVector read_marker;
  read_marker.reset(globals.num_short_reads);

  for (int fid = 0; fid < globals.num_mercy_files; ++fid) {
    FILE *fp = xfopen(FormatString("%s.mercy_cand.%d", globals.output_prefix.c_str(), fid), "rb");
    mercy_cand.clear();

    int num_read = 0;
    uint64_t buf[4096];

    while ((num_read = fread(buf, sizeof(uint64_t), 4096, fp)) > 0) {
      mercy_cand.insert(mercy_cand.end(), buf, buf + num_read);
    }
    xinfo("Mercy file: %s, %lu\n", FormatString("%s.mercy_cand.%d", globals.output_prefix.c_str(), fid),
          mercy_cand.size());

    omp_set_num_threads(globals.num_cpu_threads);
    __gnu_parallel::sort(mercy_cand.begin(), mercy_cand.end());

    // multi threading
    uint64_t avg = DivCeiling(mercy_cand.size(), globals.num_cpu_threads);
    std::vector<uint64_t> start_idx(globals.num_cpu_threads);
    std::vector<uint64_t> end_idx(globals.num_cpu_threads);

    // manually distribute threads
    for (int tid = 0; tid < globals.num_cpu_threads; ++tid) {
      if (tid == 0) {
        start_idx[tid] = 0;
      } else {
        start_idx[tid] = end_idx[tid - 1];
      }

      uint64_t this_end = std::min(start_idx[tid] + avg, (uint64_t)mercy_cand.size());
      uint64_t read_id = 0;
      if (this_end < mercy_cand.size()) {
        read_id = globals.package.GetSeqID(mercy_cand[this_end] >> 2);
      }

      while (this_end < mercy_cand.size() && globals.package.GetSeqID(mercy_cand[this_end] >> 2) == read_id) {
        ++this_end;
      }

      end_idx[tid] = this_end;
    }

#pragma omp parallel for reduction(+ : num_mercy)

    for (int tid = 0; tid < globals.num_cpu_threads; ++tid) {
      std::vector<bool> no_in(globals.max_read_length);
      std::vector<bool> no_out(globals.max_read_length);
      std::vector<bool> has_solid_kmer(globals.max_read_length);

      uint64_t i = start_idx[tid];

      // go read by read
      while (i != end_idx[tid]) {
        uint64_t read_id = globals.package.GetSeqID(mercy_cand[i] >> 2);
        assert(!read_marker.at(read_id));
        read_marker.set(read_id);
        int first_0_out = globals.max_read_length + 1;
        int last_0_in = -1;

        std::fill(no_in.begin(), no_in.end(), false);
        std::fill(no_out.begin(), no_out.end(), false);
        std::fill(has_solid_kmer.begin(), has_solid_kmer.end(), false);

        while (i != end_idx[tid] && globals.package.GetSeqID(mercy_cand[i] >> 2) == read_id) {
          int offset = (mercy_cand[i] >> 2) - globals.package.StartPos(read_id);
          if ((mercy_cand[i] & 3) == 2) {
            no_out[offset] = true;
            first_0_out = std::min(first_0_out, offset);
          } else if ((mercy_cand[i] & 3) == 1) {
            no_in[offset] = true;
            last_0_in = std::max(last_0_in, offset);
          }

          has_solid_kmer[offset] = true;
          ++i;
        }

        if (last_0_in < first_0_out) {
          continue;
        }

        int read_length = globals.package.SequenceLength(read_id);
        int last_no_out = -1;

        for (int i = 0; i + globals.kmer_k < read_length; ++i) {
          if (globals.is_solid.at(globals.package.StartPos(read_id) + i)) {
            has_solid_kmer[i] = has_solid_kmer[i + 1] = true;
          }
        }

        for (int i = 0; i + globals.kmer_k <= read_length; ++i) {
          if (no_in[i] && last_no_out != -1) {
            for (int j = last_no_out; j < i; ++j) {
              globals.is_solid.set(globals.package.StartPos(read_id) + j);
            }

            num_mercy += i - last_no_out;
          }

          if (has_solid_kmer[i]) {
            last_no_out = -1;
          }

          if (no_out[i]) {
            last_no_out = i;
          }
        }
      }
    }

    fclose(fp);
  }

  timer.stop();
  xinfo("Adding mercy Done. Time elapsed: %.4lf\n", timer.elapsed());
  xinfo("Number mercy: %llu\n", (unsigned long long)num_mercy);

  // set cx1 param
  globals.cx1.num_cpu_threads_ = globals.num_cpu_threads;
  globals.cx1.num_items_ = globals.package.Size();
}

void *s2_lv0_calc_bucket_size(void *_data) {
  readpartition_data_t &rp = *((readpartition_data_t *)_data);
  read2sdbg_global_t &globals = *(rp.globals);
  auto &bucket_sizes = rp.rp_bucket_sizes;
  std::fill(bucket_sizes.begin(), bucket_sizes.end(), 0);
  GenericKmer edge, rev_edge;  // (k+1)-mer and its rc

  for (int64_t read_id = rp.rp_start_id; read_id < rp.rp_end_id; ++read_id) {
    int read_length = globals.package.SequenceLength(read_id);

    if (read_length < globals.kmer_k + 1) {
      continue;
    }

    auto ptr_and_offset = globals.package.WordPtrAndOffset(read_id);
    int64_t offset = ptr_and_offset.second;
    const uint32_t *read_p = ptr_and_offset.first;

    edge.InitFromPtr(read_p, offset, globals.kmer_k + 1);
    rev_edge = edge;
    rev_edge.ReverseComplement(globals.kmer_k + 1);

    int last_char_offset = globals.kmer_k;
    int64_t full_offset = globals.package.StartPos(read_id);
    bool is_solid = globals.kmer_freq_threshold == 1 || read_id >= globals.num_short_reads;

    while (true) {
      if (is_solid || globals.is_solid.at(full_offset)) {
        bool is_palindrome = rev_edge.cmp(edge, globals.kmer_k + 1) == 0;
        bucket_sizes[(edge.data()[0] << 2) >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;

        if (!is_palindrome)
          bucket_sizes[(rev_edge.data()[0] << 2) >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;

        if (last_char_offset == globals.kmer_k || !(is_solid || globals.is_solid.at(full_offset - 1))) {
          bucket_sizes[edge.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;

          if (!is_palindrome)
            bucket_sizes[(rev_edge.data()[0] << 4) >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
        }

        if (last_char_offset == read_length - 1 || !(is_solid || globals.is_solid.at(full_offset + 1))) {
          bucket_sizes[(edge.data()[0] << 4) >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;

          if (!is_palindrome)
            bucket_sizes[rev_edge.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
        }
      }

      ++full_offset;

      if (++last_char_offset >= read_length) {
        break;
      } else {
        int c = globals.package.GetBase(read_id, last_char_offset);
        edge.ShiftAppend(c, globals.kmer_k + 1);
        rev_edge.ShiftPreappend(3 - c, globals.kmer_k + 1);
      }
    }
  }

  return NULL;
}

void s2_init_global_and_set_cx1(read2sdbg_global_t &globals) {
  globals.max_bucket_size = *std::max_element(globals.cx1.bucket_sizes_.begin(), globals.cx1.bucket_sizes_.end());
  globals.tot_bucket_size = 0;
  int num_non_empty = 0;

  for (int i = 0; i < kNumBuckets; ++i) {
    globals.tot_bucket_size += globals.cx1.bucket_sizes_[i];

    if (globals.cx1.bucket_sizes_[i] > 0) {
      num_non_empty++;
    }
  }

  globals.words_per_substring = DivCeiling(globals.kmer_k * kBitsPerEdgeChar + kBWTCharNumBits + 1, kBitsPerEdgeWord);
  globals.words_per_dummy_node = DivCeiling(globals.kmer_k * kBitsPerEdgeChar, kBitsPerEdgeWord);

  xinfo("%d words per substring, words per dummy node ($v): %d\n", globals.words_per_substring,
        globals.words_per_dummy_node);

  // --- calculate lv2 memory ---

  num_non_empty = std::max(1, num_non_empty);

  int64_t lv2_bytes_per_item = globals.words_per_substring * sizeof(uint32_t);

  globals.max_sorting_items =
      std::max(globals.tot_bucket_size * globals.num_cpu_threads / num_non_empty, globals.max_bucket_size);
  globals.num_output_threads = globals.num_cpu_threads;

  int64_t mem_remained = globals.host_mem - globals.mem_packed_reads -
                         globals.num_cpu_threads * 65536 * sizeof(uint64_t)  // radix sort buckets
                         - kNumBuckets * sizeof(int64_t) * (globals.num_cpu_threads * 3 + 1);
  int64_t min_lv1_items = globals.tot_bucket_size / (kMaxLv1ScanTime - 0.5);

  if (globals.mem_flag == 1) {
    // auto set memory
    globals.cx1.max_lv1_items_ = int64_t(globals.tot_bucket_size / (kDefaultLv1ScanTime - 0.5));
    globals.cx1.max_lv1_items_ = std::max(globals.cx1.max_lv1_items_, globals.max_bucket_size);
    int64_t mem_needed =
        globals.cx1.max_lv1_items_ * cx1_t::kLv1BytePerItem + globals.max_sorting_items * lv2_bytes_per_item;

    if (mem_needed > mem_remained) {
      globals.cx1.adjust_mem(mem_remained, lv2_bytes_per_item, min_lv1_items, globals.max_bucket_size,
                             globals.max_sorting_items, globals.cx1.max_lv1_items_, globals.max_sorting_items);
    }

  } else if (globals.mem_flag == 0) {
    // min memory
    globals.cx1.max_lv1_items_ = int64_t(globals.tot_bucket_size / (kMaxLv1ScanTime - 0.5));
    globals.cx1.max_lv1_items_ = std::max(globals.cx1.max_lv1_items_, globals.max_bucket_size);
    int64_t mem_needed =
        globals.cx1.max_lv1_items_ * cx1_t::kLv1BytePerItem + globals.max_sorting_items * lv2_bytes_per_item;

    if (mem_needed > mem_remained) {
      globals.cx1.adjust_mem(mem_remained, lv2_bytes_per_item, min_lv1_items, globals.max_bucket_size,
                             globals.max_sorting_items, globals.cx1.max_lv1_items_, globals.max_sorting_items);
    } else {
      globals.cx1.adjust_mem(mem_needed, lv2_bytes_per_item, min_lv1_items, globals.max_bucket_size,
                             globals.max_sorting_items, globals.cx1.max_lv1_items_, globals.max_sorting_items);
    }

  } else {
    // use all
    globals.cx1.adjust_mem(mem_remained, lv2_bytes_per_item, min_lv1_items, globals.max_bucket_size,
                           globals.max_sorting_items, globals.cx1.max_lv1_items_, globals.max_sorting_items);
  }

  if (globals.cx1.max_lv1_items_ < min_lv1_items) {
    xfatal("No enough memory to process.");
  }

  globals.cx1.max_mem_remain_ =
      globals.cx1.max_lv1_items_ * sizeof(int) + globals.max_sorting_items * lv2_bytes_per_item;
  globals.cx1.bytes_per_sorting_item_ = lv2_bytes_per_item;
  globals.lv1_items.resize(globals.cx1.max_lv1_items_ +
                           globals.max_sorting_items * lv2_bytes_per_item / sizeof(int32_t));

  xinfo("Memory for sequence: %lld\n", globals.mem_packed_reads);
  xinfo("max # lv.1 items = %lld\n", globals.cx1.max_lv1_items_);
  // --- init output ---
  globals.sdbg_writer.set_num_threads(globals.num_output_threads);
  globals.sdbg_writer.set_kmer_size(globals.kmer_k);
  globals.sdbg_writer.set_num_buckets(kNumBuckets);
  globals.sdbg_writer.set_file_prefix(globals.output_prefix);
  globals.sdbg_writer.InitFiles();
}

void *s2_lv1_fill_offset(void *_data) {
  readpartition_data_t &rp = *((readpartition_data_t *)_data);
  read2sdbg_global_t &globals = *(rp.globals);
  std::array<int64_t, kNumBuckets> prev_full_offsets{};  // temporary array for computing differentials

  for (auto b = globals.cx1.lv1_start_bucket_; b < globals.cx1.lv1_end_bucket_; ++b)
    prev_full_offsets[b] = rp.rp_lv1_differential_base;

  // this loop is VERY similar to that in PreprocessScanToFillBucketSizesThread
  GenericKmer edge, rev_edge;  // (k+1)-mer and its rc

  int key;

  for (int64_t read_id = rp.rp_start_id; read_id < rp.rp_end_id; ++read_id) {
    int read_length = globals.package.SequenceLength(read_id);

    if (read_length < globals.kmer_k + 1) {
      continue;
    }

    auto ptr_and_offset = globals.package.WordPtrAndOffset(read_id);
    int64_t offset = ptr_and_offset.second;
    const uint32_t *read_p = ptr_and_offset.first;

    edge.InitFromPtr(read_p, offset, globals.kmer_k + 1);
    rev_edge = edge;
    rev_edge.ReverseComplement(globals.kmer_k + 1);

    // ===== this is a macro to save some copy&paste ================
#define CHECK_AND_SAVE_OFFSET(offset, strand, edge_type)                                              \
  do {                                                                                                \
    if (globals.cx1.cur_lv1_buckets_[key]) {                                                          \
      int key_ = globals.cx1.bucket_rank_[key];                                                       \
      int64_t full_offset = EncodeOffset(read_id, offset, strand, globals.package, edge_type);        \
      int64_t differential = full_offset - prev_full_offsets[key_];                                   \
      if (differential > cx1_t::kDifferentialLimit) {                                                 \
        std::lock_guard<std::mutex> lk(globals.lv1_items_scanning_lock);                              \
        globals.lv1_items[rp.rp_bucket_offsets[key_]++] = -globals.cx1.lv1_items_special_.size() - 1; \
        globals.cx1.lv1_items_special_.push_back(full_offset);                                        \
      } else {                                                                                        \
        assert((int)differential >= 0);                                                               \
        globals.lv1_items[rp.rp_bucket_offsets[key_]++] = (int)differential;                          \
      }                                                                                               \
      prev_full_offsets[key_] = full_offset;                                                          \
    }                                                                                                 \
  } while (0)
    // ^^^^^ why is the macro surrounded by a do-while? please ask Google
    // =========== end macro ==========================

    // shift the key char by char
    int last_char_offset = globals.kmer_k;
    int64_t full_offset = globals.package.StartPos(read_id);
    bool is_solid = globals.kmer_freq_threshold == 1 || read_id >= globals.num_short_reads;

    while (true) {
      if (is_solid || globals.is_solid.at(full_offset)) {
        bool is_palindrome = rev_edge.cmp(edge, globals.kmer_k + 1) == 0;

        // left $
        if (last_char_offset == globals.kmer_k || !(is_solid || globals.is_solid.at(full_offset - 1))) {
          key = edge.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
          CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 0, 0);

          if (!is_palindrome) {
            key = (rev_edge.data()[0] << 4) >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
            CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 1, 0);
          }
        }

        // solid
        key = (edge.data()[0] << 2) >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
        CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 0, 1);

        if (!is_palindrome) {
          key = (rev_edge.data()[0] << 2) >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
          CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 1, 1);
        }

        // right $
        if (last_char_offset == read_length - 1 || !(is_solid || globals.is_solid.at(full_offset + 1))) {
          key = (edge.data()[0] << 4) >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
          CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 0, 2);

          if (!is_palindrome) {
            key = rev_edge.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
            CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 1, 2);
          }
        }
      }

      ++full_offset;

      if (++last_char_offset >= read_length) {
        break;
      } else {
        int c = globals.package.GetBase(read_id, last_char_offset);
        edge.ShiftAppend(c, globals.kmer_k + 1);
        rev_edge.ShiftPreappend(3 - c, globals.kmer_k + 1);
      }
    }
  }

#undef CHECK_AND_SAVE_OFFSET
  return nullptr;
}

void s2_lv2_extract_substr_(int bp_from, int bp_to, read2sdbg_global_t &globals, uint32_t *substr) {
  auto lv1_p = globals.lv1_items.begin() + globals.cx1.rp_[0].rp_bucket_offsets[bp_from];

  for (auto b = bp_from; b < bp_to; ++b) {
    for (int t = 0; t < globals.num_cpu_threads; ++t) {
      int64_t full_offset = globals.cx1.rp_[t].rp_lv1_differential_base;
      int64_t num = globals.cx1.rp_[t].rp_bucket_sizes[b];

      for (int64_t i = 0; i < num; ++i) {
        if (*lv1_p >= 0) {
          full_offset += *(lv1_p++);
        } else {
          full_offset = globals.cx1.lv1_items_special_[-1 - *(lv1_p++)];
        }

        int64_t read_id = globals.package.GetSeqID(full_offset >> 3);
        int offset = (full_offset >> 3) - globals.package.StartPos(read_id);
        int strand = full_offset & 1;
        int edge_type = (full_offset >> 1) & 3;
        int read_length = globals.package.SequenceLength(read_id);

        auto ptr_and_offset = globals.package.WordPtrAndOffset(read_id);
        int64_t start_offset = ptr_and_offset.second;
        const uint32_t *read_p = ptr_and_offset.first;
        int words_this_read = DivCeiling(start_offset + read_length, 16);

        if (strand == 0) {
          int num_chars_to_copy = globals.kmer_k;
          uint8_t prev = kSentinelValue;

          switch (edge_type) {
            case 0:
              break;

            case 1:
              prev = globals.package.GetBase(read_id, offset);
              offset++;
              break;

            case 2:
              prev = globals.package.GetBase(read_id, offset + 1);
              offset += 2;
              num_chars_to_copy--;
              break;

            default:
              assert(false);
          }

          CopySubstring(substr, read_p, offset + start_offset, num_chars_to_copy, 1, words_this_read,
                        globals.words_per_substring);

          uint32_t *last_word = substr + (globals.words_per_substring - 1) * 1;
          *last_word |= int(num_chars_to_copy == globals.kmer_k) << kBWTCharNumBits;
          *last_word |= prev;
        } else {
          int num_chars_to_copy = globals.kmer_k;
          uint8_t prev = kSentinelValue;

          switch (edge_type) {
            case 0:
              num_chars_to_copy--;
              prev = 3 - globals.package.GetBase(read_id, offset + globals.kmer_k - 1);
              break;

            case 1:
              prev = 3 - globals.package.GetBase(read_id, offset + globals.kmer_k);
              break;

            case 2:
              offset++;
              break;

            default:
              assert(false);
          }

          CopySubstringRC(substr, read_p, offset + start_offset, num_chars_to_copy, 1, words_this_read,
                          globals.words_per_substring);

          uint32_t *last_word = substr + (globals.words_per_substring - 1) * 1;
          *last_word |= int(num_chars_to_copy == globals.kmer_k) << kBWTCharNumBits;
          *last_word |= prev;
        }

        substr += globals.words_per_substring;
      }
    }
  }
}

void output_(int64_t from, int64_t to, read2sdbg_global_t &globals, uint32_t *substr, int tid) {
  int64_t start_idx, end_idx;
  int has_solid_a = 0;  // has solid (k+1)-mer aSb
  int has_solid_b = 0;  // has solid aSb
  int64_t last_a[4], outputed_b;
  uint32_t tip_label[32];
  SdbgWriter::Snapshot snapshot;

  for (start_idx = from; start_idx < to; start_idx = end_idx) {
    end_idx = start_idx + 1;
    uint32_t *item = substr + start_idx * globals.words_per_substring;

    while (end_idx < to &&
           !IsDiffKMinusOneMer(item, substr + end_idx * globals.words_per_substring, 1, globals.kmer_k)) {
      ++end_idx;
    }

    // clean marking
    has_solid_a = has_solid_b = 0;
    outputed_b = 0;

    for (int64_t i = start_idx; i < end_idx; ++i) {
      uint32_t *cur_item = substr + i * globals.words_per_substring;
      int a = Extract_a(cur_item, globals.words_per_substring, 1, globals.kmer_k);
      int b = Extract_b(cur_item, globals.words_per_substring, 1);

      if (a != kSentinelValue && b != kSentinelValue) {
        has_solid_a |= 1 << a;
        has_solid_b |= 1 << b;
      }

      if (a != kSentinelValue && (b != kSentinelValue || !(has_solid_a & (1 << a)))) {
        last_a[a] = i;
      }
    }

    for (int64_t i = start_idx, j; i < end_idx; i = j) {
      uint32_t *cur_item = substr + i * globals.words_per_substring;
      int a = Extract_a(cur_item, globals.words_per_substring, 1, globals.kmer_k);
      int b = Extract_b(cur_item, globals.words_per_substring, 1);

      j = i + 1;

      while (j < end_idx) {
        uint32_t *next_item = substr + j * globals.words_per_substring;

        if (Extract_a(next_item, globals.words_per_substring, 1, globals.kmer_k) != a ||
            Extract_b(next_item, globals.words_per_substring, 1) != b) {
          break;
        } else {
          ++j;
        }
      }

      int w, last, is_dollar = 0;
      int64_t count = std::min(j - i, int64_t(kMaxMul));

      if (a == kSentinelValue) {
        assert(b != kSentinelValue);

        if (has_solid_b & (1 << b)) {
          continue;
        }

        is_dollar = 1;
      }

      if (b == kSentinelValue) {
        assert(a != kSentinelValue);

        if (has_solid_a & (1 << a)) {
          continue;
        }
      }

      w = (b == kSentinelValue) ? 0 : ((outputed_b & (1 << b)) ? b + 5 : b + 1);
      last = (a == kSentinelValue) ? 0 : ((last_a[a] == j - 1) ? 1 : 0);
      outputed_b |= 1 << b;

      if (is_dollar) {
        for (int64_t i = 0; i < globals.words_per_dummy_node; ++i) {
          tip_label[i] = cur_item[i];
        }
      }

      globals.sdbg_writer.Write(tid, cur_item[0] >> (32 - kBucketPrefixLength * 2), w, last, is_dollar, count,
                                tip_label, &snapshot);
    }
  }
  globals.sdbg_writer.SaveSnapshot(snapshot);
}

struct kt_sort_t {
  read2sdbg_global_t *globals;
  std::vector<int64_t> thread_offset;
  std::vector<int> rank;
  int64_t acc = 0;
  int seen = 0;
  std::mutex mutex;
};

void kt_sort(void *g, long b, int tid) {
  auto kg = reinterpret_cast<kt_sort_t *>(g);
  if (kg->thread_offset[tid] == -1) {
    std::lock_guard<std::mutex> lk(kg->mutex);
    kg->thread_offset[tid] = kg->acc;
    kg->acc += kg->globals->cx1.bucket_sizes_[b];
    kg->rank[tid] = kg->seen;
    kg->seen++;
  }

  if (kg->globals->cx1.bucket_sizes_[b] == 0) {
    return;
  }

  size_t offset = kg->globals->cx1.lv1_num_items_ +
                  kg->thread_offset[tid] * kg->globals->cx1.bytes_per_sorting_item_ / sizeof(uint32_t);
  auto substr_ptr = reinterpret_cast<uint32_t *>(kg->globals->lv1_items.data() + offset);
  s2_lv2_extract_substr_(b, b + 1, *(kg->globals), substr_ptr);
  SortSubStr(substr_ptr, kg->globals->words_per_substring, kg->globals->cx1.bucket_sizes_[b]);
  output_(0, kg->globals->cx1.bucket_sizes_[b], *(kg->globals), substr_ptr, tid);
}

void s2_lv1_direct_sort_and_proc(read2sdbg_global_t &globals) {
  kt_sort_t kg;
  kg.globals = &globals;

  kg.thread_offset.resize(globals.num_cpu_threads, -1);
  kg.rank.resize(globals.num_cpu_threads, 0);
  omp_set_num_threads(globals.num_cpu_threads);
#pragma omp parallel for schedule(dynamic)
  for (auto i = globals.cx1.lv1_start_bucket_; i < globals.cx1.lv1_end_bucket_; ++i) {
    kt_sort(&kg, i, omp_get_thread_num());
  }
}

void s2_post_proc(read2sdbg_global_t &globals) {
  globals.sdbg_writer.Finalize();
  xinfo("Number of $ A C G T A- C- G- T-:\n");
  xinfo("");
  for (int i = 0; i < 9; ++i) {
    xinfoc("%lld ", (long long)globals.sdbg_writer.final_meta().w_count(i));
  }
  xinfoc("\n");
  xinfo("Total number of edges: %lld\n", (long long)globals.sdbg_writer.final_meta().item_count());
  xinfo("Total number of ONEs: %lld\n", (long long)globals.sdbg_writer.final_meta().ones_in_last());
  xinfo("Total number of $v edges: %lld\n", (long long)globals.sdbg_writer.final_meta().tip_count());
  assert(globals.sdbg_writer.final_meta().w_count(0) == globals.sdbg_writer.final_meta().tip_count());
}

}  // namespace s2

}  // namespace cx1_read2sdbg