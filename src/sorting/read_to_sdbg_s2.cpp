/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong & L3 Bioinformatics
 * Limited
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
#include <string>
#include <vector>

#include "kmlib/kmsort.h"
#include "sequence/copy_substr.h"
#include "sequence/kmer.h"
#include "utils/utils.h"

namespace {

// helper functions
inline int64_t EncodeOffset(int64_t read_id, int offset, int strand,
                            const SeqPackage &p, int edge_type) {
  // edge_type: 0 left $; 1 solid; 2 right $
  return ((p.GetSeqView(read_id).full_offset_in_pkg() + offset) << 3) |
         (edge_type << 1) | strand;
}

// helper: see whether two lv2 items have the same (k-1)-mer
inline bool IsDiffKMinusOneMer(uint32_t *item1, uint32_t *item2,
                               int64_t spacing, int kmer_k) {
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
inline int Extract_a(uint32_t *item, int num_words, int64_t spacing,
                     int kmer_k) {
  int non_dollar =
      (item[(num_words - 1) * spacing] >> Read2SdbgS2::kBWTCharNumBits) & 1;

  if (non_dollar) {
    int which_word = (kmer_k - 1) / kCharsPerEdgeWord;
    int word_index = (kmer_k - 1) % kCharsPerEdgeWord;
    return (item[which_word * spacing] >>
            (kCharsPerEdgeWord - 1 - word_index) * kBitsPerEdgeChar) &
           kEdgeCharMask;
  } else {
    return Read2SdbgS2::kSentinelValue;
  }
}

inline int Extract_b(uint32_t *item, int num_words, int64_t spacing) {
  return item[(num_words - 1) * spacing] &
         ((1 << Read2SdbgS2::kBWTCharNumBits) - 1);
}
}  // namespace

// sorting core functions
int64_t Read2SdbgS2::Lv0EncodeDiffBase(int64_t read_id) {
  return EncodeOffset(read_id, 0, 0, seq_pkg_->package, 0);
}

Read2SdbgS2::MemoryStat Read2SdbgS2::Initialize() {
  words_per_substr_ = DivCeiling(
      opt_.k * kBitsPerEdgeChar + kBWTCharNumBits + 1, kBitsPerEdgeWord);
  words_per_dummy_node_ =
      DivCeiling(opt_.k * kBitsPerEdgeChar, kBitsPerEdgeWord);
  xinfo("{} words per substring, words per dummy node ($v): {}\n",
        words_per_substr_, words_per_dummy_node_);

  // --- init output ---
  sdbg_writer_.set_num_threads(opt_.n_threads);
  sdbg_writer_.set_kmer_size(opt_.k);
  sdbg_writer_.set_num_buckets(kNumBuckets);
  sdbg_writer_.set_file_prefix(opt_.output_prefix);
  sdbg_writer_.InitFiles();

  int64_t memory_for_data = DivCeiling(seq_pkg_->is_solid.size(), 8) +
                            seq_pkg_->package.size_in_byte();

  MemoryStat ret{
      static_cast<int64_t>(seq_pkg_->package.seq_count()),
      memory_for_data,
      words_per_substr_,
      0,
  };

  if (!opt_.need_mercy || opt_.solid_threshold == 1) {
    return ret;
  }

  SimpleTimer timer;

  timer.reset();
  timer.start();
  xinfo("Adding mercy edges...\n");

  std::vector<uint64_t> mercy_cand;
  uint64_t num_mercy = 0;
  AtomicBitVector read_marker;
  read_marker.reset(seq_pkg_->package.seq_count());

  for (int fid = 0; fid < seq_pkg_->n_mercy_files; ++fid) {
    auto file_name = opt_.output_prefix + ".mercy_cand." + std::to_string(fid);
    std::ifstream mercy_file(file_name,
                             std::ifstream::binary | std::ifstream::in);
    mercy_cand.clear();

    uint64_t buf[4096];
    BufferedReader mercy_reader;
    mercy_reader.reset(&mercy_file);
    size_t size_read;

    while ((size_read = mercy_reader.read(buf, 4096)) > 0) {
      assert(size_read % sizeof(uint64_t) == 0);
      mercy_cand.insert(mercy_cand.end(), buf,
                        buf + size_read / sizeof(uint64_t));
    }

    xinfo("Mercy file: {s}, {}\n", file_name.c_str(), mercy_cand.size());
    kmlib::kmsort(mercy_cand.begin(), mercy_cand.end());

    // multi threading
    uint64_t avg = DivCeiling(mercy_cand.size(), opt_.n_threads);
    std::vector<uint64_t> start_idx(opt_.n_threads);
    std::vector<uint64_t> end_idx(opt_.n_threads);

    // manually distribute threads
    for (int tid = 0; tid < opt_.n_threads; ++tid) {
      if (tid == 0) {
        start_idx[tid] = 0;
      } else {
        start_idx[tid] = end_idx[tid - 1];
      }

      uint64_t this_end =
          std::min(start_idx[tid] + avg, (uint64_t)mercy_cand.size());
      size_t read_id = 0;
      if (this_end < mercy_cand.size()) {
        read_id =
            seq_pkg_->package.GetSeqViewByOffset(mercy_cand[this_end] >> 2)
                .id();
      }

      while (this_end < mercy_cand.size() &&
             seq_pkg_->package.GetSeqViewByOffset(mercy_cand[this_end] >> 2)
                     .id() == read_id) {
        ++this_end;
      }

      end_idx[tid] = this_end;
    }

    omp_set_num_threads(opt_.n_threads);

#pragma omp parallel for reduction(+ : num_mercy)
    for (int tid = 0; tid < opt_.n_threads; ++tid) {
      std::vector<bool> no_in(seq_pkg_->package.max_length());
      std::vector<bool> no_out(seq_pkg_->package.max_length());
      std::vector<bool> has_solid_kmer(seq_pkg_->package.max_length());

      uint64_t mercy_index = start_idx[tid];

      // go read by read
      while (mercy_index != end_idx[tid]) {
        auto seq_view =
            seq_pkg_->package.GetSeqViewByOffset(mercy_cand[mercy_index] >> 2);
        size_t read_id = seq_view.id();
        assert(!read_marker.at(read_id));
        read_marker.set(read_id);
        int first_0_out = seq_pkg_->package.max_length() + 1;
        int last_0_in = -1;

        std::fill(no_in.begin(), no_in.end(), false);
        std::fill(no_out.begin(), no_out.end(), false);
        std::fill(has_solid_kmer.begin(), has_solid_kmer.end(), false);

        while (
            mercy_index != end_idx[tid] &&
            seq_pkg_->package.GetSeqViewByOffset(mercy_cand[mercy_index] >> 2)
                    .id() == read_id) {
          int offset =
              (mercy_cand[mercy_index] >> 2) - seq_view.full_offset_in_pkg();
          if ((mercy_cand[mercy_index] & 3) == 2) {
            no_out[offset] = true;
            first_0_out = std::min(first_0_out, offset);
          } else if ((mercy_cand[mercy_index] & 3) == 1) {
            no_in[offset] = true;
            last_0_in = std::max(last_0_in, offset);
          }

          has_solid_kmer[offset] = true;
          ++mercy_index;
        }

        if (last_0_in < first_0_out) {
          continue;
        }

        auto read_length = seq_view.length();
        int last_no_out = -1;

        for (unsigned i = 0; i + opt_.k < read_length; ++i) {
          if (seq_pkg_->is_solid.at(seq_view.full_offset_in_pkg() + i)) {
            has_solid_kmer[i] = has_solid_kmer[i + 1] = true;
          }
        }

        for (unsigned i = 0; i + opt_.k <= read_length; ++i) {
          if (no_in[i] && last_no_out != -1) {
            for (unsigned j = last_no_out; j < i; ++j) {
              seq_pkg_->is_solid.set(seq_view.full_offset_in_pkg() + j);
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
  }

  timer.stop();
  xinfo("Adding mercy Done. Time elapsed: {.4}\n", timer.elapsed());
  xinfo("Number mercy: {}\n", num_mercy);

  return ret;
}

void Read2SdbgS2::Lv0CalcBucketSize(int64_t seq_from, int64_t seq_to,
                                    std::array<int64_t, kNumBuckets> *out) {
  auto &bucket_sizes = *out;
  std::fill(bucket_sizes.begin(), bucket_sizes.end(), 0);
  GenericKmer edge, rev_edge;  // (k+1)-mer and its rc

  for (int64_t read_id = seq_from; read_id < seq_to; ++read_id) {
    auto seq_view = seq_pkg_->package.GetSeqView(read_id);
    auto read_length = seq_view.length();

    if (read_length < opt_.k + 1) {
      continue;
    }

    auto ptr_and_offset = seq_view.raw_address();
    int64_t offset = ptr_and_offset.second;
    const uint32_t *read_p = ptr_and_offset.first;

    edge.InitFromPtr(read_p, offset, opt_.k + 1);
    rev_edge = edge;
    rev_edge.ReverseComplement(opt_.k + 1);

    unsigned last_char_offset = opt_.k;
    int64_t full_offset = seq_view.full_offset_in_pkg();
    bool for_sure_solid = opt_.solid_threshold == 1;

    while (true) {
      if (for_sure_solid || seq_pkg_->is_solid.at(full_offset)) {
        bool is_palindrome = rev_edge.cmp(edge, opt_.k + 1) == 0;
        bucket_sizes[(edge.data()[0] << 2) >>
                     (kCharsPerEdgeWord - kBucketPrefixLength) *
                         kBitsPerEdgeChar]++;

        if (!is_palindrome)
          bucket_sizes[(rev_edge.data()[0] << 2) >>
                       (kCharsPerEdgeWord - kBucketPrefixLength) *
                           kBitsPerEdgeChar]++;

        if (last_char_offset == opt_.k ||
            !(for_sure_solid || seq_pkg_->is_solid.at(full_offset - 1))) {
          bucket_sizes[edge.data()[0] >>
                       (kCharsPerEdgeWord - kBucketPrefixLength) *
                           kBitsPerEdgeChar]++;

          if (!is_palindrome)
            bucket_sizes[(rev_edge.data()[0] << 4) >>
                         (kCharsPerEdgeWord - kBucketPrefixLength) *
                             kBitsPerEdgeChar]++;
        }

        if (last_char_offset == read_length - 1 ||
            !(for_sure_solid || seq_pkg_->is_solid.at(full_offset + 1))) {
          bucket_sizes[(edge.data()[0] << 4) >>
                       (kCharsPerEdgeWord - kBucketPrefixLength) *
                           kBitsPerEdgeChar]++;

          if (!is_palindrome)
            bucket_sizes[rev_edge.data()[0] >>
                         (kCharsPerEdgeWord - kBucketPrefixLength) *
                             kBitsPerEdgeChar]++;
        }
      }

      ++full_offset;

      if (++last_char_offset >= read_length) {
        break;
      } else {
        int c = seq_view.base_at(last_char_offset);
        edge.ShiftAppend(c, opt_.k + 1);
        rev_edge.ShiftPreappend(3 - c, opt_.k + 1);
      }
    }
  }
}

void Read2SdbgS2::Lv1FillOffsets(OffsetFiller &filler, int64_t seq_from,
                                 int64_t seq_to) {
  GenericKmer edge, rev_edge;  // (k+1)-mer and its rc
  int key;

  for (int64_t read_id = seq_from; read_id < seq_to; ++read_id) {
    auto seq_view = seq_pkg_->package.GetSeqView(read_id);
    auto read_length = seq_view.length();

    if (read_length < opt_.k + 1) {
      continue;
    }

    auto ptr_and_offset = seq_view.raw_address();
    int64_t offset = ptr_and_offset.second;
    const uint32_t *read_p = ptr_and_offset.first;

    edge.InitFromPtr(read_p, offset, opt_.k + 1);
    rev_edge = edge;
    rev_edge.ReverseComplement(opt_.k + 1);

// ===== this is a macro to save some copy&paste ================
#define CHECK_AND_SAVE_OFFSET(offset, strand, edge_type)                       \
  do {                                                                         \
    if (filler.IsHandling(key)) {                                              \
      filler.WriteNextOffset(key, EncodeOffset(read_id, offset, strand,        \
                                               seq_pkg_->package, edge_type)); \
    }                                                                          \
  } while (0)
    // =========== end macro ==========================

    // shift the key char by char
    unsigned last_char_offset = opt_.k;
    int64_t full_offset = seq_view.full_offset_in_pkg();
    bool for_sure_solid = opt_.solid_threshold == 1;

    while (true) {
      if (for_sure_solid || seq_pkg_->is_solid.at(full_offset)) {
        bool is_palindrome = rev_edge.cmp(edge, opt_.k + 1) == 0;

        // left $
        if (last_char_offset == opt_.k ||
            !(for_sure_solid || seq_pkg_->is_solid.at(full_offset - 1))) {
          key = edge.data()[0] >>
                (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
          CHECK_AND_SAVE_OFFSET(last_char_offset - opt_.k, 0, 0);

          if (!is_palindrome) {
            key = (rev_edge.data()[0] << 4) >>
                  (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
            CHECK_AND_SAVE_OFFSET(last_char_offset - opt_.k, 1, 0);
          }
        }

        // solid
        key = (edge.data()[0] << 2) >>
              (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
        CHECK_AND_SAVE_OFFSET(last_char_offset - opt_.k, 0, 1);

        if (!is_palindrome) {
          key = (rev_edge.data()[0] << 2) >>
                (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
          CHECK_AND_SAVE_OFFSET(last_char_offset - opt_.k, 1, 1);
        }

        // right $
        if (last_char_offset == read_length - 1 ||
            !(for_sure_solid || seq_pkg_->is_solid.at(full_offset + 1))) {
          key = (edge.data()[0] << 4) >>
                (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
          CHECK_AND_SAVE_OFFSET(last_char_offset - opt_.k, 0, 2);

          if (!is_palindrome) {
            key = rev_edge.data()[0] >>
                  (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
            CHECK_AND_SAVE_OFFSET(last_char_offset - opt_.k, 1, 2);
          }
        }
      }

      ++full_offset;

      if (++last_char_offset >= read_length) {
        break;
      } else {
        int c = seq_view.base_at(last_char_offset);
        edge.ShiftAppend(c, opt_.k + 1);
        rev_edge.ShiftPreappend(3 - c, opt_.k + 1);
      }
    }
  }

#undef CHECK_AND_SAVE_OFFSET
}

void Read2SdbgS2::Lv2ExtractSubString(OffsetFetcher &fetcher,
                                      SubstrPtr substr) {
  while (fetcher.HasNext()) {
    int64_t full_offset = fetcher.Next();
    auto seq_view = seq_pkg_->package.GetSeqViewByOffset(full_offset >> 3);
    int offset = (full_offset >> 3) - seq_view.full_offset_in_pkg();
    int strand = full_offset & 1;
    int edge_type = (full_offset >> 1) & 3;
    int read_length = seq_view.length();

    auto ptr_and_offset = seq_view.raw_address();
    int64_t start_offset = ptr_and_offset.second;
    const uint32_t *read_p = ptr_and_offset.first;
    int words_this_read = DivCeiling(start_offset + read_length, 16);

    if (strand == 0) {
      unsigned num_chars_to_copy = opt_.k;
      uint8_t prev = kSentinelValue;

      switch (edge_type) {
        case 0:
          break;

        case 1:
          prev = seq_view.base_at(offset);
          offset++;
          break;

        case 2:
          prev = seq_view.base_at(offset + 1);
          offset += 2;
          num_chars_to_copy--;
          break;

        default:
          assert(false);
      }

      CopySubstring(substr, read_p, offset + start_offset, num_chars_to_copy, 1,
                    words_this_read, words_per_substr_);

      auto last_word = substr + (words_per_substr_ - 1) * 1;
      *last_word |= unsigned(num_chars_to_copy == opt_.k) << kBWTCharNumBits;
      *last_word |= prev;
    } else {
      unsigned num_chars_to_copy = opt_.k;
      uint8_t prev = kSentinelValue;

      switch (edge_type) {
        case 0: {
          num_chars_to_copy--;
          prev = 3 - seq_view.base_at(offset + opt_.k - 1);
          break;
        }
        case 1: {
          prev = 3 - seq_view.base_at(offset + opt_.k);
          break;
        }
        case 2: {
          offset++;
          break;
        }
        default: {
          assert(false);
        }
      }

      CopySubstringRC(substr, read_p, offset + start_offset, num_chars_to_copy,
                      1, words_this_read, words_per_substr_);

      auto last_word = substr + (words_per_substr_ - 1);
      *last_word |= unsigned(num_chars_to_copy == opt_.k) << kBWTCharNumBits;
      *last_word |= prev;
    }

    substr += words_per_substr_;
  }
}

void Read2SdbgS2::Lv2Postprocess(int64_t from, int64_t to, int tid,
                                 uint32_t *substr) {
  int64_t start_idx, end_idx;
  int has_solid_a = 0;  // has solid (k+1)-mer aSb
  int has_solid_b = 0;  // has solid aSb
  int64_t last_a[4], outputed_b;
  uint32_t tip_label[32];
  SdbgWriter::Snapshot snapshot;

  for (start_idx = from; start_idx < to; start_idx = end_idx) {
    end_idx = start_idx + 1;
    uint32_t *item = substr + start_idx * words_per_substr_;

    while (end_idx < to &&
           !IsDiffKMinusOneMer(item, substr + end_idx * words_per_substr_, 1,
                               opt_.k)) {
      ++end_idx;
    }

    // clean marking
    has_solid_a = has_solid_b = 0;
    outputed_b = 0;

    for (int64_t item_idx = start_idx; item_idx < end_idx; ++item_idx) {
      uint32_t *cur_item = substr + item_idx * words_per_substr_;
      int a = Extract_a(cur_item, words_per_substr_, 1, opt_.k);
      int b = Extract_b(cur_item, words_per_substr_, 1);

      if (a != kSentinelValue && b != kSentinelValue) {
        has_solid_a |= 1 << a;
        has_solid_b |= 1 << b;
      }

      if (a != kSentinelValue &&
          (b != kSentinelValue || !(has_solid_a & (1 << a)))) {
        last_a[a] = item_idx;
      }
    }

    for (int64_t item_idx = start_idx, j; item_idx < end_idx; item_idx = j) {
      uint32_t *cur_item = substr + item_idx * words_per_substr_;
      int a = Extract_a(cur_item, words_per_substr_, 1, opt_.k);
      int b = Extract_b(cur_item, words_per_substr_, 1);

      j = item_idx + 1;

      while (j < end_idx) {
        uint32_t *next_item = substr + j * words_per_substr_;

        if (Extract_a(next_item, words_per_substr_, 1, opt_.k) != a ||
            Extract_b(next_item, words_per_substr_, 1) != b) {
          break;
        } else {
          ++j;
        }
      }

      int w, last, is_dollar = 0;
      int64_t count = std::min(j - item_idx, int64_t(kMaxMul));

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
        for (int64_t i = 0; i < words_per_dummy_node_; ++i) {
          tip_label[i] = cur_item[i];
        }
      }

      sdbg_writer_.Write(tid, cur_item[0] >> (32 - kBucketPrefixLength * 2), w,
                         last, is_dollar, count, tip_label, &snapshot);
    }
  }
  sdbg_writer_.SaveSnapshot(snapshot);
}

void Read2SdbgS2::Lv0Postprocess() {
  sdbg_writer_.Finalize();
  xinfo("Number of $ A C G T A- C- G- T-:\n");
  xinfo("");
  for (int i = 0; i < 9; ++i) {
    xinfoc("{} ", sdbg_writer_.final_meta().w_count(i));
  }
  xinfoc("{}", "\n");
  xinfo("Total number of edges: {}\n", sdbg_writer_.final_meta().item_count());
  xinfo("Total number of ONEs: {}\n", sdbg_writer_.final_meta().ones_in_last());
  xinfo("Total number of $v edges: {}\n",
        sdbg_writer_.final_meta().tip_count());
  assert(sdbg_writer_.final_meta().w_count(0) ==
         sdbg_writer_.final_meta().tip_count());
}
