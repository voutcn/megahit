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

#include <algorithm>

#include "sequence/copy_substr.h"
#include "sequence/kmer.h"
#include "utils/utils.h"

namespace {  // helper functions

/**
 * @brief encode read_id and its offset in one int64_t
 */
inline int64_t EncodeOffset(int64_t read_id, int offset, int strand,
                            const SeqPackage &p) {
  return ((p.GetSeqView(read_id).full_offset_in_pkg() + offset) << 1) | strand;
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

inline uint64_t ComposeUint64(uint32_t *read_info_ptr) {
  return (static_cast<uint64_t>(*read_info_ptr) << 32u) | *(read_info_ptr + 1);
}

inline uint8_t ExtractHeadTail(const uint32_t *item, int64_t spacing,
                               int words_per_substring) {
  return *(item + spacing * (words_per_substring - 1)) &
         ((1 << 2 * Read2SdbgS1::kBWTCharNumBits) - 1);
}

inline uint8_t ExtractPrevNext(int64_t readinfo) {
  return readinfo & ((1 << 2 * Read2SdbgS1::kBWTCharNumBits) - 1);
}

}  // namespace

// sorting core functions

int64_t Read2SdbgS1::Lv0EncodeDiffBase(int64_t read_id) {
  return EncodeOffset(read_id, 0, 0, seq_pkg_->package);
}

Read2SdbgS1::MemoryStat Read2SdbgS1::Initialize() {
  bool is_reverse = true;
  int64_t num_bases, num_reads;

  SequenceLibCollection seq_collection(opt_.read_lib_file);
  auto collection_size = seq_collection.GetSize();
  num_bases = collection_size.first;
  num_reads = collection_size.second;

  seq_pkg_->package.ReserveSequences(num_reads);
  seq_pkg_->package.ReserveBases(num_bases);

  seq_collection.Read(&seq_pkg_->package, is_reverse);
  seq_pkg_->package.BuildIndex();

  xinfo("{} reads, {} max read length, {} total bases\n",
        seq_pkg_->package.seq_count(), seq_pkg_->package.max_length(),
        seq_pkg_->package.base_count());

  words_per_substr_ = DivCeiling(
      (opt_.k - 1) * kBitsPerEdgeChar + 2 * kBWTCharNumBits, kBitsPerEdgeWord);
  xinfo("{} words per substring\n", words_per_substr_);

  if (opt_.solid_threshold > 1) {
    seq_pkg_->is_solid.reset(seq_pkg_->package.base_count());
  }

  // --- initialize output mercy files ---
  seq_pkg_->n_mercy_files = 1;

  while (seq_pkg_->n_mercy_files * 10485760LL <
             static_cast<int64_t>(seq_pkg_->package.seq_count()) &&
         seq_pkg_->n_mercy_files < 64) {
    seq_pkg_->n_mercy_files <<= 1;
  }
  xinfo("Number of files for mercy candidate reads: {}\n",
        seq_pkg_->n_mercy_files);

  for (int i = 0; i < seq_pkg_->n_mercy_files; ++i) {
    auto file_name = opt_.output_prefix + ".mercy_cand." + std::to_string(i);
    mercy_files_.push_back(xfopen(file_name.c_str(), "wb"));
  }

  // --- initialize stat ---
  edge_counter_.SetNumThreads(opt_.n_threads);
  int64_t memory_for_data = DivCeiling(seq_pkg_->is_solid.size(), 8) +
                            seq_pkg_->package.size_in_byte() +
                            edge_counter_.size_in_byte();

  return {
      num_reads,
      memory_for_data,
      words_per_substr_ + 2,
      2,
  };
}

void Read2SdbgS1::Lv0CalcBucketSize(int64_t seq_from, int64_t seq_to,
                                    std::array<int64_t, kNumBuckets> *out) {
  auto &bucket_sizes = *out;
  std::fill(bucket_sizes.begin(), bucket_sizes.end(), 0);
  GenericKmer k_minus1_mer, rev_k_minus1_mer;  // (k-1)-mer and its rc

  for (int64_t read_id = seq_from; read_id < seq_to; ++read_id) {
    auto seq_view = seq_pkg_->package.GetSeqView(read_id);
    auto read_length = seq_view.length();

    if (read_length < opt_.k + 1) {
      continue;
    }

    auto ptr_and_offset = seq_view.raw_address();
    int64_t offset = ptr_and_offset.second;
    const uint32_t *read_p = ptr_and_offset.first;

    k_minus1_mer.InitFromPtr(read_p, offset, opt_.k - 1);
    rev_k_minus1_mer = k_minus1_mer;
    rev_k_minus1_mer.ReverseComplement(opt_.k - 1);

    // the first one special handling
    bucket_sizes[k_minus1_mer.data()[0] >>
                 (kCharsPerEdgeWord - kBucketPrefixLength) *
                     kBitsPerEdgeChar]++;
    bucket_sizes[rev_k_minus1_mer.data()[0] >>
                 (kCharsPerEdgeWord - kBucketPrefixLength) *
                     kBitsPerEdgeChar]++;

    unsigned last_char_offset = opt_.k - 1;
    int c = seq_view.base_at(last_char_offset);
    k_minus1_mer.ShiftAppend(c, opt_.k - 1);
    rev_k_minus1_mer.ShiftPreappend(3 - c, opt_.k - 1);

    while (last_char_offset < read_length - 1) {
      int cmp = k_minus1_mer.cmp(rev_k_minus1_mer, opt_.k - 1);

      if (cmp > 0) {
        bucket_sizes[rev_k_minus1_mer.data()[0] >>
                     (kCharsPerEdgeWord - kBucketPrefixLength) *
                         kBitsPerEdgeChar]++;
      } else {
        bucket_sizes[k_minus1_mer.data()[0] >>
                     (kCharsPerEdgeWord - kBucketPrefixLength) *
                         kBitsPerEdgeChar]++;
      }

      c = seq_view.base_at(++last_char_offset);
      k_minus1_mer.ShiftAppend(c, opt_.k - 1);
      rev_k_minus1_mer.ShiftPreappend(3 - c, opt_.k - 1);
    }

    // last one special handling
    bucket_sizes[k_minus1_mer.data()[0] >>
                 (kCharsPerEdgeWord - kBucketPrefixLength) *
                     kBitsPerEdgeChar]++;
    bucket_sizes[rev_k_minus1_mer.data()[0] >>
                 (kCharsPerEdgeWord - kBucketPrefixLength) *
                     kBitsPerEdgeChar]++;
  }
}

void Read2SdbgS1::Lv1FillOffsets(OffsetFiller &filler, int64_t seq_from,
                                 int64_t seq_to) {
  GenericKmer k_minus1_mer, rev_k_minus1_mer;  // (k+1)-mer and its rc
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

    k_minus1_mer.InitFromPtr(read_p, offset, opt_.k - 1);
    rev_k_minus1_mer = k_minus1_mer;
    rev_k_minus1_mer.ReverseComplement(opt_.k - 1);

// ===== this is a macro to save some copy&paste ================
#define CHECK_AND_SAVE_OFFSET(offset, strand)                             \
  do {                                                                    \
    if (filler.IsHandling(key)) {                                         \
      filler.WriteNextOffset(                                             \
          key, EncodeOffset(read_id, offset, strand, seq_pkg_->package)); \
    }                                                                     \
  } while (0)
    // =========== end macro ==========================

    // the first one special handling
    key = k_minus1_mer.data()[0] >>
          (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
    CHECK_AND_SAVE_OFFSET(0, 0);
    key = rev_k_minus1_mer.data()[0] >>
          (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
    CHECK_AND_SAVE_OFFSET(0, 1);

    unsigned last_char_offset = opt_.k - 1;
    int c = seq_view.base_at(last_char_offset);
    k_minus1_mer.ShiftAppend(c, opt_.k - 1);
    rev_k_minus1_mer.ShiftPreappend(3 - c, opt_.k - 1);

    // shift the key char by char
    while (last_char_offset < read_length - 1) {
      int cmp = k_minus1_mer.cmp(rev_k_minus1_mer, opt_.k - 1);

      if (cmp > 0) {
        key = rev_k_minus1_mer.data()[0] >>
              (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
        CHECK_AND_SAVE_OFFSET(last_char_offset - opt_.k + 2, 1);
      } else if (cmp < 0) {
        key = k_minus1_mer.data()[0] >>
              (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
        CHECK_AND_SAVE_OFFSET(last_char_offset - opt_.k + 2, 0);
      } else {
        // a not-that-math-correct solution if the edge is palindrome, but works
        // well enough
        int prev = seq_view.base_at(last_char_offset - (opt_.k - 1));
        int next = seq_view.base_at(last_char_offset + 1);

        if (prev <= 3 - next) {
          key = rev_k_minus1_mer.data()[0] >>
                (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
          CHECK_AND_SAVE_OFFSET(last_char_offset - opt_.k + 2, 0);
        } else {
          key = rev_k_minus1_mer.data()[0] >>
                (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
          CHECK_AND_SAVE_OFFSET(last_char_offset - opt_.k + 2, 1);
        }
      }

      c = seq_view.base_at(++last_char_offset);
      k_minus1_mer.ShiftAppend(c, opt_.k - 1);
      rev_k_minus1_mer.ShiftPreappend(3 - c, opt_.k - 1);
    }

    // the last one special handling
    key = k_minus1_mer.data()[0] >>
          (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
    CHECK_AND_SAVE_OFFSET(last_char_offset - opt_.k + 2, 0);
    key = rev_k_minus1_mer.data()[0] >>
          (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
    CHECK_AND_SAVE_OFFSET(last_char_offset - opt_.k + 2, 1);
  }

#undef CHECK_AND_SAVE_OFFSET
}

void Read2SdbgS1::Lv2ExtractSubString(OffsetFetcher &fetcher,
                                      SubstrPtr substr_ptr) {
  while (fetcher.HasNext()) {
    int64_t full_offset = fetcher.Next();
    auto seq_view = seq_pkg_->package.GetSeqViewByOffset(full_offset >> 1);
    unsigned strand = full_offset & 1;
    unsigned offset = (full_offset >> 1) - seq_view.full_offset_in_pkg();
    unsigned read_length = seq_view.length();
    unsigned num_chars_to_copy = opt_.k - 1;
    unsigned char prev, next, head,
        tail;  // (k+1)=abScd, prev=a, head=b, tail=c, next=d

    assert(offset < read_length);

    if (offset > 1) {
      head = seq_view.base_at(offset - 1);
      prev = seq_view.base_at(offset - 2);
    } else {
      prev = kSentinelValue;

      if (offset > 0) {
        head = seq_view.base_at(offset - 1);
      } else {
        head = kSentinelValue;
      }
    }

    if (offset + opt_.k < read_length) {
      tail = seq_view.base_at(offset + opt_.k - 1);
      next = seq_view.base_at(offset + opt_.k);
    } else {
      next = kSentinelValue;

      if (offset + opt_.k - 1 < read_length) {
        tail = seq_view.base_at(offset + opt_.k - 1);
      } else {
        tail = kSentinelValue;
      }
    }

    auto ptr_and_offset = seq_view.raw_address();
    int start_offset = ptr_and_offset.second;
    const uint32_t *read_p = ptr_and_offset.first;
    int words_this_read = DivCeiling(start_offset + read_length, 16);
    uint64_t read_info;

    if (strand == 0) {
      CopySubstring(substr_ptr, read_p, offset + start_offset,
                    num_chars_to_copy, 1, words_this_read, words_per_substr_);
      auto last_word = substr_ptr + int64_t(words_per_substr_ - 1) * 1;
      *last_word |= (head << kBWTCharNumBits) | tail;
      read_info = (full_offset << 6) | (prev << 3) | next;
    } else {
      CopySubstringRC(substr_ptr, read_p, offset + start_offset,
                      num_chars_to_copy, 1, words_this_read, words_per_substr_);
      auto last_word = substr_ptr + int64_t(words_per_substr_ - 1) * 1;
      *last_word |= ((tail == kSentinelValue ? kSentinelValue : 3 - tail)
                     << kBWTCharNumBits) |
                    (head == kSentinelValue ? kSentinelValue : 3 - head);
      read_info =
          (full_offset << 6) |
          ((next == kSentinelValue ? kSentinelValue : (3 - next)) << 3) |
          (prev == kSentinelValue ? kSentinelValue : (3 - prev));
    }

    DecomposeUint64(substr_ptr + words_per_substr_, read_info);
    substr_ptr += words_per_substr_ + 2;
  }
}

void Read2SdbgS1::Lv2Postprocess(int64_t from, int64_t to, int thread,
                                 uint32_t *substr) {
  int64_t end_idx;
  int64_t count_prev_head[5][5];
  int64_t count_tail_next[5][5];
  int64_t count_head_tail[(1 << 2 * kBWTCharNumBits) - 1];

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
      if (IsDiffKMinusOneMer(first_item,
                             substr + end_idx * (2 + words_per_substr_), 1,
                             opt_.k)) {
        break;
      }

      auto read_info = ComposeUint64(first_item + words_per_substr_);
      uint8_t prev_and_next = ExtractPrevNext(read_info);
      uint8_t head_and_tail = ExtractHeadTail(
          substr + end_idx * (words_per_substr_ + 2), 1, words_per_substr_);
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
        if (count_head_tail[(j << kBWTCharNumBits) | x] >=
            opt_.solid_threshold) {
          l_has_out |= 1 << j;
          r_has_in |= 1 << x;
        }
      }
    }

    while (i < end_idx) {
      uint8_t head_and_tail = ExtractHeadTail(
          substr + i * (2 + words_per_substr_), 1, words_per_substr_);
      uint8_t head = head_and_tail >> 3;
      uint8_t tail = head_and_tail & 7;

      if (head != kSentinelValue && tail != kSentinelValue) {
        edge_counter_.Add(count_head_tail[head_and_tail], thread);
      }
      if (head != kSentinelValue && tail != kSentinelValue &&
          count_head_tail[head_and_tail] >= opt_.solid_threshold) {
        for (int64_t j = 0; j < count_head_tail[head_and_tail]; ++j, ++i) {
          int64_t read_info =
              ComposeUint64(substr + i * (2 + words_per_substr_) +
                            words_per_substr_) >>
              6;
          auto seq_view = seq_pkg_->package.GetSeqViewByOffset(read_info >> 1);
          int strand = read_info & 1;
          int64_t read_id = seq_view.id();
          int offset = (read_info >> 1) - seq_view.full_offset_in_pkg() - 1;
          int l_offset = strand == 0 ? offset : offset + 1;
          int r_offset = strand == 0 ? offset + 1 : offset;

          // mark this is a solid edge
          seq_pkg_->is_solid.set((read_info >> 1) - 1);

          if (!(has_in & (1 << head))) {
            // no in
            int64_t packed_mercy_cand =
                ((seq_view.full_offset_in_pkg() + l_offset) << 2) |
                (1 + strand);
            fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1,
                   mercy_files_[read_id & (seq_pkg_->n_mercy_files - 1)]);
          }

          if (!(has_out & (1 << tail))) {
            // no out
            int64_t packed_mercy_cand =
                ((seq_view.full_offset_in_pkg() + r_offset) << 2) |
                (2 - strand);
            fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1,
                   mercy_files_[read_id & (seq_pkg_->n_mercy_files - 1)]);
          }
        }
      } else {
        // not solid, but we still need to tell whether its left/right kmer is
        // solid
        for (int64_t j = 0; j < count_head_tail[head_and_tail]; ++j, ++i) {
          int64_t read_info =
              ComposeUint64(substr + i * (2 + words_per_substr_) +
                            words_per_substr_) >>
              6;
          auto seq_view = seq_pkg_->package.GetSeqViewByOffset(read_info >> 1);
          int strand = read_info & 1;
          int64_t read_id = seq_view.id();
          int offset = (read_info >> 1) - seq_view.full_offset_in_pkg() - 1;

          int l_offset = strand == 0 ? offset : offset + 1;
          int r_offset = strand == 0 ? offset + 1 : offset;

          if (l_has_out & (1 << head)) {
            if (has_in & (1 << head)) {
              // has both in & out
              int64_t packed_mercy_cand =
                  ((seq_view.full_offset_in_pkg() + l_offset) << 2) | 0;
              fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1,
                     mercy_files_[read_id & (seq_pkg_->n_mercy_files - 1)]);
            } else {
              // has out but no in
              int64_t packed_mercy_cand =
                  ((seq_view.full_offset_in_pkg() + l_offset) << 2) |
                  (1 + strand);
              fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1,
                     mercy_files_[read_id & (seq_pkg_->n_mercy_files - 1)]);
            }
          } else {
            if (has_in & (1 << head)) {
              // has in but no out
              int64_t packed_mercy_cand =
                  ((seq_view.full_offset_in_pkg() + l_offset) << 2) |
                  (2 - strand);
              fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1,
                     mercy_files_[read_id & (seq_pkg_->n_mercy_files - 1)]);
            }
          }

          if (r_has_in & (1 << tail)) {
            if (has_out & (1 << tail)) {
              // has both in & out
              int64_t packed_mercy_cand =
                  ((seq_view.full_offset_in_pkg() + r_offset) << 2) | 0;
              fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1,
                     mercy_files_[read_id & (seq_pkg_->n_mercy_files - 1)]);
            } else {
              // has in but no out
              int64_t packed_mercy_cand =
                  ((seq_view.full_offset_in_pkg() + r_offset) << 2) |
                  (2 - strand);
              fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1,
                     mercy_files_[read_id & (seq_pkg_->n_mercy_files - 1)]);
            }
          } else {
            if (has_out & (1 << tail)) {
              // has out but no in
              int64_t packed_mercy_cand =
                  ((seq_view.full_offset_in_pkg() + r_offset) << 2) |
                  (1 + strand);
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
  // --- stat ---
  xinfo("Total number of solid edges: {}\n",
        edge_counter_.GetNumSolidEdges(opt_.solid_threshold));
  std::ofstream counting_file(std::string(opt_.output_prefix) + ".counting");
  edge_counter_.DumpStat(counting_file);
  for (int i = 0; i < seq_pkg_->n_mercy_files; ++i) {
    fclose(mercy_files_[i]);
  }
}
