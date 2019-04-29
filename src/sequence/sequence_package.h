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

#ifndef MEGAHIT_SEQUENCE_PACKAGE_H
#define MEGAHIT_SEQUENCE_PACKAGE_H

#include <cstdint>
#include <cassert>
#include <vector>
#include "kmlib/kmbit.h"
#include "kmlib/kmcompactvector.h"
#include "utils/utils.h"

/**
 * @brief hold a set of sequences
 */
template<class WordType = unsigned long>
class SequencePackage {
 public:
  using word_type = WordType;
  using vector_type = kmlib::CompactVector<2, word_type, kmlib::kBigEndian>;
  const static unsigned kCharsPerWord = vector_type::kBasesPerWord;

 public:
  SequencePackage() {
    clear();
    for (int i = 0; i < 10; ++i) {
      dna_map_[static_cast<int>("ACGTNacgtn"[i])] = "0123201232"[i] - '0';
    }
  }

  void clear() {
    sequences_.clear();
    start_pos_.clear();
    start_pos_.push_back(0);
    pos_to_id_.clear();
    max_len_ = 0;
    fixed_len_ = 0;
    num_fixed_len_ = 0;
  }

  void ReserveBases(size_t num_bases) {
    sequences_.reserve(num_bases);
  }

  void ReserveSequences(size_t num_seq) {
    start_pos_.reserve(num_seq + 1);
  }

  size_t size() const {
    return num_fixed_len_ + start_pos_.size() - 1;
  }

  size_t BaseCount() const {
    return start_pos_.back();
  }

  size_t SizeInByte() const {
    return sizeof(word_type) * sequences_.word_capacity()
        + sizeof(uint64_t) * start_pos_.capacity()
        + sizeof(uint64_t) * pos_to_id_.capacity();
  }

  unsigned MaxSequenceLength() const {
    return max_len_;
  }

  unsigned SequenceLength(size_t seq_id) const {
    if (seq_id < num_fixed_len_) {
      return fixed_len_;
    } else {
      return start_pos_[seq_id - num_fixed_len_ + 1] - start_pos_[seq_id - num_fixed_len_];
    }
  }

  uint8_t GetBase(size_t seq_id, unsigned offset) const {
    return sequences_[StartPos(seq_id) + offset];
  }

  uint64_t StartPos(size_t seq_id) const {
    if (seq_id < num_fixed_len_) {
      return seq_id * fixed_len_;
    } else {
      return start_pos_[seq_id - num_fixed_len_];
    }
  }

  std::pair<const word_type *, unsigned> WordPtrAndOffset(size_t seq_id, unsigned offset = 0) const {
    size_t index = StartPos(seq_id) + offset;
    return {sequences_.data() + index / kCharsPerWord, index % kCharsPerWord};
  }

  void BuildIndex() {
    pos_to_id_.clear();
    pos_to_id_.reserve(start_pos_.back() / kLookupStep + 4);
    size_t abs_offset = num_fixed_len_ * fixed_len_;
    size_t cur_id = num_fixed_len_;

    while (abs_offset <= start_pos_.back()) {
      while (cur_id < size() && start_pos_[cur_id - num_fixed_len_ + 1] <= abs_offset) {
        ++cur_id;
      }

      pos_to_id_.push_back(cur_id);
      abs_offset += kLookupStep;
    }

    pos_to_id_.push_back(size());
    pos_to_id_.push_back(size());
  }

  uint64_t GetSeqID(size_t full_offset) {
    if (full_offset < num_fixed_len_ * fixed_len_) {
      return full_offset / fixed_len_;
    } else {
      size_t look_up_entry = (full_offset - num_fixed_len_ * fixed_len_) / kLookupStep;
      size_t l = pos_to_id_[look_up_entry], r = pos_to_id_[look_up_entry + 1];

      while (l < r) {
        size_t mid = (l + r) / 2;
        if (start_pos_[mid - num_fixed_len_] > full_offset) {
          r = mid - 1;
        } else if (start_pos_[mid - num_fixed_len_ + 1] <= full_offset) {
          l = mid + 1;
        } else {
          return mid;
        }
      }
      return l;
    }
  }

 public:
  void AppendStringSequence(const char *s, unsigned len) {
    AppendStringSequence(s, s + len, len);
  }

  void AppendReversedStringSequence(const char *s, unsigned len) {
    AppendStringSequence(s + len - 1, s - 1, len);
  }

  void AppendCompactSequence(const word_type *s, unsigned len) {
    AppendCompactSequence(s, len, false);
  }

  void AppendReversedCompactSequence(const word_type *s, unsigned len) {
    AppendCompactSequence(s, len, true);
  }

  void FetchSequence(size_t seq_id, std::vector<word_type> *s) const {
    vector_type cvec(s);
    auto ptr_and_offset = WordPtrAndOffset(seq_id);
    auto ptr = ptr_and_offset.first;
    auto offset = ptr_and_offset.second;
    auto len = SequenceLength(seq_id);
    if (offset != 0) {
      unsigned remaining_len = std::min(len, kCharsPerWord - offset);
      cvec.push_word(*ptr, offset, remaining_len);
      len -= remaining_len;
      ++ptr;
    }
    unsigned n_full = len / kCharsPerWord;
    for (unsigned i = 0; i < n_full; ++i) {
      cvec.push_word(ptr[i]);
    }
    if (len % kCharsPerWord > 0) {
      cvec.push_word(ptr[n_full], 0, len % kCharsPerWord);
    }
  }

 private:
  bool IsFixedLength() const {
    return start_pos_.size() == 1;
  }
  void UpdateLength(unsigned len) {
    if (num_fixed_len_ == 0) {
      num_fixed_len_ = 1;
      fixed_len_ = len;
      start_pos_.back() += len;
    } else if (IsFixedLength() && len == fixed_len_) {
      num_fixed_len_++;
      start_pos_.back() += len;
    } else {
      start_pos_.push_back(start_pos_.back() + len);
    }
    if (len > max_len_) {
      max_len_ = len;
    }
  }

  void AppendStringSequence(const char *from, const char *to, unsigned len) {
    UpdateLength(len);
    std::ptrdiff_t step = from < to ? 1 : -1;
    for (auto ptr = from; ptr != to; ptr += step) {
      sequences_.push_back(dna_map_[static_cast<int>(*ptr)]);
    }
  }

  void AppendCompactSequence(const word_type *ptr, unsigned len, bool rev) {
    UpdateLength(len);
    if (rev) {
      auto rptr = ptr + DivCeiling(len, kCharsPerWord) - 1;
      unsigned bases_in_last_word = len % kCharsPerWord;
      if (bases_in_last_word > 0) {
        auto val = kmlib::bit::Reverse<2>(*rptr);
        sequences_.push_word(val, kCharsPerWord - bases_in_last_word, bases_in_last_word);
        --rptr;
      }
      for (auto p = rptr; p >= ptr; --p) {
        sequences_.push_word(kmlib::bit::Reverse<2>(*p));
      }
    } else {
      while (len >= kCharsPerWord) {
        sequences_.push_word(*ptr);
        len -= kCharsPerWord;
        ++ptr;
      }
      if (len > 0) {
        sequences_.push_word(*ptr, 0, len);
      }
    }
  }

 private:
  vector_type sequences_;
  unsigned fixed_len_{0};
  size_t num_fixed_len_{0};
  std::vector<uint64_t> start_pos_; // the index of the starting position of a sequence
  char dna_map_[256]{};
  unsigned max_len_{0};

  // for looking up the seq_id of a full offset
  std::vector<uint64_t> pos_to_id_;
  const static unsigned kLookupStep = 1024;
};

using SeqPackage = SequencePackage<uint32_t>;

#endif
