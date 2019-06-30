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

#include <cassert>
#include <cstdint>
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
  /**
   * The sequence view of a sequence
   */
  class SeqView {
   public:
    SeqView(const SequencePackage *pkg, int64_t seq_id)
        : package_(pkg), seq_id_(seq_id) {}

    unsigned length() const {
      return package_->GetSeqLength(seq_id_);
    }

    std::pair<const WordType*, unsigned> raw_address(unsigned offset = 0) const {
      return package_->GetRawAddress(seq_id_, offset);
    }

    size_t base_at(unsigned index) const {
      return package_->GetBase(seq_id_, index);
    }

    int64_t id() const { return seq_id_; }

    int64_t full_offset_in_pkg() const {
      return package_->StartPos(seq_id_);
    }
   private:
    const SequencePackage *package_;
    int64_t seq_id_;
  };

  using TWord = WordType;
  using TVector = kmlib::CompactVector<2, TWord, kmlib::kBigEndian>;
  using TView = SeqView;
  using TAddress = std::pair<const TWord*, unsigned>;
  const static unsigned kBasesPerWord = TVector::kBasesPerWord;

 public:
  SequencePackage() {
    Clear();
    for (int i = 0; i < 10; ++i) {
      dna_map_[static_cast<int>("ACGTNacgtn"[i])] = "0123201232"[i] - '0';
    }
  }

  void Clear() {
    sequences_.clear();
    start_pos_.clear();
    start_pos_.push_back(0);
    pos_to_id_.clear();
    max_len_ = 0;
    fixed_len_ = 0;
    num_fixed_len_ = 0;
  }

  void ReserveBases(size_t num_bases) { sequences_.reserve(num_bases); }

  void ReserveSequences(size_t num_seq) { start_pos_.reserve(num_seq + 1); }

  size_t seq_count() const { return num_fixed_len_ + start_pos_.size() - 1; }

  size_t base_count() const { return start_pos_.back(); }

  size_t size_in_byte() const {
    return sizeof(TWord) * sequences_.word_capacity() + sizeof(uint64_t) * start_pos_.capacity() +
        sizeof(uint64_t) * pos_to_id_.capacity();
  }

  unsigned max_length() const { return max_len_; }

  TView GetSeqView(int64_t seq_id) const {
    return TView(this, seq_id);
  }

  TView GetSeqViewByOffset(size_t offset) const {
    return TView(this, GetSeqID(offset));
  }

 private:
  unsigned GetSeqLength(size_t seq_id) const {
    if (seq_id < num_fixed_len_) {
      return fixed_len_;
    } else {
      return start_pos_[seq_id - num_fixed_len_ + 1] - start_pos_[seq_id - num_fixed_len_];
    }
  }

  uint8_t GetBase(size_t seq_id, unsigned offset) const { return sequences_[StartPos(seq_id) + offset]; }

  uint64_t StartPos(size_t seq_id) const {
    if (seq_id < num_fixed_len_) {
      return seq_id * fixed_len_;
    } else {
      return start_pos_[seq_id - num_fixed_len_];
    }
  }

  TAddress GetRawAddress(size_t seq_id, unsigned offset = 0) const {
    size_t index = StartPos(seq_id) + offset;
    return {sequences_.data() + index / kBasesPerWord, index % kBasesPerWord};
  }

  uint64_t GetSeqID(size_t full_offset) const {
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
  void AppendStringSequence(const char *s, unsigned len) { AppendStringSequence(s, s + len, len); }

  void AppendReversedStringSequence(const char *s, unsigned len) { AppendStringSequence(s + len - 1, s - 1, len); }

  void AppendCompactSequence(const TWord *s, unsigned len) { AppendCompactSequence(s, len, false); }

  void AppendReversedCompactSequence(const TWord *s, unsigned len) { AppendCompactSequence(s, len, true); }

  void FetchSequence(size_t seq_id, std::vector<TWord> *s) const {
    TVector cvec(s);
    auto ptr_and_offset = GetRawAddress(seq_id);
    auto ptr = ptr_and_offset.first;
    auto offset = ptr_and_offset.second;
    auto len = GetSeqLength(seq_id);
    if (offset != 0) {
      unsigned remaining_len = std::min(len, kBasesPerWord - offset);
      cvec.push_word(*ptr, offset, remaining_len);
      len -= remaining_len;
      ++ptr;
    }
    unsigned n_full = len / kBasesPerWord;
    for (unsigned i = 0; i < n_full; ++i) {
      cvec.push_word(ptr[i]);
    }
    if (len % kBasesPerWord > 0) {
      cvec.push_word(ptr[n_full], 0, len % kBasesPerWord);
    }
  }

  void BuildIndex() {
    pos_to_id_.clear();
    pos_to_id_.reserve(start_pos_.back() / kLookupStep + 4);
    size_t abs_offset = num_fixed_len_ * fixed_len_;
    size_t cur_id = num_fixed_len_;

    while (abs_offset <= start_pos_.back()) {
      while (cur_id < seq_count() && start_pos_[cur_id - num_fixed_len_ + 1] <= abs_offset) {
        ++cur_id;
      }

      pos_to_id_.push_back(cur_id);
      abs_offset += kLookupStep;
    }

    pos_to_id_.push_back(seq_count());
    pos_to_id_.push_back(seq_count());
  }

 private:
  bool IsFixedLength() const { return start_pos_.size() == 1; }

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

  void AppendCompactSequence(const TWord *ptr, unsigned len, bool rev) {
    UpdateLength(len);
    if (rev) {
      auto rptr = ptr + DivCeiling(len, kBasesPerWord) - 1;
      unsigned bases_in_last_word = len % kBasesPerWord;
      if (bases_in_last_word > 0) {
        auto val = kmlib::bit::Reverse<2>(*rptr);
        sequences_.push_word(val, kBasesPerWord - bases_in_last_word, bases_in_last_word);
        --rptr;
      }
      for (auto p = rptr; p >= ptr; --p) {
        sequences_.push_word(kmlib::bit::Reverse<2>(*p));
      }
    } else {
      while (len >= kBasesPerWord) {
        sequences_.push_word(*ptr);
        len -= kBasesPerWord;
        ++ptr;
      }
      if (len > 0) {
        sequences_.push_word(*ptr, 0, len);
      }
    }
  }

 private:
  TVector sequences_;
  unsigned fixed_len_{0};
  size_t num_fixed_len_{0};
  std::vector<uint64_t> start_pos_;  // the index of the starting position of a sequence
  char dna_map_[256]{};
  unsigned max_len_{0};

  // for looking up the seq_id of a full offset
  std::vector<uint64_t> pos_to_id_;
  const static unsigned kLookupStep = 1024;
};

using SeqPackage = SequencePackage<uint32_t>;

#endif
