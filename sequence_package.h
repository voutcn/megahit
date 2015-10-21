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

#ifndef SEQUENCE_PACKAGE_H__
#define SEQUENCE_PACKAGE_H__

#include <stdint.h>
#include <assert.h>
#include <vector>

#include "bit_operation.h"

/**
 * @brief hold a set of sequences
 */

struct SequencePackage {

    typedef uint32_t word_t;	// do not change
    const static unsigned kBitsPerWord = 8 * sizeof(word_t);
    const static unsigned kCharsPerWord = kBitsPerWord / 2;
    char dna_map_[256];

    std::vector<word_t> packed_seq; // packed all
    std::vector<uint64_t> start_idx_; // the index of the starting position of a sequence

    uint8_t unused_bits_; // the number of unused bits in the last word
    int max_read_len_;

    // the package's first fixed_len_size_ items have the same length of fixed_len_
    // if fixed_len_sealed_, no more fixed-length items can be pushed
    int fixed_len_;
    size_t num_fixed_len_items_;
    bool fixed_len_sealed_;

    // for looking up the seq_id of a position
    std::vector<uint64_t> pos_to_id_;
    const static int kLookupStep = 1024;

    SequencePackage() {
        start_idx_.push_back(0);
        packed_seq.push_back(word_t(0));
        unused_bits_ = kBitsPerWord;
        max_read_len_ = 0;

        fixed_len_ = 0;
        num_fixed_len_items_ = 0;
        fixed_len_sealed_ = false;

        for (int i = 0; i < 10; ++i) {
            dna_map_[(int)("ACGTNacgtn"[i])] = "0123201232"[i] - '0';
        }
    }

    ~SequencePackage() {}

    void clear() {
        packed_seq.clear();
        packed_seq.push_back(word_t(0));
        start_idx_.clear();
        start_idx_.push_back(0);
        pos_to_id_.clear();
        unused_bits_ = kBitsPerWord;
        max_read_len_ = 0;

        fixed_len_ = 0;
        num_fixed_len_items_ = 0;
        fixed_len_sealed_ = false;
    }

    size_t size() {
        return num_fixed_len_items_ + start_idx_.size() - 1;
    }

    size_t base_size() {
        return start_idx_.back();
    }

    size_t size_in_byte() {
        return sizeof(word_t) * packed_seq.capacity() + sizeof(uint64_t) * start_idx_.capacity() + sizeof(uint64_t) * pos_to_id_.capacity();
    }

    size_t max_read_len() {
        return max_read_len_;
    }

    void shrink_to_fit() {
        // packed_seq.shrink_to_fit();
        // start_idx_.shrink_to_fit();
    }

    void reserve_bases(size_t num_bases) {
        packed_seq.reserve((num_bases + kCharsPerWord) / kCharsPerWord);
    }

    void reserve_num_seq(size_t num_seq) {
        start_idx_.reserve(num_seq + 1);
    }

    size_t length(size_t seq_id) {
        if (seq_id < num_fixed_len_items_) {
            return fixed_len_;
        }
        else {
            return start_idx_[seq_id - num_fixed_len_items_ + 1] - start_idx_[seq_id - num_fixed_len_items_];
        }
    }

    uint8_t get_base(size_t seq_id, size_t offset) {
        uint64_t where = get_start_index(seq_id) + offset;
        return packed_seq[where / kCharsPerWord] >> (kCharsPerWord - 1 - where % kCharsPerWord) * 2 & 3;
    }

    uint64_t get_start_index(size_t seq_id) {
        if (seq_id < num_fixed_len_items_) {
            return seq_id * fixed_len_;
        }
        else {
            return start_idx_[seq_id - num_fixed_len_items_];
        }
    }

    void set_fixed_len(int len) {
        assert(!fixed_len_sealed_);
        fixed_len_ = len;
    }

    void BuildLookup() {
        pos_to_id_.clear();
        pos_to_id_.reserve(start_idx_.back() / kLookupStep + 4);
        size_t abs_offset = num_fixed_len_items_ * fixed_len_;
        size_t cur_id = num_fixed_len_items_;

        while (abs_offset <= start_idx_.back()) {
            while (cur_id < size() && start_idx_[cur_id - num_fixed_len_items_ + 1] <= abs_offset) {
                ++cur_id;
            }

            pos_to_id_.push_back(cur_id);
            abs_offset += kLookupStep;
        }

        pos_to_id_.push_back(size());
        pos_to_id_.push_back(size());
    }

    uint64_t get_id(size_t abs_offset) {
        if (abs_offset < num_fixed_len_items_ * fixed_len_) {
            return abs_offset / fixed_len_;
        }
        else {
            size_t look_up_entry = (abs_offset - num_fixed_len_items_ * fixed_len_) / kLookupStep;
            size_t l = pos_to_id_[look_up_entry], r = pos_to_id_[look_up_entry + 1];

            while (l < r) {
                size_t mid = (l + r) / 2;

                if (start_idx_[mid - num_fixed_len_items_] > abs_offset) {
                    r = mid - 1;
                }
                else if (start_idx_[mid - num_fixed_len_items_ + 1] <= abs_offset) {
                    l = mid + 1;
                }
                else {
                    return mid;
                }
            }

            return l;
        }
    }

    void AppendFixedLenSeq(const char *s, int len) {
        assert(len == fixed_len_);
        assert(!fixed_len_sealed_);

        AddSeqToPackedSeq_(s, len);
        ++num_fixed_len_items_;
        start_idx_.back() += len;
    }

    void AppendFixedLenRevSeq(const char *s, int len) {
        assert(len == fixed_len_);
        assert(!fixed_len_sealed_);

        AddRevSeqToPackedSeq_(s, len);
        ++num_fixed_len_items_;
        start_idx_.back() += len;
    }

    void AppendFixedLenSeq(const word_t *s, int len) {
        assert(len == fixed_len_);
        assert(!fixed_len_sealed_);

        AddSeqToPackedSeq_(s, len);
        ++num_fixed_len_items_;
        start_idx_.back() += len;
    }

    void AppendFixedLenRevSeq(const word_t *s, int len) {
        assert(len == fixed_len_);
        assert(!fixed_len_sealed_);

        AddRevSeqToPackedSeq_(s, len);
        ++num_fixed_len_items_;
        start_idx_.back() += len;
    }

    void AppendSeq(const char *s, int len) {
        fixed_len_sealed_ = true;
        AddSeqToPackedSeq_(s, len);
        uint64_t end = start_idx_.back() + len;
        start_idx_.push_back(end);
    }

    void AppendReverseSeq(const char *s, int len) {
        fixed_len_sealed_ = true;
        AddRevSeqToPackedSeq_(s, len);
        uint64_t end = start_idx_.back() + len;
        start_idx_.push_back(end);
    }

    void AppendSeq(const word_t *s, int len) {
        fixed_len_sealed_ = true;
        AddSeqToPackedSeq_(s, len);
        uint64_t end = start_idx_.back() + len;
        start_idx_.push_back(end);
    }

    void AppendRevSeq(const word_t *s, int len) {
        fixed_len_sealed_ = true;
        AddRevSeqToPackedSeq_(s, len);
        uint64_t end = start_idx_.back() + len;
        start_idx_.push_back(end);
    }

    void AddSeqToPackedSeq_(const char *s, int len) {
        for (int i = 0; i < len; ++i) {
            unused_bits_ -= 2;
            packed_seq.back() |= dna_map_[(int)s[i]] << unused_bits_;

            if (unused_bits_ == 0) {
                unused_bits_ = kBitsPerWord;
                packed_seq.push_back(word_t(0));
            }
        }

        if (len > max_read_len_) {
            max_read_len_ = len;
        }
    }

    void AddRevSeqToPackedSeq_(const char *s, int len) {
        for (int i = len - 1; i >= 0; --i) {
            unused_bits_ -= 2;
            packed_seq.back() |= dna_map_[(int)s[i]] << unused_bits_;

            if (unused_bits_ == 0) {
                unused_bits_ = kBitsPerWord;
                packed_seq.push_back(word_t(0));
            }
        }

        if (len > max_read_len_) {
            max_read_len_ = len;
        }
    }

    void AddSeqToPackedSeq_(const word_t *s, int len) {
        if (len > max_read_len_) {
            max_read_len_ = len;
        }

        if (len == 0) {
            return;
        }

        if (len * 2 <= unused_bits_) {
            unused_bits_ -= len * 2;
            packed_seq.back() |= s[0] >> (kCharsPerWord - len) * 2 << unused_bits_;

            if (unused_bits_ == 0) {
                unused_bits_ = kBitsPerWord;
                packed_seq.push_back(word_t(0));
            }
        }
        else {
            int num_words = (len + kCharsPerWord - 1) / kCharsPerWord;

            // append to unused bits
            packed_seq.back() |= s[0] >> (kBitsPerWord - unused_bits_);

            if (unused_bits_ == kBitsPerWord) {
                packed_seq.back() = s[0];

                for (int i = 1; i < num_words; ++i) {
                    packed_seq.push_back(s[i]);
                }
            }
            else {
                int num_words_to_append = (len - unused_bits_ / 2 + kCharsPerWord - 1) / kCharsPerWord;

                for (int i = 0; i < num_words_to_append; ++i) {
                    if (i + 1 < num_words)
                        packed_seq.push_back((s[i] << unused_bits_) | (s[i + 1] >> (kBitsPerWord - unused_bits_)));
                    else
                        packed_seq.push_back(s[i] << unused_bits_);
                }
            }

            int bits_in_last_word = (start_idx_.back() + len) % kCharsPerWord * 2;
            unused_bits_ = kBitsPerWord - bits_in_last_word;

            if (unused_bits_ == kBitsPerWord) {
                packed_seq.push_back(word_t(0));
            }
            else {
                packed_seq.back() >>= unused_bits_;
                packed_seq.back() <<= unused_bits_;
            }
        }
    }

    void AddRevSeqToPackedSeq_(const word_t *rs, int len) {
        if (len == 0) {
            return;
        }

        std::vector<word_t> s(rs, rs + (len + kCharsPerWord - 1) / kCharsPerWord);

        for (int j = 0; j < (int)s.size(); ++j) {
            s[j] = bit_operation::Reverse(s[j]);
        }

        for (int j = 0, k = s.size() - 1; j < k; ++j, --k) {
            std::swap(s[j], s[k]);
        }

        unsigned shift = (kCharsPerWord - len % kCharsPerWord)  * 2;

        if (shift != kBitsPerWord) {
            for (int j = 0; j < (int)s.size() - 1; ++j) {
                s[j] = (s[j] << shift) | (s[j + 1] >> (kBitsPerWord - shift));
            }

            s.back() <<= shift;
        }

        AddSeqToPackedSeq_(&s[0], len);
    }

    void get_seq(std::vector<word_t> &s, size_t seq_id, int begin = 0, int end = -1) {
        if (end == -1) {
            end = length(seq_id) - 1;
        }

        size_t first_word = (get_start_index(seq_id) + begin) / kCharsPerWord;
        size_t last_word = (get_start_index(seq_id) + end) / kCharsPerWord;
        int first_shift = (get_start_index(seq_id) + begin) % kCharsPerWord * 2;

        s.clear();

        if (end < begin) {
            return;
        }

        if (first_shift == 0) {
            for (size_t i = first_word; i <= last_word; ++i) {
                s.push_back(packed_seq[i]);
            }
        }
        else {
            for (size_t i = first_word; i < last_word; ++i) {
                s.push_back((packed_seq[i] << first_shift) | (packed_seq[i + 1] >> (kBitsPerWord - first_shift)));
            }

            if (kCharsPerWord * s.size() < (unsigned)end - begin + 1) {
                s.push_back(packed_seq[last_word] << first_shift);
            }
        }

        unsigned shift_clean = kBitsPerWord - (end - begin + 1) * 2 % kBitsPerWord;

        if (shift_clean != kBitsPerWord) {
            s.back() >>= shift_clean;
            s.back() <<= shift_clean;
        }
    }
};

#endif