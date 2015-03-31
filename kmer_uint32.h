/*
 *  MEGAHIT
 *  Copyright (C) 2014 The University of Hong Kong
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

#ifndef KMER_UINT32_H_
#define KMER_UINT32_H_

#include "definitions.h"
#include "utils.h"

/**
 * @brief a kmer struct for sdbg_builder (for efficiency). Not well OO.
 */

struct KmerUint32 {
    static const int kNumWords = 6;
    edge_word_t data_[kNumWords];
    char kmer_k;

    KmerUint32() {
        clean();
    }

    KmerUint32(edge_word_t *read_ptr, int k) {
        init(read_ptr, k);
    }

    void init(edge_word_t *read_ptr, int k) {
        this->kmer_k = k;
        int num_words = DivCeiling(k, kCharsPerEdgeWord);
        int bits_to_clean = (kCharsPerEdgeWord - k % kCharsPerEdgeWord) * kBitsPerEdgeChar;
        for (int i = 0; i < num_words; ++i) {
            data_[i] = read_ptr[i];
        }
        if (bits_to_clean > 0) {
            data_[num_words - 1] >>= bits_to_clean;
            data_[num_words - 1] <<= bits_to_clean;
        }
        for (int i = num_words; i < kNumWords; ++i) {
            data_[i] = 0;
        }
    }

    void init(edge_word_t *read_ptr, int k, int spacing) {
        this->kmer_k = k;
        int num_words = DivCeiling(k, kCharsPerEdgeWord);
        int bits_to_clean = (kCharsPerEdgeWord - k % kCharsPerEdgeWord) * kBitsPerEdgeChar;
        for (int i = 0; i < num_words; ++i) {
            data_[i] = read_ptr[i * spacing];
        }
        if (bits_to_clean != kBitsPerEdgeWord) {
            data_[num_words - 1] >>= bits_to_clean;
            data_[num_words - 1] <<= bits_to_clean;
        }
        for (int i = num_words; i < kNumWords; ++i) {
            data_[i] = 0;
        }
    }

    void clean() {
        kmer_k = 0;
        memset(data_, 0, sizeof(data_));
    }

    int cmp(const KmerUint32 &rhs) const {
        if (kmer_k != rhs.kmer_k) return kmer_k - rhs.kmer_k;
        for (int i = 0; i < kNumWords; ++i) {
            if (data_[i] != rhs.data_[i]) {
                return data_[i] < rhs.data_[i] ? -1 : 1;
            }
        }
        return 0;
    }

    bool operator< (const KmerUint32 &rhs) const {
        return cmp(rhs) < 0;
    }

    bool operator== (const KmerUint32 &rhs) const {
        return cmp(rhs) == 0;
    }

    void ShiftLeft() {
        for (int i = 0; i < kNumWords - 1; ++i) {
            data_[i] = (data_[i] << kBitsPerEdgeChar) | (data_[i + 1] >> (kCharsPerEdgeWord - 1) * kBitsPerEdgeChar);
        }
        data_[kNumWords - 1] <<= kBitsPerEdgeChar;
        kmer_k--;
    }

    void ShiftLeftAppend(int x) {
        int num_words = DivCeiling(kmer_k, kCharsPerEdgeWord);
        for (int i = 0; i < num_words - 1; ++i) {
            data_[i] = (data_[i] << kBitsPerEdgeChar) | (data_[i + 1] >> (kCharsPerEdgeWord - 1) * kBitsPerEdgeChar);
        }
        data_[num_words - 1] <<= kBitsPerEdgeChar;
        data_[num_words - 1] |= x << (kCharsPerEdgeWord - 1 - (kmer_k - 1) % kCharsPerEdgeWord) * kBitsPerEdgeChar;
    }

    void ShiftRightAppend(int x) {
        for (int i = kNumWords - 1; i > 0; --i) {
            data_[i] = (data_[i] >> kBitsPerEdgeChar) | ((data_[i - 1] & 3) << (kCharsPerEdgeWord - 1) * kBitsPerEdgeChar);
        }
        data_[0] = (data_[0] >> kBitsPerEdgeChar) | (x << (kCharsPerEdgeWord - 1) * kBitsPerEdgeChar);
        // clean the last char
        data_[kmer_k / kCharsPerEdgeWord] &= ~(3 << (kCharsPerEdgeWord - 1 - kmer_k % kCharsPerEdgeWord) * kBitsPerEdgeChar);
    }

    const KmerUint32 &ReverseComplement() {
        int which_word = (kmer_k - 1) / kCharsPerEdgeWord;
        int word_offset = (kmer_k - 1) % kCharsPerEdgeWord;
        edge_word_t new_data[kNumWords];
        if (word_offset == kCharsPerEdgeWord - 1) {
            // aligned
            for (int i = 0; i <= which_word; ++i) {
                new_data[i] = ~ mirror(data_[which_word - i]);
            }
        } else {
            // not aligned
            int bit_offset = (kCharsPerEdgeWord - 1 - word_offset) * kBitsPerEdgeChar;
            for (int i = 0; i < which_word; ++i) {
                new_data[i] = (data_[which_word - i] >> bit_offset) | (data_[which_word - i - 1] << (kBitsPerEdgeWord - bit_offset));
                new_data[i] = ~ mirror(new_data[i]);
            }
            new_data[which_word] = ~ mirror(data_[0] >> bit_offset);
        }
        init(new_data, kmer_k);
        return *this;
    }

    void Append(int x) {
        data_[kmer_k / kCharsPerEdgeWord] |= x << (kCharsPerEdgeWord - 1 - kmer_k % kCharsPerEdgeWord) * kBitsPerEdgeChar;
        kmer_k++;
    }

    int operator[] (int i) {
        return (data_[i / kCharsPerEdgeWord] >> (kCharsPerEdgeWord - 1 - i % kCharsPerEdgeWord) * kBitsPerEdgeChar) & 3;
    }
};

#endif // KMER_UINT32_H_