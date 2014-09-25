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

/**
 * @file kmer.h
 * @brief Kmer Class.  
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-02
 */

#ifndef __BASIC_KMER_H_

#define __BASIC_KMER_H_

#include <stdint.h>

#include <algorithm>
#include <cstring>

#include "bit_operation.h"


/**
 * @brief It represents a k-mer. The value of k is limited by the number of 
 * uint64 words used. The maximum value can be calculated by max_size().
 */
template <uint32_t kNumUint64 = 4>
class Kmer
{
public:
    Kmer() 
    { std::memset(data_, 0, sizeof(uint64_t) * kNumUint64); }

    Kmer(const Kmer &kmer)
    { std::memcpy(data_, kmer.data_, sizeof(uint64_t) * kNumUint64); }

    explicit Kmer(uint32_t size) 
    { std::memset(data_, 0, sizeof(uint64_t) * kNumUint64); resize(size); }

    ~Kmer() {}

    const Kmer &operator = (const Kmer &kmer)
    { std::memcpy(data_, kmer.data_, sizeof(uint64_t) * kNumUint64); return *this; }

    bool operator <(const Kmer &kmer) const
    {
        for (int i = kNumUint64-1; i >= 0; --i)
        {
            if (data_[i] != kmer.data_[i])
                return data_[i] < kmer.data_[i];
        }
        return false;
    }

    bool operator >(const Kmer &kmer) const
    {
        for (int i = kNumUint64-1; i >= 0; --i)
        {
            if (data_[i] != kmer.data_[i])
                return data_[i] > kmer.data_[i];
        }
        return false;
    }

    bool operator ==(const Kmer &kmer) const
    {
        for (unsigned i = 0; i < kNumUint64; ++i)
        {
            if (data_[i] != kmer.data_[i])
                return false;
        }
        return true;
    }

    bool operator !=(const Kmer &kmer) const
    {
        for (unsigned i = 0; i < kNumUint64; ++i)
        {
            if (data_[i] != kmer.data_[i])
                return true;
        }
        return false;
    }

    const Kmer &ReverseComplement()
    {
        uint32_t kmer_size = size();
        uint32_t used_words = (kmer_size + 31) >> 5;

        resize(0);

        for (unsigned i = 0; i < used_words; ++i)
            bit_operation::ReverseComplement(data_[i]);

        for (unsigned i = 0; i < (used_words >> 1); ++i)
            std::swap(data_[i], data_[used_words-1-i]);

        if ((kmer_size & 31) != 0)
        {
            unsigned offset = (32 - (kmer_size & 31)) << 1;
            for (unsigned i = 0; i+1 < used_words; ++i)
                data_[i] = (data_[i] >> offset) | data_[i+1] << (64 - offset);
            data_[used_words-1] >>= offset;
        }

        resize(kmer_size);

        return *this;
    }

    void ShiftAppend(uint8_t ch)
    {
        ch &= 3;
        uint32_t kmer_size = size();
        uint32_t used_words = (kmer_size + 31) >> 5;

        resize(0);

        for (unsigned i = 0; i+1 < used_words ; ++i)
            data_[i] = (data_[i] >> 2) | (data_[i+1] << 62);
        data_[used_words-1] = (data_[used_words-1] >> 2) | (uint64_t(ch) << (((kmer_size - 1) & 31) << 1));

        resize(kmer_size);
    }

    void ShiftPreappend(uint8_t ch)
    {
        ch &= 3;
        uint32_t kmer_size = size();
        uint32_t used_words = (kmer_size + 31) >> 5;

        resize(0);

        for (int i = used_words-1; i > 0; --i)
            data_[i] = (data_[i] << 2) | (data_[i-1] >> 62);
        data_[0] = (data_[0] << 2) | ch;

        if ((kmer_size & 31) != 0)
            data_[used_words-1] &= (1ULL << ((kmer_size & 31) << 1)) - 1;

        resize(kmer_size);
    }

    bool IsPalindrome() const
    {
        Kmer kmer(*this);
        return kmer.ReverseComplement() == *this;
    }

    uint64_t hash() const
    {
        uint64_t key = 0;
        for (unsigned i = 0; i < kNumUint64; ++i)
            key ^= data_[i];
        return (key * 1299709 + 104729) % 323780508946331ULL;
    }

    Kmer unique_format() const
    {
        Kmer rev_comp = *this;
        rev_comp.ReverseComplement();
        return (*this < rev_comp ? *this : rev_comp);
    }

    uint8_t operator [] (uint32_t index) const 
    { return (data_[index>>5] >> ((index & 31) << 1)) & 3; }

    uint8_t get_base(uint32_t index) const
    { return (data_[index>>5] >> ((index & 31) << 1)) & 3; }

    void set_base(uint32_t index, uint8_t ch)
    {
        ch &= 3;
        unsigned offset = (index & 31) << 1;
        data_[index>>5] = (data_[index>>5] & ~(3ULL << offset)) | (uint64_t(ch) << offset);
    }

    void swap(Kmer &kmer)
    {
        if (this != &kmer)
        {
            for (unsigned i = 0; i < kNumUint64; ++i)
                std::swap(data_[i], kmer.data_[i]);
        }
    }

    uint32_t size() const
    { return data_[kNumUint64-1] >> (64 - kBitsForSize); }
    void resize(uint32_t new_size)
    { data_[kNumUint64-1] = ((data_[kNumUint64-1] << kBitsForSize) >> kBitsForSize) | (uint64_t(new_size) << (64 - kBitsForSize)); }

    void clear()
    { 
        uint32_t kmer_size = size();
        memset(data_, 0, sizeof(uint64_t) * kNumUint64); 
        resize(kmer_size);
    }

    static uint32_t max_size()
    { return kMaxSize; }

    static const uint32_t kBitsForSize = ((kNumUint64 <= 2) ? 6 : ((kNumUint64 <= 8) ? 8 : 16));
    static const uint32_t kBitsForKmer = (kNumUint64 * 64 - kBitsForSize);
    static const uint32_t kMaxSize = kBitsForKmer / 2;

    uint64_t data_[kNumUint64];
};

namespace std
{
template <const uint32_t kNumUint64> inline void swap(Kmer<kNumUint64> &kmer1, Kmer<kNumUint64> &kmer2) { kmer1.swap(kmer2); }
}

#endif

