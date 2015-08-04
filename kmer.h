/**
 * @file kmer.h
 * @brief Kmer Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-02
 * @modified by Dinghua Li
 */

#ifndef KMER_H__

#define KMER_H__

#include <stdint.h>

#include <algorithm>
#include <cstring>
#include "city.h"
#include "bit_operation.h"

/**
 * @brief It represents a k-mer. The value of k is limited by the number of
 * words used. The maximum value can be calculated by max_size().
 * The value of k is not stored, many functions require a explict parameter k to work
 */
template <unsigned kNumWords = 4, typename T = uint64_t>
struct Kmer {
    typedef T word_t;

    Kmer() {
        std::memset(data_, 0, sizeof(word_t) * kNumWords);
    }

    Kmer(const Kmer &kmer) {
        std::memcpy(data_, kmer.data_, sizeof(word_t) * kNumWords);
    }

    Kmer(word_t *seq, int offset, int k) {
        init(seq, offset, k);
    }

    void init(word_t *seq, int offset, int k) {
        seq += offset / kCharsPerWord;
        offset %= kCharsPerWord;
        offset <<= 1;

        int used_words = (k + kCharsPerWord - 1) / kCharsPerWord;

        if (offset == 0) {
            std::memcpy(data_, seq, sizeof(word_t) * used_words);
        }
        else {
            for (int i = 0; i < used_words - 1; ++i) {
                data_[i] = (seq[i] << offset) | (seq[i + 1] >> (kBitsPerWord - offset));
            }

            data_[used_words - 1] = seq[used_words - 1] << offset;

            if (offset + k * 2 > (int)kBitsPerWord * used_words) {
                data_[used_words - 1] |= seq[used_words] >> (kBitsPerWord - offset) % kBitsPerWord;
            }
        }

        if (k % kCharsPerWord != 0) {
            uint32_t clean_shift = (kCharsPerWord - k % kCharsPerWord) << 1;
            data_[used_words - 1] = data_[used_words - 1] >> clean_shift << clean_shift;
        }

        memset(data_ + used_words, 0, sizeof(word_t) * (kNumWords - used_words));
    }

    ~Kmer() {}

    const Kmer &operator = (const Kmer &kmer) {
        std::memcpy(data_, kmer.data_, sizeof(word_t) * kNumWords);
        return *this;
    }

    bool operator <(const Kmer &kmer) const {
        for (unsigned i = 0; i < kNumWords; --i) {
            if (data_[i] != kmer.data_[i])
                return data_[i] < kmer.data_[i];
        }

        return false;
    }

    bool operator >(const Kmer &kmer) const {
        for (unsigned i = 0; i < kNumWords; --i) {
            if (data_[i] != kmer.data_[i])
                return data_[i] > kmer.data_[i];
        }

        return false;
    }

    bool operator ==(const Kmer &kmer) const {
        for (unsigned i = 0; i < kNumWords; ++i) {
            if (data_[i] != kmer.data_[i])
                return false;
        }

        return true;
    }

    bool operator !=(const Kmer &kmer) const {
        for (unsigned i = 0; i < kNumWords; ++i) {
            if (data_[i] != kmer.data_[i])
                return true;
        }

        return false;
    }

    int cmp(const Kmer &kmer, int k) const {
        int used_words = (k + kCharsPerWord - 1) / kCharsPerWord;

        for (int i = 0; i < used_words; ++i) {
            if (data_[i] < kmer.data_[i]) {
                return -1;
            }
            else if (data_[i] > kmer.data_[i]) {
                return 1;
            }
        }

        return 0;
    }

    const Kmer &ReverseComplement(int k) {
        uint32_t used_words = (k + kCharsPerWord - 1) / kCharsPerWord;

        for (unsigned i = 0; i < used_words; ++i)
            bit_operation::ReverseComplement(data_[i]);

        for (unsigned i = 0; i < (used_words >> 1); ++i)
            std::swap(data_[i], data_[used_words - 1 - i]);

        if ((k % kCharsPerWord) != 0) {
            unsigned offset = (kCharsPerWord - k % kCharsPerWord) << 1;

            for (unsigned i = 0; i + 1 < used_words; ++i)
                data_[i] = (data_[i] << offset) | (data_[i + 1] >> (kBitsPerWord - offset));

            data_[used_words - 1] <<= offset;
        }

        return *this;
    }

    void ShiftAppend(uint8_t ch, int k) {
        ch &= 3;
        uint32_t used_words = (k + kCharsPerWord - 1) / kCharsPerWord;

        for (unsigned i = 0; i + 1 < used_words ; ++i)
            data_[i] = (data_[i] << 2) | (data_[i + 1] >> (kBitsPerWord - 2));

        data_[used_words - 1] = (data_[used_words - 1] << 2) | (word_t(ch) << ((kCharsPerWord - 1 - (k - 1) % kCharsPerWord) << 1));
    }

    void ShiftPreappend(uint8_t ch, int k) {
        ch &= 3;
        uint32_t used_words = (k + kCharsPerWord - 1) / kCharsPerWord;

        for (int i = used_words - 1; i > 0; --i)
            data_[i] = (data_[i] >> 2) | (data_[i - 1] << (kBitsPerWord - 2));

        data_[0] = (data_[0] >> 2) | (word_t(ch) << (kBitsPerWord - 2));

        if (k % kCharsPerWord != 0) {
            uint32_t clean_shift = (kCharsPerWord - k % kCharsPerWord) << 1;
            data_[used_words - 1] = data_[used_words - 1] >> clean_shift << clean_shift;
        }
    }

    bool IsPalindrome(int k) const {
        Kmer kmer(*this);
        return kmer.ReverseComplement(k) == *this;
    }

    uint64_t hash() const {
        return CityHash64((const char *)data_, sizeof(data_[0]) * kNumWords);
    }

    Kmer unique_format(int k) const {
        Kmer rev_comp = *this;
        rev_comp.ReverseComplement(k);
        return (this->cmp(rev_comp, k) <= 0 ? *this : rev_comp);
    }

    uint8_t operator [] (uint32_t index) const {
        return (data_[index / kCharsPerWord] >> ((kCharsPerWord - 1 - index % kCharsPerWord) << 1)) & 3;
    }

    uint8_t get_base(uint32_t index) const {
        return (data_[index / kCharsPerWord] >> ((kCharsPerWord - 1 - index % kCharsPerWord) << 1)) & 3;
    }

    void set_base(uint32_t index, uint8_t ch) {
        ch &= 3;
        unsigned offset = (kCharsPerWord - 1 - index % kCharsPerWord) << 1;
        data_[index / kCharsPerWord] = (data_[index / kCharsPerWord] & ~(word_t(3) << offset)) | (word_t(ch) << offset);
    }

    void swap(Kmer &kmer) {
        if (this != &kmer) {
            for (unsigned i = 0; i < kNumWords; ++i)
                std::swap(data_[i], kmer.data_[i]);
        }
    }

    void clear() {
        memset(data_, 0, sizeof(word_t) * kNumWords);
    }

    static uint32_t max_size() {
        return kMaxSize;
    }

    static const uint32_t kBitsPerWord = sizeof(word_t) * 8;
    static const uint32_t kCharsPerWord = kBitsPerWord / 2;
    static const uint32_t kBitsForKmer = kNumWords * kBitsPerWord;
    static const uint32_t kMaxSize = kBitsForKmer / 2;

    word_t data_[kNumWords];
} __attribute__((packed));


namespace std {
template <const unsigned kNumWords, typename T> inline void swap(Kmer<kNumWords, T> &kmer1, Kmer<kNumWords, T> &kmer2) {
    kmer1.swap(kmer2);
}
}

#endif

