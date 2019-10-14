//
// Created by dinghua.li on 1/27/18.
//

#ifndef MEGAHIT_KMER_H
#define MEGAHIT_KMER_H

#include <algorithm>
#include <cstdint>
#include <cstring>
#include "kmlib/kmbit.h"

/**
 * @brief It represents a k-mer. The value of k is limited by the number of
 * words used. The maximum value can be calculated by max_size().
 * The value of k is not stored, many functions require a explict parameter k to
 * work
 */
template <unsigned NWords = 4, class TWord = uint64_t>
class Kmer {
 public:
  using word_type = TWord;
  static const unsigned kNumWords = NWords;

  Kmer() {
    static_assert(sizeof(*this) == sizeof(TWord) * NWords, "");
    std::memset(data_, 0, sizeof(data_));
  }

  Kmer(const Kmer &kmer) { std::memcpy(data_, kmer.data_, sizeof(data_)); }

  Kmer(const word_type *seq, unsigned offset, unsigned k) {
    InitFromPtr(seq, offset, k);
  }

  void InitFromPtr(const word_type *seq, unsigned offset, unsigned k) {
    seq += offset / kCharsPerWord;
    offset %= kCharsPerWord;
    offset <<= 1;

    unsigned used_words = (k + kCharsPerWord - 1) / kCharsPerWord;

    if (offset == 0) {
      std::memcpy(data_, seq, sizeof(word_type) * used_words);
    } else {
      for (unsigned i = 0; i < used_words - 1; ++i) {
        data_[i] = (seq[i] << offset) | (seq[i + 1] >> (kBitsPerWord - offset));
      }

      data_[used_words - 1] = seq[used_words - 1] << offset;

      if (offset + k * 2 > kBitsPerWord * used_words) {
        data_[used_words - 1] |=
            seq[used_words] >> (kBitsPerWord - offset) % kBitsPerWord;
      }
    }

    if (k % kCharsPerWord != 0) {
      uint32_t clean_shift = (kCharsPerWord - k % kCharsPerWord) << 1;
      data_[used_words - 1] = data_[used_words - 1] >> clean_shift
                                                           << clean_shift;
    }

    memset(data_ + used_words, 0, sizeof(word_type) * (kNumWords - used_words));
  }

  ~Kmer() = default;

  const word_type *data() const { return data_; }

  Kmer &operator=(const Kmer &kmer) {
    std::memcpy(data_, kmer.data_, sizeof(data_));
    return *this;
  }

  bool operator<(const Kmer &kmer) const {
    for (unsigned i = 0; i < kNumWords; ++i) {
      if (data_[i] != kmer.data_[i]) return data_[i] < kmer.data_[i];
    }

    return false;
  }

  bool operator>(const Kmer &kmer) const {
    for (unsigned i = 0; i < kNumWords; ++i) {
      if (data_[i] != kmer.data_[i]) return data_[i] > kmer.data_[i];
    }

    return false;
  }

  bool operator==(const Kmer &kmer) const {
    return memcmp(data_, kmer.data_, sizeof(data_)) == 0;
  }

  bool operator!=(const Kmer &kmer) const {
    return memcmp(data_, kmer.data_, sizeof(data_)) != 0;
  }

  int cmp(const Kmer &kmer, unsigned k) const {
    unsigned used_words = (k + kCharsPerWord - 1) / kCharsPerWord;
    for (unsigned i = 0; i < used_words; ++i) {
      if (data_[i] < kmer.data_[i]) {
        return -1;
      } else if (data_[i] > kmer.data_[i]) {
        return 1;
      }
    }
    return 0;
  }

  const Kmer &ReverseComplement(unsigned k) {
    uint32_t used_words = (k + kCharsPerWord - 1) / kCharsPerWord;
    for (unsigned i = 0; i < used_words; ++i) {
      data_[i] = kmlib::bit::ReverseComplement<2>(data_[i]);
    }

    for (unsigned i = 0; i + i < used_words; ++i) {
      std::swap(data_[i], data_[used_words - 1 - i]);
    }

    if ((k % kCharsPerWord) != 0) {
      unsigned offset = (kCharsPerWord - k % kCharsPerWord) << 1;

      for (unsigned i = 0; i + 1 < used_words; ++i) {
        data_[i] =
            (data_[i] << offset) | (data_[i + 1] >> (kBitsPerWord - offset));
      }

      data_[used_words - 1] <<= offset;
    }

    return *this;
  }

  void ShiftAppend(uint8_t ch, unsigned k) {
    ch &= 3;
    uint32_t used_words =
        std::min(NWords, (k + kCharsPerWord - 1) / kCharsPerWord);

    for (unsigned i = 0; i + 1 < used_words; ++i) {
      data_[i] = (data_[i] << 2) | (data_[i + 1] >> (kBitsPerWord - 2));
    }

    data_[used_words - 1] =
        (data_[used_words - 1] << 2) |
        (word_type(ch) << ((kCharsPerWord - 1 - (k - 1) % kCharsPerWord) << 1));
  }

  void ShiftPreappend(uint8_t ch, unsigned k) {
    ch &= 3;
    uint32_t used_words = (k + kCharsPerWord - 1) / kCharsPerWord;
    used_words = std::min(NWords, used_words);

    for (int i = used_words - 1; i > 0; --i) {
      data_[i] = (data_[i] >> 2) | (data_[i - 1] << (kBitsPerWord - 2));
    }

    data_[0] = (data_[0] >> 2) | (word_type(ch) << (kBitsPerWord - 2));

    if (k % kCharsPerWord != 0) {
      uint32_t clean_shift = (kCharsPerWord - k % kCharsPerWord) << 1;
      data_[used_words - 1] = data_[used_words - 1] >> clean_shift
                                                           << clean_shift;
    }
  }

  bool IsPalindrome(unsigned k) const {
    if (k % 2 == 1) return false;
    Kmer kmer(*this);
    return kmer.ReverseComplement(k) == *this;
  }

  Kmer unique_format(unsigned k) const {
    Kmer rev_comp = *this;
    rev_comp.ReverseComplement(k);
    return (this->cmp(rev_comp, k) <= 0 ? *this : rev_comp);
  }

  uint8_t operator[](uint32_t index) const { return GetBase(index); }

  uint8_t GetBase(uint32_t index) const {
    return (data_[index / kCharsPerWord] >>
            ((kCharsPerWord - 1 - index % kCharsPerWord) << 1)) &
           3;
  }

  void SetBase(uint32_t index, uint8_t ch) {
    ch &= 3;
    unsigned offset = (kCharsPerWord - 1 - index % kCharsPerWord) << 1;
    data_[index / kCharsPerWord] =
        (data_[index / kCharsPerWord] & ~(word_type(3) << offset)) |
        (word_type(ch) << offset);
  }

  void swap(Kmer &kmer) {
    if (this != &kmer) {
      for (unsigned i = 0; i < kNumWords; ++i) {
        std::swap(data_[i], kmer.data_[i]);
      }
    }
  }

  void clear() { memset(data_, 0, sizeof(word_type) * kNumWords); }

  static uint32_t max_size() { return kMaxSize; }

  unsigned kth_byte(unsigned k) const {
    return (data_[k / sizeof(word_type)] >> k % sizeof(word_type) * 8) & 0xFF;
  }

  static const uint32_t kBitsPerWord = sizeof(word_type) * 8;
  static const uint32_t kCharsPerWord = kBitsPerWord / 2;
  static const uint32_t kBitsForKmer = kNumWords * kBitsPerWord;
  static const uint32_t kMaxSize = kBitsForKmer / 2;
  static const unsigned n_bytes = sizeof(word_type) * kNumWords;

 private:
  word_type data_[kNumWords];
};

namespace std {
template <const unsigned NumWords, typename T>
inline void swap(Kmer<NumWords, T> &kmer1, Kmer<NumWords, T> &kmer2) {
  kmer1.swap(kmer2);
}
}  // namespace std

#endif  // MEGAHIT_KMER_H
