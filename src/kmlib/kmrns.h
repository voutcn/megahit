//
// Created by vout on 3/3/2018.
//

#ifndef KMLIB_RNS_H
#define KMLIB_RNS_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <x86intrin.h>
#include <vector>
#include <new>
#include <limits>

namespace kmlib {

namespace internal {

template<unsigned BaseSize, unsigned Index>
struct SubPopcountMask {
  static_assert(Index % BaseSize == 0, "");
  static const unsigned long long value =
      (SubPopcountMask<BaseSize, Index - BaseSize>::value << BaseSize) | 1ULL;
};

template<unsigned BaseSize>
struct SubPopcountMask<BaseSize, 0> {
  static const unsigned long long value = 0;
};

template<unsigned BaseSize>
struct PopcountMask {
  static const unsigned long long value =
      SubPopcountMask<BaseSize, sizeof(unsigned long long) * 8>::value;
};

template<unsigned BaseSize, typename T, unsigned Index = BaseSize / 2>
inline T PackToLowestBit(T value) {
  if (Index == 0) {
    return value;
  }
  value &= value >> Index;
  return PackToLowestBit<BaseSize, T, Index / 2>(value);
};

using U = unsigned int;
using UL = unsigned long int;
using ULL = unsigned long long int;

#define CHECK_TYPE(T) \
static_assert(std::is_integral<T>::value, "only integral types are supported by popcount"); \
static_assert(sizeof(T) <= sizeof(ULL), "size bigger than unsigned long long not supported"); \
static_assert(!std::is_same<T, bool>::value, "bool type not supported")

template<typename T>
inline unsigned popcount(T val) {
  CHECK_TYPE(T);
  return sizeof(T) <= sizeof(U) ? __builtin_popcount(val) :
         sizeof(T) <= sizeof(UL) ? __builtin_popcountl(val) :
         sizeof(T) <= sizeof(ULL) ? __builtin_popcountll(val) :
         0;
}

template<typename T>
inline unsigned ctz(T val) {
  CHECK_TYPE(T);
  return sizeof(T) <= sizeof(U) ? __builtin_ctz(val) :
         sizeof(T) <= sizeof(UL) ? __builtin_ctzl(val) :
         sizeof(T) <= sizeof(ULL)? __builtin_ctzll(val) :
         0;
}

template<typename Tx, typename Ty>
inline Tx pdep(Tx x, Ty y) {
  CHECK_TYPE(Tx);
  return sizeof(Tx) <= sizeof(U) ? _pdep_u32(x, y) :
         sizeof(Tx) <= sizeof(ULL) ? _pdep_u64(x, y) :
         0;
}

#undef CHECK_TYPE
}

enum rnsmode {
  kRankOnly,
  kRandAndSelect
};

template<
    unsigned BaseSize,
    unsigned AlphabetSize,
    unsigned Mode = rnsmode::kRandAndSelect,
    typename WordType = uint64_t,
    typename IntervalType = uint32_t,
    unsigned BasePerL2Interval = 65536 / BaseSize,
    unsigned BasePerL1Interval = 1024 / BaseSize,
    unsigned SelectSampleSize = 4096
>
class RankAndSelect {
 public:
  using word_type = WordType;
  using interval_type = IntervalType;
  static const unsigned kBitsPerByte = 8;
  static const unsigned kBitsPWord = sizeof(word_type) * kBitsPerByte;
  static const unsigned kBitsPerBase = BaseSize;
  static const unsigned kAlphabetSize = AlphabetSize;
  static const unsigned kBasesPerL1 = BasePerL1Interval;
  static const unsigned kBasesPerL2 = BasePerL2Interval;
  static const unsigned kSelectSampleSize = SelectSampleSize;
  static const unsigned kL1PerL2 = kBasesPerL2 / kBasesPerL1;
  static const unsigned kBasesPerWord = kBitsPWord / kBitsPerBase;

  RankAndSelect() {
    for (unsigned i = kBitsPerBase == 1 ? 1 : 0; i < kAlphabetSize; ++i) {
      xor_masks_[i] = 0;
      for (unsigned j = 0; j < kBasesPerWord; ++j) {
        xor_masks_[i] |= (word_type) i << (kBitsPerBase * j);
      }
      xor_masks_[i] = ~xor_masks_[i];
    }
  }

  ~RankAndSelect() = default;

  void build(word_type *packed_array, uint64_t size) {
    uint64_t num_l1 = DivCeiling(size, kBasesPerL1) + 1;
    uint64_t num_l2 = DivCeiling(size, kBasesPerL2) + 1;

    for (unsigned c = BaseSize == 1 ? 1 : 0; c < kAlphabetSize; ++c) {
      word_type *cur_word = packed_array;
      word_type count = 0;
      l2_occ_[c] = std::vector<uint64_t>(num_l2);
      l1_occ_[c] = std::vector<uint16_t>(num_l1);

      for (uint64_t i = 0; i < size; i += kBasesPerWord, ++cur_word) {
        if (i % kBasesPerL1 == 0) {
          if (i % kBasesPerL2 == 0) {
            l2_occ_[c][i / kBasesPerL2] = count;
          }
          l1_occ_[c][i / kBasesPerL1] = count - l2_occ_[c][i / kBasesPerL2];
        }
        count += CountCharInWord(c, *cur_word);
      }
      l2_occ_[c][num_l2 - 1] = count;
      l1_occ_[c][num_l1 - 1] = count - l2_occ_[c][(num_l1 - 1) / kL1PerL2];
      char_count_[c] = count;

      if (Mode != rnsmode::kRankOnly) {
        rank2itv_[c].reserve(DivCeiling(count, kSelectSampleSize) + 1);
        for (uint64_t i = 0; i < num_l1; ++i) {
          while (rank2itv_[c].size() * kSelectSampleSize < OccValue(c, i)) {
            rank2itv_[c].push_back(i - 1);
          }
        }
        rank2itv_[c].push_back(num_l1 - 1);
      }
    }
    packed_array_ = packed_array;
    this->size_ = size;
  }

  int64_t rank(int64_t pos) const {
    static_assert(BaseSize == 1, "");
    return InternalRank(1, pos);
  }

  int64_t rank(uint8_t c, int64_t pos) const {
    static_assert(BaseSize != 1, "");
    return InternalRank(c, pos);
  }

  int64_t select(int64_t ranking) const {
    static_assert(BaseSize == 1, "");
    return InternalSelect(1, ranking);
  }

  int64_t select(uint8_t c, int64_t ranking) const {
    static_assert(BaseSize != 1, "");
    return InternalSelect(c, ranking);
  }

  int64_t pred(uint8_t c, int64_t pos) const {
    // the last c in [0...pos]
    if (GetBaseAt(pos) == c) {
      return pos;
    }
    return InternalSelect(c, InternalRank(c, pos) - 1);
  }

  int64_t pred(int64_t pos) const {
    static_assert(BaseSize == 1, "");
    return pred(1, pos);
  }

  int64_t pred(uint8_t c, int64_t pos, int step) const {
    // the last c in [pos-step, pos], return pos-step-1 if not exist
    int64_t end = pos - step;
    if (end < 0) {
      end = 0;
    }
    while (pos >= end) {
      if (GetBaseAt(pos) == c) {
        return pos;
      }
      --pos;
    }
    return pos;
  }

  int64_t succ(uint8_t c, int64_t pos) const {
    // the first c in [pos...ReadLength]
    if (GetBaseAt(pos) == c) {
      return pos;
    }
    return InternalSelect(c, InternalRank(c, pos - 1));
  }

  int64_t succ(int64_t pos) const {
    return succ(1, pos);
  }

  int64_t succ(uint8_t c, int64_t pos, int step) const {
    // the first c in [pos, pos+step], return pos+step+1 if not exist
    int64_t end = pos - step;
    if (end >= size_) {
      end = size_;
    }
    while (pos <= end) {
      if (GetBaseAt(pos) == c) {
        return pos;
      }
      ++pos;
    }
    return pos;
  }

 private:
  unsigned CountCharInWord(uint8_t c, word_type x, word_type mask = word_type(-1)) const {
    if (BaseSize != 1) {
      x ^= xor_masks_[c];
      x = internal::PackToLowestBit<BaseSize>(x);
      x &= kPopcountMask;
    }
    return internal::popcount(x & mask);
  }

  unsigned SelectInWord(uint8_t c, int num_c, word_type x) const {
    if (BaseSize != 1) {
      x ^= xor_masks_[c];
      x = internal::PackToLowestBit<BaseSize>(x);
      x &= kPopcountMask;
    }
#ifdef __BMI2__
    return internal::ctz(internal::pdep(word_type(1) << (num_c - 1), x)) / kBitsPerBase;
#else
    unsigned trailing_zeros = 0;
    while (num_c > 0) {
      trailing_zeros = internal::ctz(x);
      x ^= word_type(1) << trailing_zeros;
      --num_c;
    }
    return trailing_zeros / kBitsPerBase;
#endif
  }

  void PrefetchOcc(uint8_t c, int64_t i) const {
    __builtin_prefetch(&l2_occ_[c][i / kL1PerL2], 0);
    __builtin_prefetch(&l1_occ_[c][i], 0);
  }

  uint64_t OccValue(uint8_t c, int64_t i) const {
    return l2_occ_[c][i / kL1PerL2] + l1_occ_[c][i];
  }

 private:
  int64_t InternalRank(uint8_t c, int64_t pos) const {
    // the number of c's in [0...pos]
    if (pos >= size_ - 1) {
      return char_count_[c];
    }
    ++pos;
    int64_t itv_idx = (pos + kBasesPerL1 / 2 - 1) / kBasesPerL1;
    int64_t sampled_pos = itv_idx * kBasesPerL1;

    if (sampled_pos >= size_) {
      sampled_pos -= kBasesPerL1;
      itv_idx--;
    }
    PrefetchOcc(c, itv_idx);

    if (sampled_pos > pos) {
      return RankFwd(c, itv_idx, sampled_pos, sampled_pos - pos);
    } else if (sampled_pos < pos) {
      return RankBwd(c, itv_idx, sampled_pos, pos - sampled_pos);
    } else {
      return OccValue(c, itv_idx);
    }
  }

  int64_t InternalSelect(uint8_t c, int64_t k) const {
    static_assert(Mode != rnsmode::kRankOnly,
                  "cannot select on rank only struct");
    // return the pos (0-based) of the kth (0-based) c
    if (k >= char_count_[c]) {
      return size_;
    } else if (k < 0) {
      return -1;
    }
    // first locate which interval Select(c, k) falls
    interval_type interval_l = rank2itv_[c][k / kSelectSampleSize];
    interval_type interval_r = rank2itv_[c][DivCeiling(k, kSelectSampleSize)];
    PrefetchOcc(c, interval_l);
    while (interval_r > interval_l) {
      interval_type interval_m = (interval_r + interval_l + 1) / 2;
      if (OccValue(c, interval_m) > uint64_t(k)) {
        interval_r = interval_m - 1;
      } else {
        interval_l = interval_m;
      }
    }
    // refined select
    __builtin_prefetch(packed_array_ + interval_l * kBasesPerL1 / kBasesPerWord);
    unsigned remain = k + 1 - OccValue(c, interval_l);
    unsigned exceed = (interval_l + 1) * kBasesPerL1 >= size_ ? kBasesPerL1 :
                      (OccValue(c, interval_l + 1) - (k + 1));
    if (remain <= exceed * 2) {
      return SelectFwd(c, interval_l, remain);
    } else {
      return SelectBwd(c, interval_l, exceed);
    }
  }

  uint64_t RankFwd(uint8_t c, interval_type itv, uint64_t sampled_pos, unsigned n_bases) const {
    unsigned n_words = n_bases / kBasesPerWord;
    word_type *p = packed_array_ + sampled_pos / kBasesPerWord - n_words - 1;
    __builtin_prefetch(p);
    unsigned n_residual = n_bases % kBasesPerWord;
    unsigned count = 0;
    if (n_residual != 0) {
      word_type mask = 1 + ~(1ULL << kBitsPerBase * (kBasesPerWord - n_residual));
      count += CountCharInWord(c, p[0], mask);
    }
    for (unsigned i = 1; i <= n_words; ++i) {
      count += CountCharInWord(c, p[i]);
    }
    return OccValue(c, itv) - count;
  }

  uint64_t RankBwd(uint8_t c, interval_type itv, uint64_t sampled_pos, unsigned n_bases) const {
    word_type *p = packed_array_ + sampled_pos / kBasesPerWord;
    __builtin_prefetch(p);
    unsigned n_words = n_bases / kBasesPerWord;
    unsigned n_residual = n_bases % kBasesPerWord;
    unsigned count = 0;
    for (unsigned i = 0; i < n_words; ++i) {
      count += CountCharInWord(c, p[i]);
    }
    if (n_residual != 0) {
      word_type mask = (1ULL << kBitsPerBase * n_residual) - 1;
      count += CountCharInWord(c, p[n_words], mask);
    }
    return OccValue(c, itv) + count;
  }

  uint64_t SelectFwd(uint8_t c, interval_type itv, unsigned remain) const {
    uint64_t pos = (uint64_t) itv * kBasesPerL1;
    word_type *begin = packed_array_ + pos / kBasesPerWord;
    word_type *p = begin;
    unsigned popcnt;
    for (; (popcnt = CountCharInWord(c, *p)) < remain; remain -= popcnt, ++p);
    return pos + (p - begin) * kBasesPerWord + SelectInWord(c, remain, *p);
  }

  uint64_t SelectBwd(uint8_t c, interval_type itv_l, unsigned exceed) const {
    uint64_t pos = (uint64_t) (itv_l + 1) * kBasesPerL1;
    word_type *end = packed_array_ + pos / kBasesPerWord - 1;
    word_type *p = end;
    unsigned popcnt;
    for (; (popcnt = CountCharInWord(c, *p)) <= exceed; exceed -= popcnt, --p);
    return pos - kBasesPerWord * (end - p)
        - (kBasesPerWord - SelectInWord(c, popcnt - exceed, *p));
  }

  uint8_t GetBaseAt(uint64_t i) const {
    return (packed_array_[i / kBasesPerWord]
        >> (i % kBasesPerWord * kBitsPerBase)) & ((1 << kBitsPerBase) - 1);
  }

  template<typename T1, typename T2>
  T1 DivCeiling(T1 x, T2 y) const {
    return (x + y - 1) / y;
  };

  static const word_type kPopcountMask = internal::PopcountMask<BaseSize>::value;
  int64_t size_;
  int64_t char_count_[kAlphabetSize];
  // main memory for the structure
  word_type *packed_array_;
  // sampled structure for rank and select
  // two level sampling for rank (occ value)
  // call the function OccValue(c, i) to get the number of c's
  // in packed_array_[0...i*kBasesPerL1-1]
  // sampling for select
  // rank_interval_lookup_[c][i]=j: the jth interval (0 based)
  // contains the (i*kSelectSampleSize)th (0 based) c
  // i.e. OccValue(c, j)<=i*kSelectSampleSize and OccValue(c, j+1)>i*kSelectSampleSize
  std::vector<interval_type> rank2itv_[kAlphabetSize];
  std::vector<uint16_t> l1_occ_[kAlphabetSize]; // level 1 OCC
  std::vector<uint64_t> l2_occ_[kAlphabetSize]; // level 2 OCC
  word_type xor_masks_[kAlphabetSize];
  // e.g. if c = 0110(2), popcount_xorers_[c] = 1001 1001 1001 1001...(2),
  // to make all c's in a word 1111
  static_assert((1 << kBitsPerBase) >= kAlphabetSize, "");
  static_assert(kBitsPWord % kBitsPerBase == 0, "");
  static_assert(kBitsPerBase <= 8, "");
  static_assert(kBasesPerL2 <= 65536, "");
};

} // namespace kmlib

using RankAndSelect4Bits = kmlib::RankAndSelect<4, 9>;
using RankAndSelect1Bit = kmlib::RankAndSelect<1, 2>;
using Rank1Bit = kmlib::RankAndSelect<1, 2, kmlib::rnsmode::kRankOnly>;

#endif // KMLIB_RNS_H
