//
// Created by vout on 3/3/2018.
//

#ifndef KMLIB_RNS_H
#define KMLIB_RNS_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
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

}

enum rnsmode {
  kRankOnly,
  kRandAndSelect
};

template<
    unsigned BaseSize,
    unsigned AlphabetSize,
    unsigned Mode = rnsmode::kRandAndSelect,
    unsigned BasePerL2Interval = 65536 / BaseSize,
    unsigned BasePerL1Interval = 1024 / BaseSize,
    unsigned SelectSampleSize = 4096
>
class RankAndSelect {
 public:
  using ull_t = uint64_t;
  using interval_t = uint32_t;
  static const unsigned kBitsPerByte = 8;
  static const unsigned kBitsPerULL = sizeof(ull_t) * kBitsPerByte;
  static const unsigned kBitsPerBase = BaseSize;
  static const unsigned kAlphabetSize = AlphabetSize;
  static const unsigned kBasesPerL1 = BasePerL1Interval;
  static const unsigned kBasesPerL2 = BasePerL2Interval;
  static const unsigned kSelectSampleSize = SelectSampleSize;
  static const unsigned kL1PerL2 = kBasesPerL2 / kBasesPerL1;
  static const unsigned kBasesPerWord = kBitsPerULL / kBitsPerBase;
  // public data, can call directly
  int64_t length;
  int64_t char_frequency[kAlphabetSize];

  RankAndSelect() {
    for (unsigned i = kBitsPerBase == 1 ? 1 : 0; i < kAlphabetSize; ++i) {
      xor_masks_[i] = 0;
      for (unsigned j = 0; j < kBasesPerWord; ++j) {
        xor_masks_[i] |= (ull_t) i << (kBitsPerBase * j);
      }
      xor_masks_[i] = ~xor_masks_[i];
    }
  }

  ~RankAndSelect() = default;

  void Build(ull_t *packed_array, uint64_t length) {
    uint64_t num_l1 = DivCeiling(length, kBasesPerL1) + 1;
    uint64_t num_l2 = DivCeiling(length, kBasesPerL2) + 1;

    for (unsigned c = BaseSize == 1 ? 1 : 0; c < kAlphabetSize; ++c) {
      ull_t *cur_word = packed_array;
      ull_t count = 0;
      l2_occ_[c] = std::vector<uint64_t>(num_l2);
      l1_occ_[c] = std::vector<uint16_t>(num_l1);

      for (uint64_t i = 0; i < length; i += kBasesPerWord, ++cur_word) {
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
      char_frequency[c] = count;

      if (Mode != rnsmode::kRankOnly) {
        rank2itv_[c].reserve(DivCeiling(count, kSelectSampleSize) + 1);
        for (uint64_t i = 0; i < num_l1; ++i) {
          while (rank2itv_[c].size() * kSelectSampleSize < OccValue(c, i)) {
            rank2itv_[c].push_back(i - 1);
          }
        }
        rank2itv_[c].push_back(num_l1 - 1);
        assert(rank2itv_[c].size() <= DivCeiling(count, kSelectSampleSize) + 1);
      }
    }
    packed_array_ = packed_array;
    this->length = length;
  }

  int64_t Rank(int64_t pos) const {
    static_assert(BaseSize == 1, "");
    return InternalRank(1, pos);
  }

  int64_t Rank(uint8_t c, int64_t pos) const {
    static_assert(BaseSize != 1, "");
    return InternalRank(c, pos);
  }

  int64_t Select(int64_t ranking) const {
    static_assert(BaseSize == 1, "");
    return InternalSelect(1, ranking);
  }

  int64_t Select(uint8_t c, int64_t ranking) const {
    static_assert(BaseSize != 1, "");
    return InternalSelect(c, ranking);
  }

  int64_t Pred(uint8_t c, int64_t pos) const {
    // the last c in [0...pos]
    if (GetBaseAt(pos) == c) {
      return pos;
    }
    return InternalSelect(c, InternalRank(c, pos) - 1);
  }

  int64_t Pred(int64_t pos) const {
    return Pred(1, pos);
  }

  int64_t PredLimitedStep(uint8_t c, int64_t pos, int step) const {
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

  int64_t Succ(uint8_t c, int64_t pos) const {
    // the first c in [pos...ReadLength]
    if (GetBaseAt(pos) == c) {
      return pos;
    }
    return InternalSelect(c, InternalRank(c, pos - 1));
  }

  int64_t Succ(int64_t pos) const {
    return Succ(1, pos);
  }

  int64_t SuccLimitedStep(uint8_t c, int64_t pos, int step) const {
    // the first c in [pos, pos+step], return pos+step+1 if not exist
    int64_t end = pos - step;
    if (end >= length) {
      end = length;
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
  unsigned CountCharInWord(uint8_t c, ull_t x, ull_t mask = ull_t(-1)) const {
    if (BaseSize != 1) {
      x ^= xor_masks_[c];
      x = internal::PackToLowestBit<BaseSize>(x);
      x &= kPopcountMask;
    }
    return __builtin_popcountll(x & mask);
  }

  unsigned SelectInWord(uint8_t c, int num_c, ull_t x) const {
    if (BaseSize != 1) {
      x ^= xor_masks_[c];
      x = internal::PackToLowestBit<BaseSize>(x);
      x &= kPopcountMask;
    }
    return __builtin_ctzll(_pdep_u64(1ULL << (num_c - 1), x)) / kBitsPerBase;
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
    if (pos >= length - 1) {
      return char_frequency[c];
    }
    ++pos;
    int64_t itv_idx = (pos + kBasesPerL1 / 2 - 1) / kBasesPerL1;
    int64_t sampled_pos = itv_idx * kBasesPerL1;

    if (sampled_pos >= length) {
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
    if (k >= char_frequency[c]) {
      return length;
    } else if (k < 0) {
      return -1;
    }
    // first locate which interval Select(c, k) falls
    interval_t itv_l = rank2itv_[c][k / kSelectSampleSize];
    interval_t interval_r = rank2itv_[c][DivCeiling(k, kSelectSampleSize)];
    PrefetchOcc(c, itv_l);
    while (interval_r > itv_l) {
      interval_t interval_m = (interval_r + itv_l + 1) / 2;
      if (OccValue(c, interval_m) > uint64_t(k)) {
        interval_r = interval_m - 1;
      } else {
        itv_l = interval_m;
      }
    }
    // refined select
    __builtin_prefetch(packed_array_ + itv_l * kBasesPerL1 / kBasesPerWord);
    unsigned remain = k + 1 - OccValue(c, itv_l);
    unsigned exceed = (itv_l + 1) * kBasesPerL1 >= length ? kBasesPerL1 :
                      (OccValue(c, itv_l + 1) - (k + 1));
    if (remain <= exceed * 2) {
      return SelectFwd(c, itv_l, remain);
    } else {
      return SelectBwd(c, itv_l, exceed);
    }
  }

  uint64_t RankFwd(uint8_t c, interval_t itv, uint64_t sampled_pos, unsigned n_bases) const {
    unsigned n_words = n_bases / kBasesPerWord;
    ull_t *p = packed_array_ + sampled_pos / kBasesPerWord - n_words - 1;
    __builtin_prefetch(p);
    unsigned n_residual = n_bases % kBasesPerWord;
    unsigned count = 0;
    if (n_residual != 0) {
      ull_t mask = 1 + ~(1ULL << kBitsPerBase * (kBasesPerWord - n_residual));
      count += CountCharInWord(c, p[0], mask);
    }
    for (unsigned i = 1; i <= n_words; ++i) {
      count += CountCharInWord(c, p[i]);
    }
    return OccValue(c, itv) - count;
  }

  uint64_t RankBwd(uint8_t c, interval_t itv, uint64_t sampled_pos, unsigned n_bases) const {
    ull_t *p = packed_array_ + sampled_pos / kBasesPerWord;
    __builtin_prefetch(p);
    unsigned n_words = n_bases / kBasesPerWord;
    unsigned n_residual = n_bases % kBasesPerWord;
    unsigned count = 0;
    for (unsigned i = 0; i < n_words; ++i) {
      count += CountCharInWord(c, p[i]);
    }
    if (n_residual != 0) {
      ull_t mask = (1ULL << kBitsPerBase * n_residual) - 1;
      count += CountCharInWord(c, p[n_words], mask);
    }
    return OccValue(c, itv) + count;
  }

  uint64_t SelectFwd(uint8_t c, interval_t itv, unsigned remain) const {
    uint64_t pos = (uint64_t) itv * kBasesPerL1;
    ull_t *begin = packed_array_ + pos / kBasesPerWord;
    ull_t *p = begin;
    unsigned popcnt;
    for (; (popcnt = CountCharInWord(c, *p)) < remain; remain -= popcnt, ++p);
    return pos + (p - begin) * kBasesPerWord + SelectInWord(c, remain, *p);
  }

  uint64_t SelectBwd(uint8_t c, interval_t itv_l, unsigned exceed) const {
    uint64_t pos = (uint64_t) (itv_l + 1) * kBasesPerL1;
    ull_t *end = packed_array_ + pos / kBasesPerWord - 1;
    ull_t *p = end;
    unsigned popcnt;
    for (; (popcnt = CountCharInWord(c, *p)) <= exceed; exceed -= popcnt, --p);
    return pos - kBasesPerWord * (end - p)
        - (kBasesPerWord - SelectInWord(c, popcnt - exceed, *p));
  }

  uint8_t GetBaseAt(uint64_t i) const {
    return (*(packed_array_ + i / kBasesPerWord)
            >> (i % kBasesPerWord * kBitsPerBase)) & ((1 << kBitsPerBase) - 1);
  }

  template<typename T1, typename T2>
  T1 DivCeiling(T1 x, T2 y) const {
    return (x + y - 1) / y;
  };

  static const ull_t kPopcountMask = internal::PopcountMask<BaseSize>::value;
  // main memory for the structure
  ull_t *packed_array_;
  // sampled structure for rank and select
  // two level sampling for rank (occ value)
  // call the function OccValue(c, i) to get the number of c's
  // in packed_array_[0...i*kBasesPerL1-1]
  // sampling for select
  // rank_interval_lookup_[c][i]=j: the jth interval (0 based)
  // contains the (i*kSelectSampleSize)th (0 based) c
  // i.e. OccValue(c, j)<=i*kSelectSampleSize and OccValue(c, j+1)>i*kSelectSampleSize
  std::vector<interval_t> rank2itv_[kAlphabetSize];
  std::vector<uint16_t> l1_occ_[kAlphabetSize]; // level 1 OCC
  std::vector<uint64_t> l2_occ_[kAlphabetSize]; // level 2 OCC
  ull_t xor_masks_[kAlphabetSize];
  // e.g. if c = 0110(2), popcount_xorers_[c] = 1001 1001 1001 1001...(2),
  // to make all c's in a word 1111
  static_assert((1 << kBitsPerBase) >= kAlphabetSize, "");
  static_assert(kBitsPerULL % kBitsPerBase == 0, "");
  static_assert(kBitsPerBase <= 8, "");
  static_assert(kBasesPerL2 <= 65536, "");
};

} // namespace kmlib

using RankAndSelect4Bits = kmlib::RankAndSelect<4, 9>;
using RankAndSelect1Bit = kmlib::RankAndSelect<1, 2>;
using Rank1Bit = kmlib::RankAndSelect<1, 2, kmlib::rnsmode::kRankOnly>;

#endif // KMLIB_RNS_H
