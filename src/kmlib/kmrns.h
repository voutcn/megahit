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

namespace kmlib {

namespace internal {

template <unsigned BaseSize, unsigned Index, typename T>
struct SubPopcountMask {
  static_assert(Index % BaseSize == 0, "");
  static const T value =
      (SubPopcountMask<BaseSize, Index - BaseSize, T>::value << BaseSize) |
      1ULL;
};

template <unsigned BaseSize, typename T>
struct SubPopcountMask<BaseSize, 0, T> {
  static const T value = 0;
};

template <unsigned BaseSize, typename T>
struct PopcountMask {
  static const T value = SubPopcountMask<BaseSize, sizeof(T) * 8, T>::value;
};

template <unsigned BaseSize, typename T, unsigned Index = BaseSize / 2>
inline T PackToLowestBit(T value) {
  if (Index == 0) {
    return value;
  }
  value &= value >> Index;
  return PackToLowestBit<BaseSize, T, Index / 2>(value);
}

using U = unsigned int;
using UL = unsigned long int;
using ULL = unsigned long long int;

#define CHECK_TYPE(T)                                                 \
  static_assert(std::is_integral<T>::value,                           \
                "only integral types are supported by popcount");     \
  static_assert(sizeof(T) <= sizeof(ULL),                             \
                "size bigger than unsigned long long not supported"); \
  static_assert(!std::is_same<T, bool>::value, "bool type not supported")

template <typename T>
inline unsigned Popcount(T val) {
  CHECK_TYPE(T);
  return sizeof(T) <= sizeof(U)
             ? __builtin_popcount(val)
             : sizeof(T) <= sizeof(UL)
                   ? __builtin_popcountl(val)
                   : sizeof(T) <= sizeof(ULL) ? __builtin_popcountll(val) : 0;
}

template <typename T>
inline unsigned Ctz(T val) {
  CHECK_TYPE(T);
  return sizeof(T) <= sizeof(U)
             ? __builtin_ctz(val)
             : sizeof(T) <= sizeof(UL)
                   ? __builtin_ctzl(val)
                   : sizeof(T) <= sizeof(ULL) ? __builtin_ctzll(val) : 0;
}

template <typename Tx, typename Ty>
inline Tx Pdep(Tx x, Ty y) {
  CHECK_TYPE(Tx);
  return sizeof(Tx) <= sizeof(U)
             ? _pdep_u32(x, y)
             : sizeof(Tx) <= sizeof(ULL) ? _pdep_u64(x, y) : 0;
}

#undef CHECK_TYPE
}  // namespace internal

enum rnsmode { kRankOnly, kRandAndSelect };

template <unsigned BaseSize, unsigned AlphabetSize,
          unsigned Mode = rnsmode::kRandAndSelect, typename TWord = uint64_t,
          typename TInterval = uint32_t,
          unsigned BasePerL2Interval = 65536 / BaseSize,
          unsigned BasePerL1Interval = 1024 / BaseSize,
          unsigned SelectSampleSize = 4096>
class RankAndSelect {
 public:
  static const unsigned kBitsPerByte = 8;
  static const unsigned kBitsPerWord = sizeof(TWord) * kBitsPerByte;
  static const unsigned kBitsPerBase = BaseSize;
  static const unsigned kAlphabetSize = AlphabetSize;
  static const unsigned kBasesPerL1 = BasePerL1Interval;
  static const unsigned kBasesPerL2 = BasePerL2Interval;
  static const unsigned kSelectSampleSize = SelectSampleSize;
  static const unsigned kL1PerL2 = kBasesPerL2 / kBasesPerL1;
  static const unsigned kBasesPerWord = kBitsPerWord / kBitsPerBase;
  using size_type = int64_t;
  using word_type = TWord;
  static const size_type kNullID = static_cast<size_type>(-1);

  RankAndSelect() {
    for (unsigned i = kBitsPerBase == 1 ? 1 : 0; i < kAlphabetSize; ++i) {
      xor_masks_[i] = 0;
      for (unsigned j = 0; j < kBasesPerWord; ++j) {
        xor_masks_[i] |= (word_type)i << (kBitsPerBase * j);
      }
      xor_masks_[i] = ~xor_masks_[i];
    }
  }

  ~RankAndSelect() = default;

  void from_packed_array(const word_type *packed_array, size_type size) {
    size_type num_l1 = DivCeiling(size, kBasesPerL1) + 1;
    size_type num_l2 = DivCeiling(size, kBasesPerL2) + 1;

    for (unsigned c = BaseSize == 1 ? 1 : 0; c < kAlphabetSize; ++c) {
      const word_type *cur_word = packed_array;
      size_type count = 0;
      l2_occ_[c] = std::vector<size_type>(num_l2);
      l1_occ_[c] = std::vector<uint16_t>(num_l1);

      size_type size_rd = size - size % kBasesPerL1;
      for (size_type i = 0; i < size_rd;
           i += kBasesPerL1, cur_word += kBasesPerL1 / kBasesPerWord) {
        if (i % kBasesPerL2 == 0) {
          l2_occ_[c][i / kBasesPerL2] = count;
        }
        l1_occ_[c][i / kBasesPerL1] = count - l2_occ_[c][i / kBasesPerL2];
        count += CountCharInWords(c, cur_word, kBasesPerL1 / kBasesPerWord);
      }
      for (size_type i = size_rd; i < size; i += kBasesPerWord, ++cur_word) {
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
        for (size_type i = 0; i < num_l1; ++i) {
          while (static_cast<size_type>(rank2itv_[c].size() *
                                        kSelectSampleSize) < OccValue(c, i)) {
            rank2itv_[c].push_back(i - 1);
          }
        }
        rank2itv_[c].push_back(num_l1 - 1);
      }
    }
    packed_array_ = packed_array;
    this->size_ = size;
  }

  size_type rank(size_type pos) const {
    static_assert(BaseSize == 1, "");
    return InternalRank(1, pos);
  }

  size_type rank(uint8_t c, size_type pos) const {
    static_assert(BaseSize != 1, "");
    return InternalRank(c, pos);
  }

  size_type select(size_type ranking) const {
    static_assert(BaseSize == 1, "");
    return InternalSelect(1, ranking);
  }

  size_type select(uint8_t c, size_type ranking) const {
    static_assert(BaseSize != 1, "");
    return InternalSelect(c, ranking);
  }

  size_type pred(uint8_t c, size_type pos) const {
    // the last c in [0...pos]
    if (GetBaseAt(pos) == c) {
      return pos;
    }
    return InternalSelect(c, InternalRank(c, pos) - 1);
  }

  size_type pred(size_type pos) const {
    static_assert(BaseSize == 1, "");
    return pred(1, pos);
  }

  size_type pred(uint8_t c, size_type pos, int step) const {
    // the last c in [pos-step, pos], return pos-step-1 if not exist
    size_type end = pos >= step ? pos - step : 0;
    while (pos >= end) {
      if (GetBaseAt(pos) == c) {
        return pos;
      }
      --pos;
    }
    return pos;
  }

  size_type succ(uint8_t c, size_type pos) const {
    // the first c in [pos...ReadLength]
    if (GetBaseAt(pos) == c) {
      return pos;
    }
    return InternalSelect(c, InternalRank(c, pos - 1));
  }

  size_type succ(size_type pos) const { return succ(1, pos); }

  size_type succ(uint8_t c, size_type pos, int step) const {
    // the first c in [pos, pos+step], return pos+step+1 if not exist
    size_type end = pos + step;
    if (end >= size_) end = size_;
    while (pos <= end) {
      if (GetBaseAt(pos) == c) {
        return pos;
      }
      ++pos;
    }
    return pos;
  }

 private:
  unsigned CountCharInWord(uint8_t c, word_type x,
                           word_type mask = word_type(-1)) const {
    if (BaseSize != 1) {
      x ^= xor_masks_[c];
      x = internal::PackToLowestBit<BaseSize>(x);
      x &= kPopcntMask;
    }
    return internal::Popcount(x & mask);
  }

  unsigned CountCharInWords(uint8_t c, const word_type *ptr,
                            unsigned n_words) const {
    unsigned count = 0;
    for (unsigned i = 0; i < n_words; ++i) {
      count += CountCharInWord(c, ptr[i]);
    }
    return count;
  }

  unsigned SelectInWord(uint8_t c, int num_c, word_type x) const {
    if (BaseSize != 1) {
      x ^= xor_masks_[c];
      x = internal::PackToLowestBit<BaseSize>(x);
      x &= kPopcntMask;
    }
#if defined(__BMI2__) && defined(USE_BMI2)
    return internal::Ctz(internal::Pdep(word_type(1) << (num_c - 1), x)) /
           kBitsPerBase;
#else
    unsigned trailing_zeros = 0;
    while (num_c > 0) {
      trailing_zeros = internal::Ctz(x);
      x ^= TWord{1} << trailing_zeros;
      --num_c;
    }
    return trailing_zeros / kBitsPerBase;
#endif
  }

  void PrefetchOcc(uint8_t c, int64_t i) const {
    __builtin_prefetch(&l2_occ_[c][i / kL1PerL2], 0);
    __builtin_prefetch(&l1_occ_[c][i], 0);
  }

  size_type OccValue(uint8_t c, size_type i) const {
    return l2_occ_[c][i / kL1PerL2] + l1_occ_[c][i];
  }

 private:
  size_type InternalRank(uint8_t c, size_type pos) const {
    // the number of c's in [0...pos]
    if (pos >= size_)
      return kNullID;
    else if (pos == size_ - 1)
      return char_count_[c];

    ++pos;
    size_type itv_idx = (pos + kBasesPerL1 / 2 - 1) / kBasesPerL1;
    size_type sampled_pos = itv_idx * kBasesPerL1;

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

  size_type InternalSelect(uint8_t c, size_type k) const {
    static_assert(Mode != rnsmode::kRankOnly,
                  "cannot select on rank only struct");
    // return the pos (0-based) of the kth (0-based) c
    if (k > char_count_[c])
      return kNullID;
    else if (k == char_count_[c])
      return size_;
    // first locate which interval Select(c, k) falls
    size_type interval_l = rank2itv_[c][k / kSelectSampleSize];
    size_type interval_r = rank2itv_[c][DivCeiling(k, kSelectSampleSize)];
    PrefetchOcc(c, interval_l);
    while (interval_r > interval_l) {
      size_type interval_m = (interval_r + interval_l + 1) / 2;
      if (OccValue(c, interval_m) > k) {
        interval_r = interval_m - 1;
      } else {
        interval_l = interval_m;
      }
    }
    // refined select
    __builtin_prefetch(packed_array_ +
                       interval_l * kBasesPerL1 / kBasesPerWord);
    unsigned remain = k + 1 - OccValue(c, interval_l);
    unsigned exceed = (interval_l + 1) * kBasesPerL1 >= size_
                          ? kBasesPerL1
                          : (OccValue(c, interval_l + 1) - (k + 1));
    if (remain <= exceed * 2) {
      return SelectFwd(c, interval_l, remain);
    } else {
      return SelectBwd(c, interval_l, exceed);
    }
  }

  size_type RankFwd(uint8_t c, TInterval itv, size_type sampled_pos,
                    unsigned n_bases) const {
    unsigned n_words = n_bases / kBasesPerWord;
    const word_type *p =
        packed_array_ + sampled_pos / kBasesPerWord - n_words - 1;
    __builtin_prefetch(p);
    unsigned n_residual = n_bases % kBasesPerWord;
    unsigned count = 0;
    if (n_residual != 0) {
      word_type mask =
          1 + ~(1ULL << kBitsPerBase * (kBasesPerWord - n_residual));
      count += CountCharInWord(c, p[0], mask);
    }
    count += CountCharInWords(c, p + 1, n_words);
    return OccValue(c, itv) - count;
  }

  size_type RankBwd(uint8_t c, TInterval itv, size_type sampled_pos,
                    unsigned n_bases) const {
    const word_type *p = packed_array_ + sampled_pos / kBasesPerWord;
    __builtin_prefetch(p);
    unsigned n_words = n_bases / kBasesPerWord;
    unsigned n_residual = n_bases % kBasesPerWord;
    unsigned count = 0;
    count += CountCharInWords(c, p, n_words);
    if (n_residual != 0) {
      word_type mask = (1ULL << kBitsPerBase * n_residual) - 1;
      count += CountCharInWord(c, p[n_words], mask);
    }
    return OccValue(c, itv) + count;
  }

  size_type SelectFwd(uint8_t c, TInterval itv, unsigned remain) const {
    size_type pos = static_cast<size_type>(itv) * kBasesPerL1;
    const word_type *begin = packed_array_ + pos / kBasesPerWord;
    const word_type *p = begin;
    for (unsigned popcnt; (popcnt = CountCharInWord(c, *p)) < remain;
         remain -= popcnt, ++p)
      ;
    return pos + (p - begin) * kBasesPerWord + SelectInWord(c, remain, *p);
  }

  size_type SelectBwd(uint8_t c, TInterval itv_l, unsigned exceed) const {
    size_type pos = (static_cast<size_type>(itv_l) + 1) * kBasesPerL1;
    const word_type *end = packed_array_ + pos / kBasesPerWord - 1;
    const word_type *p = end;
    unsigned popcnt;
    for (; (popcnt = CountCharInWord(c, *p)) <= exceed; exceed -= popcnt, --p)
      ;
    return pos - kBasesPerWord * (end - p) -
           (kBasesPerWord - SelectInWord(c, popcnt - exceed, *p));
  }

  uint8_t GetBaseAt(size_type i) const {
    return (packed_array_[i / kBasesPerWord] >>
            (i % kBasesPerWord * kBitsPerBase)) &
           ((1 << kBitsPerBase) - 1);
  }

  template <typename T1, typename T2>
  static T1 DivCeiling(T1 x, T2 y) {
    return (x + y - 1) / y;
  };

  static const word_type kPopcntMask =
      internal::PopcountMask<BaseSize, word_type>::value;
  static const uint64_t kPopcntMask64 =
      internal::PopcountMask<BaseSize, uint64_t>::value;
  size_type size_{};
  size_type char_count_[kAlphabetSize]{};
  // main memory for the structure
  const word_type *packed_array_;
  // sampled structure for rank and select
  // two level sampling for rank (occ value)
  // call the function OccValue(c, i) to get the number of c's
  // in packed_array_[0...i*kBasesPerL1-1]
  // sampling for select
  // rank_interval_lookup_[c][i]=j: the jth interval (0 based)
  // contains the (i*kSelectSampleSize)th (0 based) c
  // i.e. OccValue(c, j)<=i*kSelectSampleSize and OccValue(c,
  // j+1)>i*kSelectSampleSize
  std::vector<TInterval> rank2itv_[kAlphabetSize];
  std::vector<uint16_t> l1_occ_[kAlphabetSize];   // level 1 OCC
  std::vector<size_type> l2_occ_[kAlphabetSize];  // level 2 OCC
  word_type xor_masks_[kAlphabetSize];
  // e.g. if c = 0110(2), popcount_xorers_[c] = 1001 1001 1001 1001...(2),
  // to make all c's in a word 1111
  static_assert((1ull << kBitsPerBase) >= kAlphabetSize, "");
  static_assert(kBitsPerWord % kBitsPerBase == 0, "");
  static_assert(kBitsPerBase <= 8, "");
  static_assert(kBasesPerL2 <= 65536, "");
};

}  // namespace kmlib

#endif  // KMLIB_RNS_H
