//
// Created by vout on 7/3/2018.
//

#ifndef KMLIB_KMSORT_H
#define KMLIB_KMSORT_H

#include <algorithm>
#include <type_traits>

namespace kmlib {

namespace kmsortconst {
static const int kRadixBits = 8;
static const int kRadixMask = (1 << kRadixBits) - 1;
static const int kNumBins = 1 << kRadixBits;
static const int kInsertSortThreshold = 64;
}  // namespace kmsortconst

namespace internal {

template <class RandomIt, class ValueType, class RadixTraits>
inline void insert_sort_core(RandomIt s, RandomIt e, RadixTraits rt) {
  for (RandomIt i = s + 1; i != e; ++i) {
    if (rt(*i, *(i - 1))) {
      RandomIt j;
      ValueType tmp = *i;
      *i = *(i - 1);
      for (j = i - 1; j > s && rt(tmp, *(j - 1)); --j) {
        *j = *(j - 1);
      }
      *j = tmp;
    }
  }
}

template <class T, class RadixTrait, int ByteIndex, int ByteIndexEnd>
inline int kth_byte(const T &val, RadixTrait rt, int byte_index) {
  return ByteIndex >= ByteIndexEnd ? rt.kth_byte(val, ByteIndex)
                                   : rt.kth_byte(val, byte_index);
}

template <class RandomIt, class ValueType, class RadixTrait, int ByteIndex,
          int ByteIndexEnd>
inline void radix_sort_core(RandomIt s, RandomIt e, RadixTrait rt,
                            int byte_index) {
  RandomIt last_[kmsortconst::kNumBins + 1];
  RandomIt *last = last_ + 1;
  size_t count[kmsortconst::kNumBins] = {0};

  for (RandomIt i = s; i < e; ++i) {
    ++count[kth_byte<ValueType, RadixTrait, ByteIndex, ByteIndexEnd>(
        *i, rt, byte_index)];
  }

  last_[0] = last_[1] = s;

  for (int i = 1; i < kmsortconst::kNumBins; ++i) {
    last[i] = last[i - 1] + count[i - 1];
  }

  for (int i = 0; i < kmsortconst::kNumBins; ++i) {
    RandomIt end = last[i - 1] + count[i];
    if (end == e) {
      last[i] = e;
      break;
    }
    while (last[i] != end) {
      ValueType swapper = *last[i];
      int tag = kth_byte<ValueType, RadixTrait, ByteIndex, ByteIndexEnd>(
          swapper, rt, byte_index);
      if (tag != i) {
        do {
          std::swap(swapper, *last[tag]++);
        } while (
            (tag = kth_byte<ValueType, RadixTrait, ByteIndex, ByteIndexEnd>(
                 swapper, rt, byte_index)) != i);
        *last[i] = swapper;
      }
      ++last[i];
    }
  }

  if (ByteIndex > ByteIndexEnd) {
    const int kNextIndex = ByteIndex > 0 ? ByteIndex - 1 : 0;
    for (int i = 0; i < kmsortconst::kNumBins; ++i) {
      if (count[i] > kmsortconst::kInsertSortThreshold) {
        radix_sort_core<RandomIt, ValueType, RadixTrait, kNextIndex,
                        ByteIndexEnd>(last[i - 1], last[i], rt, kNextIndex);
      } else if (count[i] > 1) {
        insert_sort_core<RandomIt, ValueType, RadixTrait>(last[i - 1], last[i],
                                                          rt);
      }
    }
  } else if (byte_index > 0) {
    for (int i = 0; i < kmsortconst::kNumBins; ++i) {
      if (count[i] > kmsortconst::kInsertSortThreshold) {
        radix_sort_core<RandomIt, ValueType, RadixTrait, 0, ByteIndexEnd>(
            last[i - 1], last[i], rt, byte_index - 1);
      } else if (count[i] > 1) {
        insert_sort_core<RandomIt, ValueType, RadixTrait>(last[i - 1], last[i],
                                                          rt);
      }
    }
  }
}

template <class RandomIt, class ValueType, class RadixTraits>
inline void radix_sort_entry(RandomIt s, RandomIt e, ValueType *,
                             RadixTraits radix_traits) {
  if (std::distance(s, e) <= 1) {
    return;
  } else if (std::distance(s, e) <= kmsortconst::kInsertSortThreshold) {
    insert_sort_core<RandomIt, ValueType, RadixTraits>(s, e, radix_traits);
  } else {
    const int kByteIndexEnd =
        RadixTraits::n_bytes > 8 ? RadixTraits::n_bytes - 8 : 0;
    radix_sort_core<RandomIt, ValueType, RadixTraits, RadixTraits::n_bytes - 1,
                    kByteIndexEnd>(s, e, radix_traits,
                                   RadixTraits::n_bytes - 1);
  }
}

}  // namespace internal

/*!
 * @brief the default radix trait for kmsort
 * requires the class T implements operator<(), kth_byte() and n_bytes interface
 * @tparam T type of the items to be sorted
 */
template <class T, bool = std::is_integral<T>::value>
struct RadixTraits {
  static const int n_bytes = T::n_bytes;
  int kth_byte(const T &x, int k) { return x.kth_byte(k); }
  bool operator()(const T &x, const T &y) const { return x < y; }
};

/*!
 * @brief the default radix trait for integers
 * @tparam T the type of integer
 */
template <class T>
struct RadixTraits<T, true> {
  static const int n_bytes = sizeof(T);
  static const T kMSBMask = T(T(-1) < 0) << (sizeof(T) * 8 - 1);
  int kth_byte(const T &x, int k) {
    return ((x ^ kMSBMask) >> (kmsortconst::kRadixBits * k)) &
           kmsortconst::kRadixMask;
  }
  bool operator()(const T &x, const T &y) const { return x < y; }
};

/*!
 * @brief sort [s, e) using explicitly defined radix traits
 * @tparam RandomIt the type of random iterator
 * @tparam RadixTraits the type of radix traits
 * @param s the first position
 * @param e the last position
 * @param radix_traits explicitly defined radix traits
 */
template <class RandomIt, class RadixTraits>
inline void kmsort(RandomIt s, RandomIt e, RadixTraits radix_traits) {
  typename std::iterator_traits<RandomIt>::value_type *_(0);
  internal::radix_sort_entry(s, e, _, radix_traits);
}

/*!
 * @brief sort [s, e) using default RadixTraits
 * @tparam RandomIt the type of random iterator
 * @param s the first position
 * @param e the last position
 */
template <class RandomIt>
inline void kmsort(RandomIt s, RandomIt e) {
  using ValueType = typename std::iterator_traits<RandomIt>::value_type;
  kmsort(s, e, RadixTraits<ValueType>());
}

}  // namespace kmlib

#endif  // KMLIB_KMSORT_H