//
// Created by vout on 27/2/2018.
//

#ifndef KMLIB_BIT_MAGIC_H
#define KMLIB_BIT_MAGIC_H

#include <type_traits>
namespace kmlib {
namespace bit {

using U = unsigned int;
using UL = unsigned long int;
using ULL = unsigned long long int;

#define KMLIB_CHECK_TYPE(T)                                           \
  static_assert(std::is_integral<T>::value,                           \
                "only integral types are supported by popcount");     \
  static_assert(sizeof(T) <= sizeof(ULL),                             \
                "size bigger than unsigned long long not supported"); \
  static_assert(!std::is_same<T, bool>::value, "bool type not supported")

namespace internal {

template <typename T, unsigned MaskSize, unsigned MaskIndex>
struct SubSwapMask {
 private:
  static const T prev =
      SubSwapMask<T, MaskSize, MaskIndex - MaskSize * 2>::value;

 public:
  static const T value = prev
                             ? (prev << MaskSize * 2) | ((T(1) << MaskSize) - 1)
                             : ((T(1) << MaskSize) - 1);
};

template <typename T, unsigned MaskSize>
struct SubSwapMask<T, MaskSize, 0> {
  static const T value = T(0);
};

template <typename T, unsigned MaskSize>
struct SwapMask {
  static const T value = SubSwapMask<T, MaskSize, sizeof(T) * 8>::value;
};

template <typename T, unsigned BaseSize, unsigned MaskSize>
struct SwapMaskedBitsHelper {
  T operator()(T value) {
    value = ((value & SwapMask<T, MaskSize>::value) << MaskSize) |
            ((value & ~SwapMask<T, MaskSize>::value) >> MaskSize);
    return SwapMaskedBitsHelper<T, BaseSize, MaskSize / 2>()(value);
  }
};

template <typename T, unsigned BaseSize>
struct SwapMaskedBitsHelper<T, BaseSize, BaseSize> {
  T operator()(T value) {
    return ((value & SwapMask<T, BaseSize>::value) << BaseSize) |
           ((value & ~SwapMask<T, BaseSize>::value) >> BaseSize);
  }
};

template <typename T, unsigned BaseSize, unsigned MaskSize>
inline T SwapMaskedBits(T value) {
  return SwapMaskedBitsHelper<T, BaseSize, MaskSize>()(value);
};

}  // namespace internal

/*!
 * @brief reverse an integer at a resolution base of BaseSize
 * @details e.g. for BaseSize = 2 and T = uint64_t, this function swap the 1st
 * base
 * with the 32th, the 2nd with the 31th and so on. The 1st base the first two
 * bits, and the last base is the last 2 bits
 * @tparam BaseSize the base size, i.e. the number of bits per base
 * @tparam T type of integer
 * @param value the value to reverse
 * @return the reversed value
 */
template <unsigned BaseSize, typename T>
inline T Reverse(T value) {
  KMLIB_CHECK_TYPE(T);
  static_assert(sizeof(T) * 8 % BaseSize == 0,
                "Reverse only support base size of power of 2");
  return internal::SwapMaskedBits<T, BaseSize, sizeof(T) * 4>(
      value);  // this will be optimized to bswap
}

/*!
 * @brief the same as ~Reverse(value)
 * @tparam BaseSize the base size, i.e. the number of bits per base
 * @tparam T the type of integer
 * @param value the value to reverse and then complement
 * @return the reverse-complemented value
 */
template <unsigned BaseSize, typename T>
inline T ReverseComplement(T value) {
  return ~Reverse<BaseSize>(value);
}

/*!
 * the value of floor(log2(x)) for a const integer x
 * @tparam number the number
 */
template <unsigned number>
struct FloorLog2 {
  static const int value = 1 + FloorLog2<(number >> 1)>::value;
};

template <>
struct FloorLog2<1> {
  static const int value = 0;
};

template <typename T>
inline unsigned Popcount(T val) {
  KMLIB_CHECK_TYPE(T);
  return sizeof(T) <= sizeof(U)
             ? __builtin_popcount(val)
             : sizeof(T) <= sizeof(UL)
                   ? __builtin_popcountl(val)
                   : sizeof(T) <= sizeof(ULL) ? __builtin_popcountll(val) : 0;
}

}  // namespace bit
}  // namespace kmlib

#endif  // KMLIB_BIT_MAGIC_H
