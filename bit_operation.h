/**
 * @file bit_operation.h
 * @brief Bit Operations.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-02
 * @modified by Dinghua Li
 * @date 2014-10-02
 */

#ifndef __BASIC_BIT_OPERATION_H_

#define __BASIC_BIT_OPERATION_H_

#include <stdint.h>


namespace bit_operation {

const uint64_t kSwap64Mask32 = 0x00000000FFFFFFFFULL;
const uint64_t kSwap64Mask16 = 0x0000FFFF0000FFFFULL;
const uint64_t kSwap64Mask8 = 0x00FF00FF00FF00FFULL;
const uint64_t kSwap64Mask4 = 0x0F0F0F0F0F0F0F0FULL;
const uint64_t kSwap64Mask2 = 0x3333333333333333ULL;
const uint64_t kSwap64Mask1 = 0x5555555555555555ULL;

const uint8_t kSwap8Mask4 = 0x0FU;
const uint8_t kSwap8Mask2 = 0x33U;
const uint8_t kSwap8Mask1 = 0x55U;

inline void ReverseComplement(uint64_t &value) {
    value = ((value & kSwap64Mask32) << 32) | ((value & ~kSwap64Mask32) >> 32);
    value = ((value & kSwap64Mask16) << 16) | ((value & ~kSwap64Mask16) >> 16);
    value = ((value & kSwap64Mask8) << 8) | ((value & ~kSwap64Mask8) >> 8);
    value = ((value & kSwap64Mask4) << 4) | ((value & ~kSwap64Mask4) >> 4);
    value = ((value & kSwap64Mask2) << 2) | ((value & ~kSwap64Mask2) >> 2);
    value = ~value;
}

inline void ReverseComplement(uint32_t &value) {
    value = ((value & kSwap64Mask16) << 16) | ((value & ~kSwap64Mask16) >> 16);
    value = ((value & kSwap64Mask8) << 8) | ((value & ~kSwap64Mask8) >> 8);
    value = ((value & kSwap64Mask4) << 4) | ((value & ~kSwap64Mask4) >> 4);
    value = ((value & kSwap64Mask2) << 2) | ((value & ~kSwap64Mask2) >> 2);
    value = ~value;
}

inline void ReverseComplement(uint16_t &value) {
    value = ((value & kSwap64Mask8) << 8) | ((value & ~kSwap64Mask8) >> 8);
    value = ((value & kSwap64Mask4) << 4) | ((value & ~kSwap64Mask4) >> 4);
    value = ((value & kSwap64Mask2) << 2) | ((value & ~kSwap64Mask2) >> 2);
    value = ~value;
}

inline void ReverseComplement(uint8_t &value) {
    value = ((value & kSwap64Mask4) << 4) | ((value & ~kSwap64Mask4) >> 4);
    value = ((value & kSwap64Mask2) << 2) | ((value & ~kSwap64Mask2) >> 2);
    value = ~value;
}

inline uint8_t ReverseComplement(uint8_t value) {
    value = ((value & kSwap8Mask4) << 4) | ((value & ~kSwap8Mask4) >> 4);
    value = ((value & kSwap8Mask2) << 2) | ((value & ~kSwap8Mask2) >> 2);
    return ~value;
}

inline uint32_t Reverse(uint32_t value) {
    value = ((value & kSwap64Mask16) << 16) | ((value & ~kSwap64Mask16) >> 16);
    value = ((value & kSwap64Mask8) << 8) | ((value & ~kSwap64Mask8) >> 8);
    value = ((value & kSwap64Mask4) << 4) | ((value & ~kSwap64Mask4) >> 4);
    value = ((value & kSwap64Mask2) << 2) | ((value & ~kSwap64Mask2) >> 2);
    return value;
}

inline uint8_t Reverse(uint8_t value) {
    value = ((value & kSwap8Mask4) << 4) | ((value & ~kSwap8Mask4) >> 4);
    value = ((value & kSwap8Mask2) << 2) | ((value & ~kSwap8Mask2) >> 2);
    return value;
}

inline int BitCount(uint8_t x) {
    x = (x & kSwap8Mask1) + ((x >> 1) & kSwap8Mask1);
    x = (x & kSwap8Mask2) + ((x >> 2) & kSwap8Mask2);
    x = (x & kSwap8Mask4) + ((x >> 4) & kSwap8Mask4);
    return x;
}

inline int BitCount(uint64_t x) {
    x = (x & kSwap64Mask1) + ((x >> 1) & kSwap64Mask1);
    x = (x & kSwap64Mask2) + ((x >> 2) & kSwap64Mask2);
    x = (x & kSwap64Mask4) + ((x >> 4) & kSwap64Mask4);
    x = (x & kSwap64Mask8) + ((x >> 8) & kSwap64Mask8);
    x = (x & kSwap64Mask16) + ((x >> 16) & kSwap64Mask16);
    x = (x & kSwap64Mask32) + ((x >> 32) & kSwap64Mask32);
    return x;
}

inline int BaseCount(uint64_t x) {
    x = (x & kSwap64Mask1) | ((x >> 1) & kSwap64Mask1);
    x = (x & kSwap64Mask2) + ((x >> 2) & kSwap64Mask2);
    x = (x & kSwap64Mask4) + ((x >> 4) & kSwap64Mask4);
    x = (x & kSwap64Mask8) + ((x >> 8) & kSwap64Mask8);
    x = (x & kSwap64Mask16) + ((x >> 16) & kSwap64Mask16);
    x = (x & kSwap64Mask32) + ((x >> 32) & kSwap64Mask32);
    return x;
}

inline int BitToIndex(uint8_t x) {
    static int bit_to_index[] = {0, 0, 1, 0, 2, 0, 0, 0, 3,};
    return bit_to_index[x];
}

template<uint32_t number>
struct MostSignificantBit {
    static const int value = 1 + MostSignificantBit < (number >> 1) >::value;
};

template<>
struct MostSignificantBit<0> {
    static const int value = 0;
};


}

#endif

