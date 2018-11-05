//
// Created by vout on 11/5/18.
//

#ifndef MEGAHIT_SDBG_DEF_H
#define MEGAHIT_SDBG_DEF_H

#include <cstdint>
#include <limits>

typedef uint16_t mul_t;
static const int kBitsPerMul = sizeof(mul_t) * 8;
static const int kMaxMul = std::numeric_limits<mul_t>::max();

typedef uint8_t small_mul_t;
static const small_mul_t kMaxSmallMul = std::numeric_limits<small_mul_t>::max() - 1;
static const small_mul_t kSmallMulSentinel = std::numeric_limits<small_mul_t>::max();

static const unsigned kMaxK = 255;
static const unsigned kAlphabetSize = 4;
static const unsigned kWAlphabetSize = 9;
static const unsigned kBitsPerChar = 2;

using LabelWordType = uint32_t;
static const unsigned kCharsPerLabelWord = sizeof(LabelWordType) * 8 / kBitsPerChar;

#endif //MEGAHIT_SDBG_DEF_H
