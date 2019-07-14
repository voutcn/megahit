//
// Created by vout on 11/5/18.
//

#ifndef MEGAHIT_SDBG_DEF_H
#define MEGAHIT_SDBG_DEF_H

#include <cstdint>
#include <limits>

using mul_t = uint16_t;
static const int kBitsPerMul = sizeof(mul_t) * 8;
static const int kMaxMul = std::numeric_limits<mul_t>::max();

using small_mul_t = uint8_t;
static const small_mul_t kMaxSmallMul =
    std::numeric_limits<small_mul_t>::max() - 1;
static const small_mul_t kSmallMulSentinel =
    std::numeric_limits<small_mul_t>::max();

static const unsigned kMaxK = 255;
static const unsigned kAlphabetSize = 4;
static const unsigned kBitsPerChar = 2;
static const unsigned kWAlphabetSize = 9;
static const unsigned kBitsPerWChar = 4;

using label_word_t = uint32_t;
static const unsigned kCharsPerLabelWord =
    sizeof(label_word_t) * 8 / kBitsPerChar;

#endif  // MEGAHIT_SDBG_DEF_H
