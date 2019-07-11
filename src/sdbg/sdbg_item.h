//
// Created by vout on 11/5/18.
//

#ifndef MEGAHIT_SDBG_ITEM_H
#define MEGAHIT_SDBG_ITEM_H

#include <cstdint>
#include "sdbg_def.h"

/**
 * Pack w, last, tip and multiplicity (small) into one struct of 16 bits
 */
struct SdbgItem {
  SdbgItem() = default;
  SdbgItem(uint8_t w, uint8_t last, uint8_t tip, small_mul_t mul)
      : w(w), last(last), tip(tip), mul(mul) {}
  uint8_t w : kBitsPerWChar;
  uint8_t last : 1;
  uint8_t tip : 1;
  uint8_t : 6 - kBitsPerWChar;
  small_mul_t mul : sizeof(small_mul_t) * 8;
  static_assert(sizeof(small_mul_t) <= 2, "");
};

#endif  // MEGAHIT_SDBG_ITEM_H
