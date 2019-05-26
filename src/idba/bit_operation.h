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
#include "kmlib/kmbit.h"

namespace bit_operation {

template <typename T>
inline void ReverseComplement(T &value) {
  value = kmlib::bit::ReverseComplement<2, T>(value);
}

inline int BitToIndex(uint8_t x) {
  const static int bit_to_index[] = {
      0, 0, 1, 0, 2, 0, 0, 0, 3,
  };
  return bit_to_index[x];
}

}  // namespace bit_operation

#endif
