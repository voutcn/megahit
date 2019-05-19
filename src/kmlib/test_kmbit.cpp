//
// Created by vout on 11/9/18.
//

#include <bitset>
#include <cassert>
#include <iostream>
#include "kmbit.h"
#include "kmcompactvector.h"

int main() {
  uint32_t a = 24128923;
  uint64_t b = 131231231123123ULL;
  auto rev_a1 = kmlib::bit::Reverse<1>(a);
  auto rev_a = kmlib::bit::Reverse<2>(a);
  auto rev_b = kmlib::bit::Reverse<4>(b);

  std::cout << std::bitset<32>(a) << std::endl;
  std::cout << std::bitset<32>(rev_a1) << std::endl;

  for (size_t i = 0; i < 32; ++i) {
    auto b1 = kmlib::CompactVector<1, uint32_t>::at(&a, i);
    auto b2 = kmlib::CompactVector<1, uint32_t>::at(&rev_a1, 31 - i);
    assert(b1 == b2);
  }

  for (size_t i = 0; i < 16; ++i) {
    auto b1 = kmlib::CompactVector<2, uint32_t>::at(&a, i);
    auto b2 = kmlib::CompactVector<2, uint32_t>::at(&rev_a, 15 - i);
    assert(b1 == b2);
  }

  for (size_t i = 0; i < 16; ++i) {
    auto b1 = kmlib::CompactVector<4, uint64_t>::at(&b, i);
    auto b2 = kmlib::CompactVector<4, uint64_t>::at(&rev_b, 15 - i);
    assert(b1 == b2);
  }
  return 0;
}