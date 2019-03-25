//
// Created by vout on 11/4/18.
//

#include "kmcompactvector.h"

#include <iostream>

template <unsigned BaseSize>
void print_v(const kmlib::CompactVector<BaseSize> &v) {
  std::cout << v.size() << ", " << v.capacity() << ": ";
  for (size_t i = 0; i < v.size(); ++i) {
    std::cout << v[i] << ' ';
  }
  std::cout << std::endl;
}

int main() {
  kmlib::CompactVector<2> v;
  for (size_t i = 0; i < 17; ++i) {
    v.push_back(i % 4);
    print_v(v);
    v[i] = 3 - v[i];
    print_v(v);
  }
  while (v.size() > 0) {
    v.pop_back();
    print_v(v);
  }
  kmlib::CompactVector<4> v4;
  for (size_t i = 0; i < 17; ++i) {
    v4.push_back(i % 16);
    print_v(v4);
    v4[i] = 15 - v4[i];
    print_v(v4);
  }
  while (v4.size() > 0) {
    v4.pop_back();
    print_v(v4);
  }
  return 0;
}