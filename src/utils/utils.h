//
// Created by dinghua.li on 1/27/18.
//

#ifndef MEGAHIT_UTILS_H
#define MEGAHIT_UTILS_H

#include <cassert>
#include <iostream>

template<typename T1, typename T2>
inline T1 DivCeiling(T1 x, T2 y) {
  return (x + y - 1) / y;
};

#define MERR std::cerr
#define MINFO std::cerr
#define MDEBUG std::cerr
#define MWARN std::cerr
#define MASSERT assert

#endif //MEGAHIT_UTILS_H
