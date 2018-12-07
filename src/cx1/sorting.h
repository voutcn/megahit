/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong & L3 Bioinformatics Limited
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* contact: Dinghua Li <dhli@cs.hku.hk> */

#include <assert.h>
#include "kmlib/kmsort.h"

template <int NWords>
struct Helper {
  uint32_t data[NWords];
  static const int n_bytes = sizeof(data);
  bool operator< (const Helper &rhs) const {
      for (int i = 0; i < NWords; ++i) {
          if (data[i] < rhs.data[i]) {
              return true;
          } else if (data[i] > rhs.data[i]) {
              return false;
          }
      }
      return true;
  }
  int kth_byte(int k) const {
      return data[NWords - 1 - k / sizeof(uint32_t)] >> (k % sizeof(uint32_t) * 8) & 0xFF;
  }
} __attribute((packed));

template <int NWords>
inline bool helper(uint32_t *substr, int words_per_substr, int64_t n) {
    if (words_per_substr != NWords) return false;
    auto ptr = reinterpret_cast<Helper<NWords>*>(substr);
    kmlib::kmsort(ptr, ptr + n);
    return true;
}

inline void lv2_cpu_radix_sort_st2(uint32_t *substr, int words_per_substr, int64_t n) {
    if (helper<1>(substr, words_per_substr, n)) return;
    if (helper<2>(substr, words_per_substr, n)) return;
    if (helper<3>(substr, words_per_substr, n)) return;
    if (helper<4>(substr, words_per_substr, n)) return;
    if (helper<5>(substr, words_per_substr, n)) return;
    if (helper<6>(substr, words_per_substr, n)) return;
    if (helper<7>(substr, words_per_substr, n)) return;
    if (helper<8>(substr, words_per_substr, n)) return;
    if (helper<9>(substr, words_per_substr, n)) return;
    if (helper<10>(substr, words_per_substr, n)) return;
    if (helper<11>(substr, words_per_substr, n)) return;
    if (helper<12>(substr, words_per_substr, n)) return;
    assert(false);
}