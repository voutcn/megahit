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

template <unsigned NWords>
struct Substr {
  uint32_t data[NWords];
  static const int n_bytes = sizeof(data) - 2;  // prefix len = 16 bits = 2 bytes
  bool operator< (const Substr &rhs) const {
      for (unsigned i = 0; i < NWords; ++i) {
          if (data[i] < rhs.data[i]) {
              return true;
          } else if (data[i] > rhs.data[i]) {
              return false;
          }
      }
      return false;
  }
  int kth_byte(int k) const {
      return data[NWords - 1 - k / sizeof(uint32_t)] >> (k % sizeof(uint32_t) * 8) & 0xFF;
  }
} __attribute((packed));

template <unsigned NWords>
inline bool SortSubstrHelper(uint32_t *substr, int words_per_substr, int64_t n) {
    if (words_per_substr != NWords) return false;
    auto ptr = reinterpret_cast<Substr<NWords>*>(substr);
    kmlib::kmsort(ptr, ptr + n);
    return true;
}

inline void SortSubStr(uint32_t *substr, int words_per_substr, int64_t n) {
    if (SortSubstrHelper<1>(substr, words_per_substr, n)) return;
    if (SortSubstrHelper<2>(substr, words_per_substr, n)) return;
    if (SortSubstrHelper<3>(substr, words_per_substr, n)) return;
    if (SortSubstrHelper<4>(substr, words_per_substr, n)) return;
    if (SortSubstrHelper<5>(substr, words_per_substr, n)) return;
    if (SortSubstrHelper<6>(substr, words_per_substr, n)) return;
    if (SortSubstrHelper<7>(substr, words_per_substr, n)) return;
    if (SortSubstrHelper<8>(substr, words_per_substr, n)) return;
    if (SortSubstrHelper<9>(substr, words_per_substr, n)) return;
    if (SortSubstrHelper<10>(substr, words_per_substr, n)) return;
    if (SortSubstrHelper<11>(substr, words_per_substr, n)) return;
    if (SortSubstrHelper<12>(substr, words_per_substr, n)) return;
    assert(false);
}