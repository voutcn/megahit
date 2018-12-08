//
// Created by vout on 12/8/18.
//

#include "sorting.h"

#include <cassert>
#include "kmlib/kmsort.h"

namespace {
template <int NWords, int NExtraWord>
struct Substr {
  uint32_t data[NWords + NExtraWord];
  static const int n_bytes = sizeof(data[0]) * NWords - 2;  // prefix len = 16 bits = 2 bytes
  bool operator< (const Substr &rhs) const {
    for (int i = 0; i < NWords; ++i) {
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

template <int NWords, int NExtraWords>
inline void SortSubstrHelper(uint32_t *substr, int words_per_substr, int64_t n, int extra_words) {
  assert(words_per_substr > 0 && words_per_substr <= NWords);
  assert(extra_words >= 0 && extra_words <= NExtraWords);
  if (words_per_substr < NWords) {
    SortSubstrHelper<(NWords - 1 > 1 ? NWords - 1 : 1), NExtraWords>(substr, words_per_substr, n, extra_words);
    return;
  }
  if (extra_words < NExtraWords) {
    SortSubstrHelper<NWords, (NExtraWords - 1 > 0 ? NExtraWords - 1: 0)>(substr, words_per_substr, n, extra_words);
    return;
  }
  auto ptr = reinterpret_cast<Substr<NWords, NExtraWords>*>(substr);
  kmlib::kmsort(ptr, ptr + n);
}
}

void SortSubStr(uint32_t *substr, int words_per_substr, int64_t n, int extra_words) {
  SortSubstrHelper<12, 2>(substr, words_per_substr, n, extra_words);
}