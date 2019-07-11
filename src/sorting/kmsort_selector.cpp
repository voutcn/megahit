//
// Created by vout on 12/8/18.
//

#include "kmsort_selector.h"

#include <definitions.h>
#include <utils/utils.h>
#include <cassert>
#include "kmlib/kmsort.h"

namespace {
template <int NWords, int NExtraWord>
struct Substr {
  uint32_t data[NWords + NExtraWord];
  static const int n_bytes =
      sizeof(data[0]) * NWords - 2;  // prefix len = 16 bits = 2 bytes
  bool operator<(const Substr &rhs) const {
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
    return data[NWords - 1 - k / sizeof(uint32_t)] >>
               (k % sizeof(uint32_t) * 8) &
           0xFF;
  }
} __attribute((packed));

constexpr int kMaxWords =
    (kMaxK * kBitsPerEdgeChar + 3 + 1 + kBitsPerMul + kBitsPerEdgeWord - 1) /
    kBitsPerEdgeWord;

template <int NWords, int NExtraWords>
std::function<void(uint32_t *, int64_t)> SelectSortingFuncHelper(
    int words_per_substr, int extra_words) {
  assert(words_per_substr > 0 && words_per_substr <= NWords);
  assert(extra_words >= 0 && extra_words <= NExtraWords);
  if (words_per_substr < NWords) {
    return SelectSortingFuncHelper<(NWords - 1 > 1 ? NWords - 1 : 1),
                                   NExtraWords>(words_per_substr, extra_words);
  }
  if (extra_words < NExtraWords) {
    return SelectSortingFuncHelper<NWords,
                                   (NExtraWords - 1 > 0 ? NExtraWords - 1 : 0)>(
        words_per_substr, extra_words);
  }
  return [](uint32_t *substr_ptr, int64_t n) {
    auto ptr = reinterpret_cast<Substr<NWords, NExtraWords> *>(substr_ptr);
    kmlib::kmsort(ptr, ptr + n);
  };
}

}  // namespace

std::function<void(uint32_t *, int64_t)> SelectSortingFunc(int words_per_substr,
                                                           int extra_words) {
  return SelectSortingFuncHelper<kMaxWords, 2>(words_per_substr, extra_words);
}