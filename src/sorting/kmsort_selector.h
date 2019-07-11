#ifndef MEGAHIT_KMSORT_SELECTOR_H
#define MEGAHIT_KMSORT_SELECTOR_H

#include <cstdint>
#include <functional>

std::function<void(uint32_t*, int64_t)> SelectSortingFunc(int words_per_substr,
                                                          int extra_words);

#endif