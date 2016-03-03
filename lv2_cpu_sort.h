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

#include <parallel/algorithm>
#include <algorithm>
#include <assert.h>
#include "definitions.h"

inline void lv2_cpu_sort(uint32_t *lv2_substrings, uint32_t *permutation, uint64_t *cpu_sort_space, int words_per_substring, int64_t lv2_num_items) {
    #pragma omp parallel for

    for (uint32_t i = 0; i < lv2_num_items; ++i) {
        permutation[i] = i;
    }

    for (int64_t iteration = words_per_substring - 1; iteration >= 0; --iteration) {
        uint32_t *lv2_substr_p = lv2_substrings + lv2_num_items * iteration;
        #pragma omp parallel for

        for (uint32_t i = 0; i < lv2_num_items; ++i) {
            cpu_sort_space[i] = uint64_t(*(lv2_substr_p + permutation[i])) << 32;
            cpu_sort_space[i] |= i;
        }

        // pss::parallel_stable_sort(cpu_sort_space, cpu_sort_space + lv2_num_items, CompareHigh32Bits());
        __gnu_parallel::sort(cpu_sort_space, cpu_sort_space + lv2_num_items);

        #pragma omp parallel for

        for (uint32_t i = 0; i < lv2_num_items; ++i) {
            cpu_sort_space[i] &= 0xFFFFFFFFULL;
            cpu_sort_space[i] |= uint64_t(permutation[cpu_sort_space[i]]) << 32;
        }

        #pragma omp parallel for

        for (uint32_t i = 0; i < lv2_num_items; ++i) {
            permutation[i] = cpu_sort_space[i] >> 32;
        }
    }
}

static const int kRadixBits = 8;
static const int kThreN = 64;
static const int kRadixMask = (1 << kRadixBits) - 1;
static const int kRadixBin = 1 << kRadixBits;

inline bool less_than_(uint32_t *substr, int64_t lv2_num_items, int which_word, uint32_t i, uint32_t j) {
    while (which_word-- >= 0) {
        if (substr[i] < substr[j]) return true;
        else if (substr[i] > substr[j]) return false;
        substr += lv2_num_items;
    }
    return false;
}

template<int kShift>
inline void lv2_cpu_radix_sort_st_core_(uint32_t *substr, uint32_t *permutation, int64_t n, int64_t lv2_num_items, int which_word) {
    if (n > kThreN) {
        int64_t count[kRadixBin] = {0};
        uint32_t *bin_s[kRadixBin], *bin_e[kRadixBin];
        for (int64_t i = 0; i < n; ++i) ++count[substr[permutation[i]] >> kShift & kRadixMask];
        bin_s[0] = permutation; bin_e[kRadixBin - 1] = permutation + n;
        for (int i = 1; i < kRadixBin; ++i) bin_e[i-1] = bin_s[i] = bin_s[i-1] + count[i-1];

        for (int i = 0; i < kRadixBin; ++i) {
            while (bin_s[i] < bin_e[i]) {
                uint32_t swapper = *bin_s[i];
                int bin_tag = substr[swapper] >> kShift & kRadixMask;
                if (bin_tag != i) {
                    do {
                        std::swap(swapper, *(bin_s[bin_tag]++));
                    } while ((bin_tag = substr[swapper] >> kShift & kRadixMask) != i);
                    *bin_s[i] = swapper;
                }
                ++bin_s[i];
            }
        }

        if (kShift > 0) {
            for (int i = 0; i < kRadixBin; ++i) {
                lv2_cpu_radix_sort_st_core_<(kShift > kRadixBits ? kShift - kRadixBits : 0)>(substr, 
                    bin_e[i] - count[i], count[i], lv2_num_items, which_word);
            }
        } else if (which_word > 0) {
            for (int i = 0; i < kRadixBin; ++i) {
                lv2_cpu_radix_sort_st_core_<24>(substr + lv2_num_items, 
                    bin_e[i] - count[i], count[i], lv2_num_items, which_word - 1);
            }
        }
    } else {
        for (int64_t i = 1; i < n; ++i) {
            if (less_than_(substr, lv2_num_items, which_word, permutation[i], permutation[i - 1])) {
                uint32_t tmp = permutation[i];
                permutation[i] = permutation[i-1];
                int64_t j;
                for (j = i - 1; j > 0 && less_than_(substr, lv2_num_items, which_word, tmp, permutation[j - 1]); --j) {
                    permutation[j] = permutation[j-1];
                }
                permutation[j] = tmp;
            }
        }
    }
}

inline void lv2_cpu_radix_sort_st(uint32_t *lv2_substrings, uint32_t *permutation, uint32_t *cpu_sort_space, uint64_t *buckets, int words_per_substring, int64_t lv2_num_items) {
    for (uint32_t i = 0; i < lv2_num_items; ++i) {
        permutation[i] = i;
    }
    lv2_cpu_radix_sort_st_core_<24>(lv2_substrings, permutation, lv2_num_items, lv2_num_items, words_per_substring - 1);
}
