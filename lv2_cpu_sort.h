/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong
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

inline void lv2_cpu_sort_st(uint32_t *lv2_substrings, uint32_t *permutation, uint64_t *cpu_sort_space, int words_per_substring, int64_t lv2_num_items) {
    for (uint32_t i = 0; i < lv2_num_items; ++i) {
        permutation[i] = i;
    }

    for (int64_t iteration = words_per_substring - 1; iteration >= 0; --iteration) {
        uint32_t *lv2_substr_p = lv2_substrings + lv2_num_items * iteration;

        for (uint32_t i = 0; i < lv2_num_items; ++i) {
            cpu_sort_space[i] = uint64_t(*(lv2_substr_p + permutation[i])) << 32;
            cpu_sort_space[i] |= i;
        }
        // pss::parallel_stable_sort(cpu_sort_space, cpu_sort_space + lv2_num_items, CompareHigh32Bits());
        std::sort(cpu_sort_space, cpu_sort_space + lv2_num_items);

        for (uint32_t i = 0; i < lv2_num_items; ++i) {
            cpu_sort_space[i] &= 0xFFFFFFFFULL;
            cpu_sort_space[i] |= uint64_t(permutation[cpu_sort_space[i]]) << 32;
        }

        for (uint32_t i = 0; i < lv2_num_items; ++i) {
            permutation[i] = cpu_sort_space[i] >> 32;
        }
    }
}

struct CmpSubStr {
    uint32_t *substr;
    int64_t num_items;
    int words_per_substring;
    CmpSubStr(uint32_t *substr, int64_t num_items, int words_per_substring):
        substr(substr), num_items(num_items), words_per_substring(words_per_substring) {}

    bool operator() (uint32_t x, uint32_t y) {
        int64_t idx = 0;
        for (int i = 0; i < words_per_substring; ++i, idx += num_items) {
            uint32_t xx = substr[x + idx];
            uint32_t yy = substr[y + idx];
            if (xx > yy) { return false; }
            else if (yy > xx) { return true; }
        }
        return false;
    }
};

inline void sort_digit(uint32_t *arr, uint32_t *permutation, uint32_t *buf, uint64_t *buckets, int64_t num_items, int shift_bits) {
    memset(buckets, 0, sizeof(buckets[0]) * (1 << 16));
    for (int64_t i = 0; i < num_items; ++i) {
        buckets[(arr[permutation[i]] >> shift_bits) & 0xFFFF]++;
    }

    int64_t acc = 0;
    for (unsigned i = 0; i < (1 << 16); ++i) {
        int64_t tmp = acc;
        acc += buckets[i];
        buckets[i] = tmp;
    }

    for (int64_t i = 0; i < num_items; ++i) {
        buf[buckets[(arr[permutation[i]] >> shift_bits) & 0xFFFF]++] = permutation[i];
    }
}

inline void lv2_cpu_radix_sort_st(uint32_t *lv2_substrings, uint32_t *permutation, uint32_t *cpu_sort_space, uint64_t *buckets, int words_per_substring, int64_t lv2_num_items) {
    for (uint32_t i = 0; i < lv2_num_items; ++i) {
        permutation[i] = i;
    }

    if (lv2_num_items < 65536 * 2) {
        std::sort(permutation, permutation + lv2_num_items, CmpSubStr(lv2_substrings, lv2_num_items, words_per_substring));
        return;
    }

    for (int64_t iteration = words_per_substring - 1; iteration >= 0; --iteration) {
        uint32_t *lv2_substr_p = lv2_substrings + lv2_num_items * iteration;

        // 1 pass low  16 bits
        sort_digit(lv2_substr_p, permutation, cpu_sort_space, buckets, lv2_num_items, 0);
        // 2 pass high 16 bits
        sort_digit(lv2_substr_p, cpu_sort_space, permutation, buckets, lv2_num_items, 16);
    }
}