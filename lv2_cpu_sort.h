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
// #include "parallel_stable_sort/parallel_stable_sort.h"

struct CompareHigh32Bits {
    bool operator() (uint64_t a, uint64_t b) {
        return (a >> 32) < (b >> 32);
    }
};

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