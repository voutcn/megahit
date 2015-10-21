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

#ifndef SDBG_BUILDER_LV2_GPU_FUNCTIONS_H__
#define SDBG_BUILDER_LV2_GPU_FUNCTIONS_H__

#include "definitions.h"
// single thread
void cuda_init();
void get_cuda_memory(size_t &free_mem, size_t &total_mem);
void alloc_gpu_buffers(void *&gpu_key_buffer1, void *&gpu_key_buffer2, void *&gpu_value_buffer1, void *&gpu_value_buffer2, size_t max_num_items);
void free_gpu_buffers(void *gpu_key_buffer1, void *gpu_key_buffer2, void *gpu_value_buffer1, void *gpu_value_buffer2);
void lv2_gpu_sort(uint32_t *lv2_substrings, uint32_t *permutation, int words_per_substring, int64_t lv2_num_items,
                  void *gpu_key_buffer1, void *gpu_key_buffer2, void *gpu_value_buffer1, void *gpu_value_buffer2);

#endif // SDBG_BUILDER_LV2_GPU_FUNCTIONS_H__