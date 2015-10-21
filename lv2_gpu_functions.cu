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

#include "lv2_gpu_functions.h"
#include <stdio.h>
#include <assert.h>
#include "cub/util_allocator.cuh"
#include "cub/device/device_radix_sort.cuh"
#include "utils.h"

using namespace cub;

static CachingDeviceAllocator g_allocator(true);
const int kGPUThreadPerBlock = 256;

void cuda_init() {
    cudaFree(0);
}

void get_cuda_memory(size_t &free_mem, size_t &total_mem) {
    assert(cudaMemGetInfo(&free_mem, &total_mem) == cudaSuccess);
}

void alloc_gpu_buffers(void* &gpu_key_buffer1,
                       void* &gpu_key_buffer2,
                       void* &gpu_value_buffer1,
                       void* &gpu_value_buffer2,
                       size_t max_num_items) {
    CubDebugExit(g_allocator.DeviceAllocate((void**)&gpu_key_buffer1, sizeof(uint32_t) * max_num_items));
    CubDebugExit(g_allocator.DeviceAllocate((void**)&gpu_key_buffer2, sizeof(uint32_t) * max_num_items));
    CubDebugExit(g_allocator.DeviceAllocate((void**)&gpu_value_buffer1, sizeof(uint32_t) * max_num_items));
    CubDebugExit(g_allocator.DeviceAllocate((void**)&gpu_value_buffer2, sizeof(uint32_t) * max_num_items));
}

void free_gpu_buffers(void* gpu_key_buffer1,
                      void* gpu_key_buffer2,
                      void* gpu_value_buffer1,
                      void* gpu_value_buffer2) {
    // free device memory
    CubDebugExit(g_allocator.DeviceFree(gpu_key_buffer1));
    CubDebugExit(g_allocator.DeviceFree(gpu_key_buffer2));
    CubDebugExit(g_allocator.DeviceFree(gpu_value_buffer1));
    CubDebugExit(g_allocator.DeviceFree(gpu_value_buffer2));
}


// device function for permuting an array
__global__ void permutation_kernel(uint32_t *index, uint32_t *val, uint32_t *new_val, uint32_t num_items) {
    int tid = blockIdx.x * kGPUThreadPerBlock + threadIdx.x;
    if (tid < num_items)
        new_val[tid] = val[index[tid]];
}

// device function for reset permutation
__global__ void reset_permutation_kernel(uint32_t *permutation, uint32_t num_items) {
    int tid = blockIdx.x * kGPUThreadPerBlock + threadIdx.x;
    if (tid < num_items)
        permutation[tid] = tid;
}

// single thread
void lv2_gpu_sort(uint32_t *lv2_substrings,
                  uint32_t *permutation,
                  int words_per_substring,
                  int64_t lv2_num_items,
                  void* gpu_key_buffer1,
                  void* gpu_key_buffer2,
                  void* gpu_value_buffer1,
                  void* gpu_value_buffer2) {
    DoubleBuffer<uint32_t> d_keys;
    DoubleBuffer<uint32_t> d_values;
    d_keys.d_buffers[0] = static_cast<__typeof(d_keys.d_buffers[0])>(gpu_key_buffer1);
    d_keys.d_buffers[1] = static_cast<__typeof(d_keys.d_buffers[1])>(gpu_key_buffer2);
    d_values.d_buffers[0] = static_cast<__typeof(d_values.d_buffers[0])>(gpu_value_buffer1);
    d_values.d_buffers[1] = static_cast<__typeof(d_values.d_buffers[1])>(gpu_value_buffer2);

    // Initialize permutation array
    int num_gpu_blocks = DivCeiling(lv2_num_items, kGPUThreadPerBlock);
    reset_permutation_kernel<<<num_gpu_blocks, kGPUThreadPerBlock>>>(d_values.d_buffers[d_values.selector], lv2_num_items);

    // Allocate temporary storage
    size_t  temp_storage_bytes  = 0;
    void *gpu_temp_storage     = NULL;
    CubDebugExit(DeviceRadixSort::SortPairs(gpu_temp_storage, temp_storage_bytes, d_keys, d_values, lv2_num_items));
    CubDebugExit(g_allocator.DeviceAllocate(&gpu_temp_storage, temp_storage_bytes));

    for (int64_t iteration = words_per_substring - 1; iteration >= 0; --iteration) {
        if (iteration == words_per_substring - 1) { // first iteration
            CubDebugExit(cudaMemcpy(d_keys.d_buffers[d_keys.selector], lv2_substrings + (iteration * lv2_num_items),
                                    sizeof(uint32_t) * lv2_num_items, cudaMemcpyHostToDevice));
        } else {
            CubDebugExit(cudaMemcpy(d_keys.d_buffers[1 - d_keys.selector], lv2_substrings + (iteration * lv2_num_items),
                                    sizeof(uint32_t) * lv2_num_items, cudaMemcpyHostToDevice));

            permutation_kernel<<<num_gpu_blocks, kGPUThreadPerBlock>>>(d_values.d_buffers[d_values.selector],
                    d_keys.d_buffers[1 - d_keys.selector], d_keys.d_buffers[d_keys.selector],
                    lv2_num_items);
        }

        // Run
        CubDebugExit(DeviceRadixSort::SortPairs(gpu_temp_storage, temp_storage_bytes, d_keys, d_values, lv2_num_items));
    }

    // copy answer back to host
    CubDebugExit(cudaMemcpy(permutation, d_values.d_buffers[d_values.selector], sizeof(uint32_t) * lv2_num_items, cudaMemcpyDeviceToHost));

    CubDebugExit(g_allocator.DeviceFree(gpu_temp_storage));
}