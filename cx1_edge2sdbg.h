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


#ifndef CX1_EDGE2SDBG_H__
#define CX1_EDGE2SDBG_H__

#include <pthread.h>
#include <stdint.h>
#include <vector>
#include <string>
#include "mac_pthread_barrier.h"
#include "definitions.h"
#include "cx1.h"
#include "sdbg_builder_writers.h"

struct edge2sdbg_opt_t {
    bool need_mercy;
    double host_mem;
    double gpu_mem;
    int num_edge_files;
    int num_cpu_threads;
    int num_output_threads;
    int max_read_length;
    std::string input_prefix;
    std::string output_prefix;
    int mem_flag;

    edge2sdbg_opt_t() {
        need_mercy = false;
        host_mem = 0;
        gpu_mem = 0;
        num_edge_files = 0;
        num_cpu_threads = 0;
        num_output_threads = 0;
        max_read_length = 120;
        mem_flag = 1;
        input_prefix = "";
        output_prefix = "out";
    }
};

namespace cx1_edge2sdbg {

static const int kBucketPrefixLength = 8; // less than 16 (chars per word)
static const int kBucketBase = 5;
static const int kNumBuckets = 390625; // pow(5, 8)
// binary search look up table
static const int kLookUpPrefixLength = 12;
static const int kLookUpShift = 32 - kLookUpPrefixLength * 2;
static const int kLookUpSize = 1 << (2 * kLookUpPrefixLength);

static const int64_t kMaxDummyEdges = 4294967294LL;

static const int64_t kMinLv2BatchSize = 2 * 1024 * 1024;
static const int64_t kMinLv2BatchSizeGPU = 64 * 1024 * 1024;
static const int64_t kDefaultLv1ScanTime = 8;
static const int64_t kMaxLv1ScanTime = 64;
static const int kSentinelValue = 4;
static const int kBWTCharNumBits = 3;
static const int kTopCharShift = kBitsPerEdgeWord - kBitsPerEdgeChar;

struct edge2sdbg_global_t {
    CX1<edge2sdbg_global_t, kNumBuckets> cx1;

    // input options
    int max_read_length;
    int kmer_k;
    int num_edge_files;
    int num_cpu_threads;
    int num_output_threads;
    int64_t host_mem;
    int64_t gpu_mem;
    int mem_flag;
    bool need_mercy;

    const char *input_prefix;
    const char *output_prefix;

    int64_t num_edges;
    int words_per_edge; // number of (32-bit) words needed to represent a (k+1)-mer
    int64_t words_per_substring; // substrings to be sorted by GPU
    int64_t capacity;
    int64_t max_bucket_size;
    int64_t tot_bucket_size;
    int words_per_dummy_node;
    int mult_mem_type; // 0: compact with (k+1)-mer; 1: use extra 8 bits; 2: use extra 16 bits

    int k_num_bits; // the number of bits needed to store the position of the first $ in the kmer (i.e. log2(kmer_k+1))

    // big arrays
    edge_word_t* packed_edges;
    uint8_t *multiplicity8;    // store multiplicity if 8 additional bits are needed
    uint16_t *multiplicity16;   // store multiplicity if 16 additional bits are needed

    int32_t* lv1_items; // each item is an offset (read ID and position) in differential representation
    edge_word_t* lv2_substrings; // stripped format
    edge_word_t* lv2_substrings_db; // double buffer
    uint32_t* permutation; // permutation of { 1, ..., lv2_num_items }. for sorting (as value in a key-value pair)
    uint32_t* permutation_db;    // double buffer

#ifndef USE_GPU
    uint64_t *cpu_sort_space;
#else
    void *gpu_key_buffer1;
    void *gpu_key_buffer2;
    void *gpu_value_buffer1;
    void *gpu_value_buffer2;
#endif

    pthread_mutex_t lv1_items_scanning_lock;
    int64_t lv2_num_items_db;
    // memory usage
    int64_t mem_packed_edges;

    // statistics
    int64_t num_chars_in_w[9];
    int64_t num_ones_in_last;
    int64_t total_number_edges;
    int64_t num_dollar_nodes;
    int64_t num_dummy_edges;
    int cur_prefix;
    int cur_suffix_first_char;

    // output
    DBG_BinaryWriter sdbg_writer;
    WordWriter dummy_nodes_writer;
    FILE *output_f_file;
    FILE *output_multiplicity_file;
    FILE *output_multiplicity_file2;

    unsigned char *lv2_aux;
    pthread_barrier_t output_barrier;
};

int64_t encode_lv1_diff_base(int64_t read_id, edge2sdbg_global_t &g);
void    read_edges_and_prepare(edge2sdbg_global_t &g); // num_items_, num_cpu_threads_ and num_output_threads_ must be set here
void*   lv0_calc_bucket_size(void*); // pthread working function
void    init_global_and_set_cx1(edge2sdbg_global_t &g);
void*   lv1_fill_offset(void*); // pthread working function
void*   lv2_extract_substr(void*); // pthread working function
void    lv2_sort(edge2sdbg_global_t &g);
void    lv2_pre_output_partition(edge2sdbg_global_t &g);
void*   lv2_output(void*); // pthread working function
void    lv2_post_output(edge2sdbg_global_t &g);
void    post_proc(edge2sdbg_global_t &g);

} // end of namespace cx1_edge2sdbg
#endif // CX1_EDGE2SDBG_H__