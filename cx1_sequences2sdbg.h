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


#ifndef CX1_SEQUENCES2SDBG_H__
#define CX1_SEQUENCES2SDBG_H__

#include <pthread.h>
#include <stdint.h>
#include <vector>
#include <string>
#include "mac_pthread_barrier.h"
#include "definitions.h"
#include "cx1.h"
#include "sdbg_builder_writers.h"
#include "sequence_manager.h"

struct sequences2sdbg_opt_t {
    double host_mem;
    double gpu_mem;
    int num_cpu_threads;
    int num_output_threads;
    int kmer_k;
    int kmer_from;
    std::string contig_file_name;
    std::string add_contig_file_name;
    std::string input_prefix;
    std::string output_prefix;
    int mem_flag;

    sequences2sdbg_opt_t() {
        host_mem = 0;
        gpu_mem = 0;
        num_cpu_threads = 0;
        num_output_threads = 0;
        kmer_k = 0;
        kmer_from = 0;
        mem_flag = 1;
        contig_file_name = "";
        add_contig_file_name = "";
        input_prefix = "";
        output_prefix = "out";
    }
};

namespace cx1_sequences2sdbg {

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

struct sequences2sdbg_global_t {
    CX1<sequences2sdbg_global_t, kNumBuckets> cx1;

    // input options
    int kmer_k;
    int kmer_from;
    int num_cpu_threads;
    int num_output_threads;
    int64_t host_mem;
    int64_t gpu_mem;
    int mem_flag;

    std::string contig_file_name;
    std::string add_contig_file_name;
    std::string input_prefix;
    std::string output_prefix;

    int64_t num_seq;
    int64_t words_per_substring; // substrings to be sorted by GPU
    int64_t max_bucket_size;
    int64_t tot_bucket_size;
    int words_per_dummy_node;

    int bits_for_offset; // the number of bits needed to store the position in a contig

    // big arrays
    SequencePackage package;
    SequenceManager seq_manager;

    int32_t* lv1_items; // each item is an offset (read ID and position) in differential representation
    uint32_t* lv2_substrings; // stripped format
    uint32_t* lv2_substrings_db; // double buffer
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
    int64_t mem_packed_seq;

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

int64_t encode_lv1_diff_base(int64_t read_id, sequences2sdbg_global_t &g);
void    read_seq_and_prepare(sequences2sdbg_global_t &g); // num_items_, num_cpu_threads_ and num_output_threads_ must be set here
void*   lv0_calc_bucket_size(void*); // pthread working function
void    init_global_and_set_cx1(sequences2sdbg_global_t &g);
void*   lv1_fill_offset(void*); // pthread working function
void*   lv2_extract_substr(void*); // pthread working function
void    lv2_sort(sequences2sdbg_global_t &g);
void    lv2_pre_output_partition(sequences2sdbg_global_t &g);
void*   lv2_output(void*); // pthread working function
void    lv2_post_output(sequences2sdbg_global_t &g);
void    post_proc(sequences2sdbg_global_t &g);

} // end of namespace cx1_sequences2sdbg
#endif // CX1_SEQUENCES2SDBG_H__