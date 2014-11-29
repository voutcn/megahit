/*
 *  MEGAHIT
 *  Copyright (C) 2014 The University of Hong Kong
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

#ifndef SDBG_BUILDER_UTIL_H_
#define SDBG_BUILDER_UTIL_H_

#include <vector>

#include "MAC_pthread_barrier.h"
#include "definitions.h"
#include "timer.h"
#include "sdbg_builder_writers.h"

struct global_data_t;
struct readpartition_data_t;
struct bucketpartition_data_t;
struct outputpartition_data_t;

// local data for each read partition (i.e. a subrange of input reads)
struct readpartition_data_t {
    struct global_data_t* globals;
    int rp_id; // ID of this read partition, in [ 0, num_cpu_threads ).
    pthread_t thread;
    int64_t rp_start_id, rp_end_id; // start and end IDs of this read partition (end is exclusive)
    int64_t* rp_bucket_sizes; // bucket sizes for this read partition
    int64_t* rp_bucket_offsets;
    int64_t rp_lv1_differential_base; // the initial offset globals.lv1_items
};

// local data for each bucket partition (i.e. a range of buckets), used in lv.2 (extract substring)
struct bucketpartition_data_t {
    struct global_data_t* globals;
    int bp_id;
    pthread_t thread;
    int bp_start_bucket, bp_end_bucket;
};

struct outputpartition_data_t {
    struct global_data_t *globals;
    int op_id;
    int64_t op_start_index, op_end_index;
};

// global data and structures
struct global_data_t {
    // input options
    int run_mode;
    int max_read_length;
    int kmer_k;
    int kmer_freq_threshold;
    int num_cpu_threads;
    int64_t host_mem;
    int64_t gpu_mem;

    const char *input_file;
    const char *phase2_input_prefix;
    const char *output_prefix;

    //----------common variables for two phases--------------------
    int words_per_edge; // number of (32-bit) words needed to represent a (k+1)-mer
    int64_t words_per_substring; // substrings to be sorted by GPU
    int offset_num_bits; // the number of bits needed to store the offset of a base in the read/(k+1)-mer (i.e. log(read_length))
    size_t capacity;

    struct readpartition_data_t readpartitions[kMaxNumCPUThreads];
    struct bucketpartition_data_t bucketpartitions[kMaxNumCPUThreads];
    struct outputpartition_data_t outputpartitions[kMaxNumCPUThreads];
    pthread_t output_threads[kMaxNumCPUThreads];

    // big arrays
    int64_t* bucket_sizes; // the number of items of each bucket
    int* lv1_items; // each item is an offset (read ID and position) in differential representation
    std::vector<int64_t> lv1_items_special; // if the differential > 2^31, store full offset in this vector
    int64_t* lv2_items; // each item is an offset (read ID and position), full expressed
    edge_word_t* lv2_substrings; // stripped format
    uint32_t* permutation; // permutation of { 1, ..., lv2_num_items }. for sorting (as value in a key-value pair)
    edge_word_t* lv2_substrings_to_output; // dump for double buffer
    uint32_t* permutation_to_output;    // dump for double buffer
    uint64_t *cpu_sort_space;

    // memory resources used. computational limits.
    int64_t max_lv1_items;
    int64_t max_lv2_items;

    // Lv.1 variables
    int lv1_start_bucket, lv1_end_bucket; // end is exclusive
    int64_t lv1_num_items;
    pthread_mutex_t lv1_items_scanning_lock;  

    // Lv.2 variables
    int lv2_start_bucket, lv2_end_bucket; // end is exclusive
    int64_t lv2_num_items;
    int64_t lv2_num_items_to_output;

    // output
    int64_t lv2_output_start_index[kMaxNumCPUThreads];
    int64_t lv2_output_end_index[kMaxNumCPUThreads];

    //-------------end of common parameters for two phases--------------------

    //----------------phase1-----------------------
    int words_per_read; // number of (32-bit) words needed to represent a read in 2-bit-per-char format
    int read_length_mask;
    int k_num_bits; // the number of bits needed to store the position of the first $ in the kmer (i.e. log(kmer_k+1))
    int64_t num_reads; // total number of reads

    // big arrays
    edge_word_t* packed_reads;
    unsigned char *first_0_out;
    unsigned char *last_0_in; // first potential 0-out-degree k-mer and last potential 0-in-degree k-mer
    int64_t *lv2_read_info; // to store where this lv2_item (k+1)-mer come from
    int64_t *lv2_read_info_to_output;

    // memory usage
    int64_t mem_packed_reads;

    // output
    WordWriter word_writer[kMaxNumCPUThreads];

    // stat
    int64_t *edge_counting; // count the number of (k+1)mer with occurs i times
    int64_t *thread_edge_counting;
    int64_t num_dummy_edges;
    int64_t num_incoming_zero_nodes;
    int64_t num_outgoing_zero_nodes;
    int phase1_num_output_threads;

    //-----------end-of-phase1-----------------------

    //--------------phase2---------------------------
    int words_per_dummy_node;
    int mult_mem_type; // 0: compact with (k+1)-mer; 1: use extra 8 bits; 2: use extra 16 bits
    int64_t num_edges;

    // large array
    edge_word_t *packed_edges;
    uint8_t *multiplicity8;    // store multiplicity if 8 additional bits are needed
    uint16_t *multiplicity16;   // store multiplicity if 16 additional bits are needed

    // memory usage
    int64_t mem_packed_edges;

    // output
    DBG_BinaryWriter sdbg_writer;
    WordWriter dummy_nodes_writer;
    FILE *output_f_file;
    FILE *output_multiplicity_file;

    // statistics
    int64_t num_chars_in_w[9];
    int64_t num_ones_in_last;
    int64_t total_number_edges;
    int64_t num_dollar_nodes;
    int cur_prefix;
    int cur_suffix_first_char;

    // for output thread
    unsigned char *lv2_aux;
    int phase2_num_output_threads;
    xtimer_t phase2_output_timer;
    pthread_barrier_t output_barrier;

    // for lookup binary search on sorted edges
    int64_t *edge_lookup;
    bool need_mercy;
    //---------------end-of-phase2-------------------
};

// do not modify these constants
static const int kSentinelValue = 4; // 4 = $
static const int kBWTCharNumBits = 3; // need 3 bits to represent ACGT$
static const int kTopCharShift = kBitsPerEdgeWord - kBitsPerEdgeChar; // bits >> to get the most significant char
static const uint32_t kDifferentialLimit = 2147483647; // 32-bit signed max int
static const int kSignBitMask = 0x80000000; // the MSB of 32-bit

namespace phase1 {
// definitions
static const int kBucketPrefixLength = 8;
static const int kBucketBase = 4;
static const int kNumBuckets = 65536;

void Phase1Entry(struct global_data_t &globals);

}

namespace phase2 {
// definitions
static const int kBucketPrefixLength = 8; // less than 16 (chars per word)
static const int kBucketBase = 5;
static const int kNumBuckets = 390625;
// binary search look up table
static const int kLookUpPrefixLength = 12;
static const int kLookUpShift = 32 - kLookUpPrefixLength * 2;
static const int kLookUpSize = 1 << (2 * kLookUpPrefixLength);

static const int64_t kMaxDummyEdges = 4294967294LL;

void Phase2Entry(struct global_data_t &globals);

}

#endif // SDBG_BUILDER_UTIL_H_