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

#ifndef CX1_KMER_COUNT_H__
#define CX1_KMER_COUNT_H__

#include <pthread.h>
#include <stdint.h>
#include <vector>
#include <string>
#include "definitions.h"
#include "cx1.h"
#include "edge_io.h"
#include "sequence_package.h"
#include "lib_info.h"

struct count_opt_t {
    int kmer_k;
    int kmer_freq_threshold;
    double host_mem;
    double gpu_mem;
    int num_cpu_threads;
    int num_output_threads;
    std::string read_lib_file;
    std::string assist_seq_file;
    std::string output_prefix;
    int mem_flag;
    bool need_mercy;

    count_opt_t() {
        kmer_k = 21;
        kmer_freq_threshold = 2;
        host_mem = 0;
        gpu_mem = 0;
        num_cpu_threads = 0;
        num_output_threads = 0;
        read_lib_file = "";
        output_prefix = "out";
        mem_flag = 1;
        need_mercy = true;
    }
};

namespace cx1_kmer_count {

static const int kBucketPrefixLength = 8;
static const int kBucketBase = 4;
static const int kNumBuckets = 65536; // pow(4, 8)
static const int64_t kMinLv2BatchSize = 64 * 1024 * 1024;
static const int64_t kMinLv2BatchSizeGPU = 64 * 1024 * 1024;
static const int64_t kDefaultLv1ScanTime = 8;
static const int64_t kMaxLv1ScanTime = 64;
static const int kSentinelValue = 4;

#define LONG_READS
#ifdef LONG_READS
    static const uint32_t kSentinelOffset = 4294967295U;
#else
    static const uint32_t kSentinelOffset = 255;
#endif

struct count_global_t {
    CX1<count_global_t, kNumBuckets> cx1;

    // input options
    int max_read_length;
    int kmer_k;
    int kmer_freq_threshold;
    int num_cpu_threads;
    int num_output_threads;
    int64_t host_mem;
    int64_t gpu_mem;
    int mem_flag;
    std::string read_lib_file;
    std::string assist_seq_file;
    std::string output_prefix;

    int words_per_edge; // number of (32-bit) words needed to represent a (k+1)-mer
    int64_t words_per_substring; // substrings to be sorted by GPU
    int offset_num_bits; // the number of bits needed to store the offset of a base in the read/(k+1)-mer (i.e. log(read_length))
    int64_t max_bucket_size;
    int64_t tot_bucket_size;
    int read_length_mask;
    int64_t num_reads; // total number of reads

    // big arrays
    SequencePackage package;
    std::vector<lib_info_t> lib_info;

    // lv1 new sorting scheme
    int64_t max_sorting_items;
    int64_t mem_sorting_items;

    uint32_t *substr_all;
    uint32_t *permutations_all;
    uint32_t *cpu_sort_space_all;
    int64_t *readinfo_all;

#ifndef LONG_READS
    uint8_t *first_0_out;
    uint8_t *last_0_in; // first potential 0-out-degree k-mer and last potential 0-in-degree k-mer
#else
    uint32_t *first_0_out;
    uint32_t *last_0_in;
#endif
    int32_t *lv1_items; // each item is an offset (read ID and position) in differential representation

    int64_t *lv2_read_info; // to store where this lv2_item (k+1)-mer come from
    int64_t *lv2_read_info_db; // double buffer
    uint32_t *lv2_substrings; // stripped format
    uint32_t *lv2_substrings_db; // double buffer
    uint32_t *permutation; // permutation of { 1, ..., lv2_num_items }. for sorting (as value in a key-value pair)
    uint32_t *permutation_db;    // double buffer

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
    int64_t mem_packed_reads;

    // stat
    int64_t *edge_counting; // count the number of (k+1)mer with occurs i times
    int64_t *thread_edge_counting;

    // output
    EdgeWriter edge_writer;
};

int64_t encode_lv1_diff_base(int64_t read_id, count_global_t &g);
void    read_input_prepare(count_global_t &g); // num_items_, num_cpu_threads_ and num_output_threads_ must be set here
void   *lv0_calc_bucket_size(void *); // pthread working function
void    init_global_and_set_cx1(count_global_t &g);
void   *lv1_fill_offset(void *); // pthread working function
void    lv1_direct_sort_and_count(count_global_t &g);
void   *lv2_extract_substr(void *); // pthread working function
void    lv2_sort(count_global_t &g);
void    lv2_pre_output_partition(count_global_t &g);
void   *lv2_output(void *); // pthread working function
void    lv2_post_output(count_global_t &g);
void    post_proc(count_global_t &g);

} // end of namespace cx1_kmer_count
#endif // CX1_KMER_COUNT_H__