#ifndef CX1_EDGE2SDBG_H__
#define CX1_EDGE2SDBG_H__

#include <pthread.h>
#include <stdint.h>
#include <vector>
#include <string>
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

struct edge2sdbg_global_t {
    CX1<edge2sdbg_global_t, kNumBuckets> cx1;

    // input options
    int max_read_length;
    int kmer_k;
    int kmer_freq_threshold;
    int num_cpu_threads;
    int num_output_threads;
    int64_t host_mem;
    int64_t gpu_mem;
    int mem_flag;
    int num_k1_per_read; // = max_read_length - kmer_k

    const char *input_file;
    const char *output_prefix;

    int words_per_edge; // number of (32-bit) words needed to represent a (k+1)-mer
    int64_t words_per_substring; // substrings to be sorted by GPU
    int offset_num_bits; // the number of bits needed to store the offset of a base in the read/(k+1)-mer (i.e. log(read_length))
    int64_t capacity;
    int64_t max_bucket_size;
    int64_t tot_bucket_size;
    int words_per_read; // number of (32-bit) words needed to represent a read in 2-bit-per-char format
    int read_length_mask;
    int k_num_bits; // the number of bits needed to store the position of the first $ in the kmer (i.e. log2(kmer_k+1))
    int64_t num_reads; // total number of reads

    // big arrays
    edge_word_t* packed_reads;
#ifndef LONG_READS
    unsigned char *first_0_out;
    unsigned char *last_0_in; // first potential 0-out-degree k-mer and last potential 0-in-degree k-mer
#else
    uint16_t *first_0_out;
    uint16_t *last_0_in; 
#endif
    int32_t* lv1_items; // each item is an offset (read ID and position) in differential representation

    int64_t *lv2_read_info; // to store where this lv2_item (k+1)-mer come from
    int64_t *lv2_read_info_db; // double buffer
    edge_word_t* lv2_substrings; // stripped format
    edge_word_t* lv2_substrings_db; // double buffer
    uint32_t* permutation; // permutation of { 1, ..., lv2_num_items }. for sorting (as value in a key-value pair)
    uint32_t* permutation_db;    // double buffer

#ifdef DISABLE_GPU
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
    int output_threads;

    // output
    WordWriter *word_writer;
};

int64_t encode_lv1_diff_base(int64_t read_id, edge2sdbg_global_t &g);
void    read_input_prepare(edge2sdbg_global_t &g); // num_items_, num_cpu_threads_ and num_output_threads_ must be set here
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