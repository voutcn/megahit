#ifndef SDBG_BUILDER_COUNT_H__
#define SDBG_BUILDER_COUNT_H__

#include <vector>
#include <string>
#include "cx1.h"

class std::string;
class std::vector;

struct count_opt_t {
    int kmer_k;
    int min_edge_freq;
    double host_mem;
    double gpu_mem;
    int max_read_length;
    int num_cpu_threads;
    int num_output_threads;
    std::string input_file;
    std::string output_prefix;
    int mem_flag;
    bool need_mercy;

    count_opt_t() {
        kmer_k = 21;
        min_edge_freq = 2;
        host_mem = 0;
        gpu_mem = 0;
        max_read_length = 120;
        num_cpu_threads = 0;
        num_output_threads = 0;
        input_file = "";
        output_prefix = "out";
        mem_flag = 1;
        need_mercy = true;
    }
};

struct count_global_t {
    // input options
    int max_read_length;
    int kmer_k;
    int kmer_freq_threshold;
    int num_cpu_threads;
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
    int64_t *lv2_read_info; // to store where this lv2_item (k+1)-mer come from
    int64_t *lv2_read_info_db;
    int64_t* bucket_sizes; // the number of items of each bucket
    int* lv1_items; // each item is an offset (read ID and position) in differential representation
    std::vector<int64_t> lv1_items_special; // if the differential > 2^31, store full offset in this vector
    int64_t* lv2_items; // each item is an offset (read ID and position), full expressed
    edge_word_t* lv2_substrings; // stripped format
    edge_word_t* lv2_substrings_db; // double buffer
    uint32_t* permutation; // permutation of { 1, ..., lv2_num_items }. for sorting (as value in a key-value pair)
    uint32_t* permutation_db;    // dump for double buffer

#ifdef DISABLE_GPU
    uint64_t *cpu_sort_space;
#else
    void *gpu_key_buffer1;
    void *gpu_key_buffer2;
    void *gpu_value_buffer1;
    void *gpu_value_buffer2;
#endif

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
    int64_t lv2_num_items_db;

    // memory usage
    int64_t mem_packed_reads;

    // stat
    int64_t *edge_counting; // count the number of (k+1)mer with occurs i times
    int64_t *thread_edge_counting;
    int output_threads;

    // output
    WordWriter word_writer[kMaxNumCPUThreads];
};

#endif // SDBG_BUILDER_COUNT_H__