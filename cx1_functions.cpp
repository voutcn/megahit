/*
 *  cx1_functions.cpp
 *  This file is a part of MEGAHIT
 *  
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

#include <pthread.h>
#include <zlib.h>
#include <assert.h>
#include <algorithm>
#include <string>
#include <omp.h>

#include "timer.h"
#include "definitions.h"
#include "io-utility.h"
#include "helper_functions-inl.h"
#include "mem_file_checker-inl.h"
#include "sdbg_builder_util.h"
#include "sdbg_builder_writers.h"
#include "kmer_uint32.h"
#include "lv2_cpu_sort.h"
#include "MAC_pthread_barrier.h"
#include "kseq.h"

#ifndef DISABLE_GPU
#include "lv2_gpu_functions.h"
#endif

// *** func & def commonly used in both phases ***
#define PACKED_READS(i, globals) ((globals).packed_reads + (i) * (globals).words_per_read)
#define GPU_BYTES_PER_ITEM 16 // key & value, 4 bytes each. double for radix sort internal buffer
#define LV1_BYTES_PER_ITEM 4 // 32-bit differential offset

int sdbg_builder_verbose = 3;

/**
 * @brief extract the nth char in a packed read/edge
 */
inline int ExtractNthChar(edge_word_t *read_ptr, int n) {
    int which_word = n / kCharsPerEdgeWord;
    int index_in_word = n % kCharsPerEdgeWord;
    return (read_ptr[which_word] >> (kBitsPerEdgeChar * (kCharsPerEdgeWord - 1 - index_in_word))) & kEdgeCharMask;
}

/**
 * @brief map ACGT to 0123
 */
char dna_map[256];
void InitDNAMap() {
    dna_map['A'] = dna_map['a'] = 0;
    dna_map['C'] = dna_map['c'] = 1;
    dna_map['G'] = dna_map['g'] = 2;
    dna_map['T'] = dna_map['t'] = 3;
    dna_map['N'] = dna_map['n'] = 2; // -> G
}

/**
 * @brief auto adjust lv1 and lv2 number items according to available memory
 */
void GetNumItems(struct global_data_t &globals, int64_t mem_remained, int64_t lv2_bytes_per_item) {
    int64_t min_lv1_items = globals.tot_bucket_size / (kMaxLv1ScanTime - 0.5);
    int64_t min_lv2_items = std::max(globals.max_bucket_size, kMinLv2BatchSize);

    // --- adjust max_lv2_items to fit memory ---
    while (globals.max_lv2_items >= min_lv2_items) {
        int64_t mem_lv2 = lv2_bytes_per_item * globals.max_lv2_items;
        if (mem_remained <= mem_lv2) {
            globals.max_lv2_items *= 0.9;
            continue;
        }

        globals.max_lv1_items = (mem_remained - mem_lv2) / LV1_BYTES_PER_ITEM;
        if (globals.max_lv1_items < min_lv1_items || 
            globals.max_lv1_items < globals.max_lv2_items) {
            globals.max_lv2_items *= 0.9;
        } else {
            break;
        }
    }

    if (globals.max_lv2_items < min_lv2_items) {
        err("[ERROR B::%s] No enough memory to process.\n", __func__);
        exit(1);
    }

    // --- adjust max_lv2_items to fit more lv1 item ---
    // TODO: 4 is arbitrary chosen, not fine tune
    while (globals.max_lv2_items * 4 > globals.max_lv1_items) {
        if (globals.max_lv2_items * 0.95 >= min_lv2_items) {
            globals.max_lv2_items *= 0.95;
            globals.max_lv1_items = (mem_remained - lv2_bytes_per_item * globals.max_lv2_items) / LV1_BYTES_PER_ITEM;
        } else {
            break;
        }
    }
}

/**
 * @param end_limit the max allowed value for end_bucket
 * @param num_items total num of items written to num_items by reference
 * @return end_bucket such that bucket_sizes[start_bucket..end_bucket-1] sums up to at most item_limit
 */
int FindEndBucket(int64_t * bucket_sizes, int start_bucket, int end_limit, int64_t item_limit, int64_t &num_items) {    
    num_items = 0;
    int end_bucket = start_bucket;
    while (end_bucket < end_limit) { // simple linear scan
        if (num_items + bucket_sizes[end_bucket] > item_limit) {
            return end_bucket;
        }
        num_items += bucket_sizes[end_bucket];
        end_bucket++;
    }
    return end_limit;
}

// single thread helper function
void Lv1ComputeBucketOffset(struct global_data_t &globals) {
    // compute "global" (thread 0) offsets first
    int64_t *offsets = globals.readpartitions[0].rp_bucket_offsets;
    offsets[globals.lv1_start_bucket] = 0;
    for (int b = globals.lv1_start_bucket + 1; b < globals.lv1_end_bucket; ++b) {
        offsets[b] = offsets[b-1] + globals.bucket_sizes[b-1]; // accumulate
    }
    // then for each read partition
    for (int t = 1; t < globals.num_cpu_threads; ++t) {
        int64_t *this_offsets = globals.readpartitions[t].rp_bucket_offsets;
        int64_t *prev_offsets = globals.readpartitions[t-1].rp_bucket_offsets;
        int64_t *sizes = globals.readpartitions[t-1].rp_bucket_sizes;
        for (int b = globals.lv1_start_bucket; b < globals.lv1_end_bucket; ++b) {
            this_offsets[b] = prev_offsets[b] + sizes[b];
        }
    }
}

// single thread helper function
void Lv2DistributeBucketPartitions(struct global_data_t &globals, int num_output_threads) {
    int64_t average = globals.lv2_num_items / (globals.num_cpu_threads - num_output_threads); // recall: we only have globals.num_cpu_threads-globals.phase1_num_output_threads bucketpartitions
    int bucket = globals.lv2_start_bucket;
    for (int t = 0; t < globals.num_cpu_threads-num_output_threads-1; ++t) {
        int64_t num_items = 0;
        globals.bucketpartitions[t].bp_start_bucket = bucket;
        while (bucket < globals.lv2_end_bucket) {
            num_items += globals.bucket_sizes[bucket++];
            if (num_items >= average) {
                break;
            }
        }
        globals.bucketpartitions[t].bp_end_bucket = bucket;
    }
    // last
    globals.bucketpartitions[globals.num_cpu_threads-num_output_threads-1].bp_start_bucket = bucket;
    globals.bucketpartitions[globals.num_cpu_threads-num_output_threads-1].bp_end_bucket = globals.lv2_end_bucket;
}

namespace phase1 {
// debug
void DumpPackedReads(edge_word_t* p, int read_length) {
    for (int i=0; i<read_length; ++i) {
        int j = i;
        log("%c", dna_chars[ (*(p+(j/kCharsPerEdgeWord)) >> ((kCharsPerEdgeWord-1-j%kCharsPerEdgeWord)*kBitsPerEdgeChar)) & kEdgeCharMask ]);
    }
}

/*
 * Packs an ASCII read into 2-bit per base form. The last base goes to the MSB of the first word.
 * -Params-
 * read: the read, in ASCII ACGT
 * p: a pointer to the first edge_word_t of the packed sequence to be written to
 * read_length: number of bases in the read
 * last_word_shift: the number of empty bits in the last word. we need to shift the bits up by this amount. solely determined by read_length.
 */
inline void PackReadFromAscii(char* read, edge_word_t* p, int read_length, int words_per_read) {
    // for de Bruijn graph construction, packing the reverse is more convenient
    edge_word_t w = 0;
    int i, j;
    for (i = 0, j = 0; j < read_length; ++j) {
        if (j % kCharsPerEdgeWord == 0 && j) { // TODO bitwise?
            *(p++) = w;
            w = 0;
            ++i;
        }
        while (read[j] == 'N') {
            break;
        }
        w = (w << kBitsPerEdgeChar) | dna_map[ (int)read[ j ] ];
    }

    int last_word_shift = j % kCharsPerEdgeWord;
    last_word_shift = last_word_shift ? (kCharsPerEdgeWord - last_word_shift) * kBitsPerEdgeChar : 0;
    *p = w << last_word_shift;

    while (++i < words_per_read) {
        *(++p) = 0;
    }

    *p |= read_length;
}

/**
 * @brief encode read_id and its offset in one int64_t
 */
inline int64_t EncodeOffset(int64_t read_id, int offset, int strand, int length_num_bits) {
    return (read_id << (length_num_bits + 1)) | (offset << 1) | strand;
}

inline int GetReadLength(edge_word_t* read_p, int words_per_read, int mask) {
    return *(read_p + words_per_read - 1) & mask;
}

inline int GetReadLengthByID(int64_t id, global_data_t &globals) {
    return *(globals.packed_reads + (id + 1) * globals.words_per_read - 1) & globals.read_length_mask;
}

/**
 * @brief read fastx queries and pack
 */
void ReadInputFile(struct global_data_t &globals) {
    int64_t num_reads = 0;
    int bits_read_length = 1; // bit needed to store read_length
    while ((1 << bits_read_length) - 1 < globals.max_read_length) {
        ++bits_read_length;
    }
    globals.words_per_read = DivCeiling(globals.max_read_length * kBitsPerEdgeChar + bits_read_length, kBitsPerEdgeWord);
    int64_t max_num_reads = globals.host_mem / (sizeof(edge_word_t) * globals.words_per_read) * 3 / 4; //TODO: more accurate
    int read_length;
    edge_word_t *packed_reads;
    edge_word_t *packed_reads_p; // current pointer
    globals.read_length_mask = (1 << bits_read_length) - 1;
    globals.capacity = std::min(max_num_reads, int64_t(1048576)); // initial capacity 1M
    gzFile fp = strcmp(globals.input_file, "-") ? gzopen(globals.input_file, "r") : gzdopen(fileno(stdin), "r");
    kseq_t *seq = kseq_init(fp); // kseq to read files
    packed_reads_p = packed_reads = (edge_word_t*) MallocAndCheck(globals.capacity * globals.words_per_read * sizeof(edge_word_t), __FILE__, __LINE__);
    if (sdbg_builder_verbose >= 2) {
        log("[B::%s] Max read length is %d; words per read: %d\n", __func__, globals.max_read_length, globals.words_per_read);
    }

    // --- main reading loop ---
    bool stop_reading = false;
    while ((read_length = kseq_read(seq)) >= 0 && !stop_reading) {
        std::reverse(seq->seq.s, seq->seq.s + read_length);
        char *next_p = seq->seq.s;
        while (read_length > globals.kmer_k) {
            int scan_len = 0;
            while (scan_len < read_length && next_p[scan_len] != 'N') {
                ++scan_len;
            }

            if (scan_len > globals.kmer_k && scan_len <= globals.max_read_length) {
                if (num_reads >= globals.capacity) {
                    if (globals.capacity == max_num_reads) {
                        err("[B::%s WRANING] No enough memory to hold all the reads. Only the first %llu reads are kept.\n", __func__, num_reads);
                        stop_reading = true;
                        break;
                    } 
                    globals.capacity = std::min(globals.capacity * 2, max_num_reads);
                    edge_word_t *new_ptr = (edge_word_t*) realloc(packed_reads, globals.capacity * globals.words_per_read * sizeof(edge_word_t));
                    if (new_ptr != NULL) {
                        packed_reads = new_ptr;
                        packed_reads_p = packed_reads + globals.words_per_read * num_reads;
                        globals.capacity = globals.capacity;
                    } else {
                        err("[B::%s WRANING] No enough memory to hold all the reads. Only the first %llu reads are kept.\n", __func__, num_reads);
                        stop_reading = true;
                        break;
                    }
                }
                // read length is ok! compress and store the packed read
                PackReadFromAscii(next_p, packed_reads_p, scan_len, globals.words_per_read);
                packed_reads_p += globals.words_per_read;
                ++num_reads;
            } else if (scan_len > globals.max_read_length) { // this read length is wrong
                err("[B::%s WARNING] Found a read of length %d > max read length = %d\n, it will be discarded.", __func__, scan_len, globals.max_read_length);
            }

            while (scan_len < read_length && next_p[scan_len] == 'N') {
                ++scan_len;
            }
            read_length -= scan_len;
            next_p += scan_len;
        }
    }

    globals.num_reads = num_reads;
    globals.mem_packed_reads = globals.num_reads * globals.words_per_read * sizeof(edge_word_t);
    globals.packed_reads = (edge_word_t*) ReAllocAndCheck(packed_reads, globals.mem_packed_reads, __FILE__, __LINE__);
    if (!globals.packed_reads) {
        err("[B::%s ERROR] Cannot reallocate memory for packed reads!\n", __func__);
        exit(1);
    }

    kseq_destroy(seq);
    gzclose(fp);
}

/**
 * @brief init memory for bucket scan
 */
void PrepareBucketScan(struct global_data_t &globals) {
    // init read partitions
    for (int t = 0; t < globals.num_cpu_threads; ++t) {
        struct readpartition_data_t &rp = globals.readpartitions[t];
        rp.rp_id = t;
        rp.globals = &globals;
        rp.rp_bucket_sizes = (int64_t *) MallocAndCheck(phase1::kNumBuckets * sizeof(int64_t), __FILE__, __LINE__);
        rp.rp_bucket_offsets = (int64_t *) MallocAndCheck(phase1::kNumBuckets * sizeof(int64_t), __FILE__, __LINE__);
        // distribute reads to partitions
        int64_t average = globals.num_reads / globals.num_cpu_threads;
        rp.rp_start_id = t * average;
        rp.rp_end_id = t < globals.num_cpu_threads - 1 ? (t + 1) * average : globals.num_reads;
        rp.rp_lv1_differential_base = EncodeOffset(rp.rp_start_id, 0, 0, globals.offset_num_bits);
    }

    // init bucket partitions
    for (int t = 0; t < globals.num_cpu_threads - globals.phase1_num_output_threads; ++t) {
        struct bucketpartition_data_t &bp = globals.bucketpartitions[t];
        bp.bp_id = t;
        bp.globals = &globals;
    }

    globals.bucket_sizes = (int64_t *) MallocAndCheck(phase1::kNumBuckets * sizeof(int64_t), __FILE__, __LINE__);
}

void* PreprocessScanToFillBucketSizesThread(void *_data) {
    struct readpartition_data_t &rp = *((struct readpartition_data_t*) _data);
    struct global_data_t &globals = *(rp.globals);
    int64_t *bucket_sizes = rp.rp_bucket_sizes;
    memset(bucket_sizes, 0, phase1::kNumBuckets * sizeof(int64_t));
    edge_word_t *read_p = PACKED_READS(rp.rp_start_id, globals);
    KmerUint32 edge, rev_edge; // (k+1)-mer and its rc
    for (int64_t read_id = rp.rp_start_id; read_id < rp.rp_end_id; ++read_id, read_p += globals.words_per_read) {
        int read_length = GetReadLength(read_p, globals.words_per_read, globals.read_length_mask);
        if (read_length < globals.kmer_k + 1) { continue; }
        edge.init(read_p, globals.kmer_k + 1);
        rev_edge.clean();
        for (int i = 0; i <= globals.kmer_k; ++i) {
            rev_edge.Append(3 - ExtractNthChar(read_p, globals.kmer_k - i));
        }

        int last_char_offset = globals.kmer_k;
        while (true) {
            if (rev_edge < edge) {
                bucket_sizes[rev_edge.data_[0] >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar]++;
            } else {
                bucket_sizes[edge.data_[0] >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar]++;
            }

            if (++last_char_offset >= read_length) {
                break;
            } else {
                int c = ExtractNthChar(read_p, last_char_offset);
                edge.ShiftLeftAppend(c);
                rev_edge.ShiftRightAppend(3 - c);
            }
        }
    }
    return NULL;
}

// multithread
void PreprocessScanToFillBucketSizes(struct global_data_t &globals) {
    // create threads
    for (int t = 1; t < globals.num_cpu_threads; ++t) {
        pthread_create(&(globals.readpartitions[t].thread), NULL, PreprocessScanToFillBucketSizesThread, &globals.readpartitions[t]);
    }
    PreprocessScanToFillBucketSizesThread(&globals.readpartitions[0]);
    for (int t = 1; t < globals.num_cpu_threads; ++t) {
        pthread_join(globals.readpartitions[t].thread, NULL);
    }
    // sum up readpartitions bucketsizes to form global bucketsizes
    int64_t *bucket_sizes = globals.bucket_sizes;
    memset(bucket_sizes, 0, phase1::kNumBuckets * sizeof(int64_t));
    // the array accesses in this loop are optimized by the compiler??
    for (int t = 0; t < globals.num_cpu_threads; ++t) {
        for (int b = 0; b < phase1::kNumBuckets; ++b) {
            bucket_sizes[b] += globals.readpartitions[t].rp_bucket_sizes[b];
        }
    }
}

void InitGlobalData(struct global_data_t &globals) {
    // --- initialize writer ---
    for (int t = 0; t < globals.phase1_num_output_threads; ++t) {
        char edges_file_name[10240];
        sprintf(edges_file_name, "%s.edges.%d", globals.output_prefix, t);
        globals.word_writer[t].init(edges_file_name);
    }

    // --- compute offset bits ---
    {
        globals.offset_num_bits = 0;
        int len = 1;
        while (len - 1 < globals.max_read_length) {
            globals.offset_num_bits++;
            len *= 2;
        }
    }

    // --- initialize stat ---
    globals.edge_counting = (int64_t *) MallocAndCheck((kMaxMulti_t + 1) * sizeof(int64_t), __FILE__, __LINE__);
    globals.thread_edge_counting = (int64_t *) MallocAndCheck((kMaxMulti_t + 1) * globals.phase1_num_output_threads * sizeof(int64_t), __FILE__, __LINE__);
    memset(globals.edge_counting, 0, (kMaxMulti_t + 1) * sizeof(int64_t));

    // --- Fill bucket size ---
    xtimer_t timer;
    timer.reset();
    timer.start();
    if (sdbg_builder_verbose >= 3) {
        log("[B::%s] Filling read partition buckets...\n", __func__);
    }
    pthread_mutex_init(&globals.lv1_items_scanning_lock, NULL); // init lock
    PrepareBucketScan(globals);
    PreprocessScanToFillBucketSizes(globals); // Multithread: fill the read partition buckets, then sum up into the global buckets
    globals.max_bucket_size = *std::max_element(globals.bucket_sizes, globals.bucket_sizes + phase1::kNumBuckets);
    globals.tot_bucket_size = 0;
    for (int i = 0; i < phase1::kNumBuckets; ++i) { globals.tot_bucket_size += globals.bucket_sizes[i]; }
    timer.stop();
    if (sdbg_builder_verbose >= 3) {
        log("[B::%s] Done. Time elapsed: %.4lfs\n", __func__, timer.elapsed());
    }

    // --- calculate lv2 memory ---
#ifdef DISABLE_GPU
    globals.max_lv2_items = std::max(globals.max_bucket_size, kMinLv2BatchSize);
#else
    int64_t lv2_mem = globals.gpu_mem - 1073741824; // should reserver ~1G for GPU sorting
    globals.max_lv2_items = lv2_mem / GPU_BYTES_PER_ITEM;
    globals.max_lv2_items = std::min(globals.max_lv2_items, globals.tot_bucket_size);
    if (globals.max_bucket_size > globals.max_lv2_items) {
        err("[ERROR B::%s] Bucket too large for GPU: contains %lld items. Please try CPU version.\n", __func__, globals.max_bucket_size);
        // TODO: auto switch to CPU version
        exit(1);
    }
#endif
    globals.words_per_substring = DivCeiling((globals.kmer_k + 1) * kBitsPerEdgeChar, kBitsPerEdgeWord);
    globals.words_per_edge = DivCeiling((globals.kmer_k + 1) * kBitsPerEdgeChar + kBitsPerMulti_t, kBitsPerEdgeWord);
    // lv2 bytes: substring, permutation, readinfo
    int64_t lv2_bytes_per_item = (globals.words_per_substring) * sizeof(edge_word_t) + sizeof(uint32_t) + sizeof(int64_t);
    lv2_bytes_per_item = lv2_bytes_per_item * 2; // double buffering
#ifdef DISABLE_GPU
    lv2_bytes_per_item += sizeof(uint64_t) * 2; // CPU memory is used to simulate GPU
#endif

    if (sdbg_builder_verbose >= 2) {
        log("[B::%s] %d words per read, %d words per substring, %d words per edge\n", __func__, globals.words_per_read, globals.words_per_substring, globals.words_per_edge);
    }
    // --- memory stuff ---
    int64_t mem_remained = globals.host_mem 
                         - globals.mem_packed_reads
                         - globals.num_reads * sizeof(unsigned char) * 2 // first_in0 & last_out0
                         - phase1::kNumBuckets * sizeof(int64_t) * (globals.num_cpu_threads * 3 + 1)
                         - (kMaxMulti_t + 1) * (globals.phase1_num_output_threads + 1) * sizeof(int64_t);
    if (globals.mem_flag == 1) {
        // auto set memory
        globals.max_lv1_items = std::max(globals.max_lv2_items, int64_t(globals.tot_bucket_size / (kDefaultLv1ScanTime - 0.5)));
        int64_t mem_needed = globals.max_lv1_items * LV1_BYTES_PER_ITEM + globals.max_lv2_items * lv2_bytes_per_item;
        if (mem_needed > mem_remained) {
            GetNumItems(globals, mem_remained, lv2_bytes_per_item);
        }
    } else if (globals.mem_flag == 0) {
        // min memory
        globals.max_lv1_items = std::max(globals.max_lv2_items, int64_t(globals.tot_bucket_size / (kMaxLv1ScanTime - 0.5)));
        int64_t mem_needed = globals.max_lv1_items * LV1_BYTES_PER_ITEM + globals.max_lv2_items * lv2_bytes_per_item;
        if (mem_needed > mem_remained) {
            GetNumItems(globals, mem_remained, lv2_bytes_per_item);
        } else {
            GetNumItems(globals, mem_needed, lv2_bytes_per_item);
        }
    } else {
        // use all
        GetNumItems(globals, mem_remained, lv2_bytes_per_item);
    }

    if (sdbg_builder_verbose >= 2) {
        log("[B::%s] Memory for reads: %lld\n", __func__, globals.mem_packed_reads);
        log("[B::%s] max # lv.1 items = %lld\n", __func__, globals.max_lv1_items);
        log("[B::%s] max # lv.2 items = %lld\n", __func__, globals.max_lv2_items);
    }

    // --- alloc memory ---
    globals.lv1_items = (int*) MallocAndCheck(globals.max_lv1_items * sizeof(int), __FILE__, __LINE__);
    globals.lv2_substrings = (edge_word_t*) MallocAndCheck(globals.max_lv2_items * globals.words_per_substring * sizeof(edge_word_t), __FILE__, __LINE__);
    globals.permutation = (uint32_t *) MallocAndCheck(globals.max_lv2_items * sizeof(uint32_t), __FILE__, __LINE__);
    globals.lv2_substrings_to_output = (edge_word_t*) MallocAndCheck(globals.max_lv2_items * globals.words_per_substring * sizeof(edge_word_t), __FILE__, __LINE__);
    globals.permutation_to_output = (uint32_t *) MallocAndCheck(globals.max_lv2_items * sizeof(uint32_t), __FILE__, __LINE__);
    globals.lv2_read_info = (int64_t *) MallocAndCheck(globals.max_lv2_items * sizeof(int64_t), __FILE__, __LINE__);
    globals.lv2_read_info_to_output = (int64_t *) MallocAndCheck(globals.max_lv2_items * sizeof(int64_t), __FILE__, __LINE__);
#ifdef LONG_READS
    globals.first_0_out = (uint16_t*) MallocAndCheck(globals.num_reads * sizeof(uint16_t), __FILE__, __LINE__);
    globals.last_0_in = (uint16_t*) MallocAndCheck(globals.num_reads * sizeof(uint16_t), __FILE__, __LINE__);
#else
    globals.first_0_out = (unsigned char*) MallocAndCheck(globals.num_reads * sizeof(unsigned char), __FILE__, __LINE__);
    globals.last_0_in = (unsigned char*) MallocAndCheck(globals.num_reads * sizeof(unsigned char), __FILE__, __LINE__);
#endif
#ifdef DISABLE_GPU
    globals.cpu_sort_space = (uint64_t*) MallocAndCheck(sizeof(uint64_t) * globals.max_lv2_items, __FILE__, __LINE__);
#endif
    memset(globals.first_0_out, 0xFF, globals.num_reads * sizeof(globals.first_0_out[0]));
    memset(globals.last_0_in, 0xFF, globals.num_reads * sizeof(globals.last_0_in[0]));

    // --- write the edge file header ---
    globals.word_writer[0].output(globals.kmer_k);
    globals.word_writer[0].output(globals.words_per_edge);
}

/**
 * @brief worker thread for Lv1ScanToFillOffests
 */
void* Lv1ScanToFillOffsetsThread(void *_data) {
    struct readpartition_data_t &rp = *((struct readpartition_data_t*) _data);
    struct global_data_t &globals = *(rp.globals);
    int64_t *prev_full_offsets = (int64_t *)MallocAndCheck(phase1::kNumBuckets * sizeof(int64_t), __FILE__, __LINE__); // temporary array for computing differentials
    assert(prev_full_offsets != NULL);
    for (int b = globals.lv1_start_bucket; b < globals.lv1_end_bucket; ++b)
        prev_full_offsets[b] = rp.rp_lv1_differential_base;
    // this loop is VERY similar to that in PreprocessScanToFillBucketSizesThread
    edge_word_t *read_p = PACKED_READS(rp.rp_start_id, globals);
    KmerUint32 edge, rev_edge; // (k+1)-mer and its rc
    int key;
    for (int64_t read_id = rp.rp_start_id; read_id < rp.rp_end_id; ++read_id, read_p += globals.words_per_read) {
        int read_length = GetReadLength(read_p, globals.words_per_read, globals.read_length_mask);
        if (read_length < globals.kmer_k + 1) { continue; }
        edge.init(read_p, globals.kmer_k + 1);
        rev_edge.clean();
        for (int i = 0; i <= globals.kmer_k; ++i) {
            rev_edge.Append(3 - ExtractNthChar(read_p, globals.kmer_k - i));
        }

        // ===== this is a macro to save some copy&paste ================
#define CHECK_AND_SAVE_OFFSET(offset, strand)                                   \
    do {                                                                \
      assert(offset + globals.kmer_k < read_length); \
      if (((key - globals.lv1_start_bucket) ^ (key - globals.lv1_end_bucket)) & kSignBitMask) { \
        int64_t full_offset = EncodeOffset(read_id, offset, strand, globals.offset_num_bits); \
        int64_t differential = full_offset - prev_full_offsets[key];      \
        if (differential > kDifferentialLimit) {                      \
          pthread_mutex_lock(&globals.lv1_items_scanning_lock); \
          globals.lv1_items[ rp.rp_bucket_offsets[key]++ ] = -globals.lv1_items_special.size() - 1; \
          globals.lv1_items_special.push_back(full_offset);                  \
          pthread_mutex_unlock(&globals.lv1_items_scanning_lock); \
        } else {                                                              \
          assert((int) differential >= 0); \
          globals.lv1_items[ rp.rp_bucket_offsets[key]++ ] = (int) differential; \
        } \
        prev_full_offsets[key] = full_offset;                           \
      }                                                                 \
    } while (0)
        // ^^^^^ why is the macro surrounded by a do-while? please ask Google
        // =========== end macro ==========================

        // shift the key char by char
        int last_char_offset = globals.kmer_k;
        while (true) {
            if (rev_edge < edge) {
                key = rev_edge.data_[0] >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar;
                CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 1);
            } else {
                key = edge.data_[0] >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar;
                CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 0);
            }

            if (++last_char_offset >= read_length) {
                break;
            } else {
                int c = ExtractNthChar(read_p, last_char_offset);
                edge.ShiftLeftAppend(c);
                rev_edge.ShiftRightAppend(3 - c);
            }
        }
    }

#undef CHECK_AND_SAVE_OFFSET

    free(prev_full_offsets);
    return NULL;
}

/**
 * @brief file LV1 items into host mem
 */
void Lv1ScanToFillOffests(struct global_data_t &globals) {
    globals.lv1_items_special.clear();
    Lv1ComputeBucketOffset(globals);
    // create threads
    for (int t = 1; t < globals.num_cpu_threads; ++t) {
        pthread_create(&(globals.readpartitions[t].thread), NULL, Lv1ScanToFillOffsetsThread, &globals.readpartitions[t]);
    }
    Lv1ScanToFillOffsetsThread(&globals.readpartitions[0]);
    for (int t = 1; t < globals.num_cpu_threads; ++t) {
        pthread_join(globals.readpartitions[t].thread, NULL);
    }
    // revert rp_bucket_offsets
    Lv1ComputeBucketOffset(globals);
}

// single thread helper function
// 'spacing' is the strip length for read-word "coalescing"
void CopySubstring(edge_word_t* dest, edge_word_t* src_read, int offset, int num_chars_to_copy, global_data_t &globals) {
    int64_t spacing = globals.lv2_num_items;
    int words_per_read = globals.words_per_read;
    int words_per_substring = globals.words_per_substring;

    // copy words of the suffix to the suffix pool
    int which_word = offset / kCharsPerEdgeWord;
    int word_offset = offset % kCharsPerEdgeWord;
    edge_word_t *src_p = src_read + which_word;
    edge_word_t *dest_p = dest;
    int num_words_copied = 0;
    if (!word_offset) { // special case (word aligned), easy
        while (which_word < words_per_read && num_words_copied < words_per_substring) {
            *dest_p = *src_p; // write out
            dest_p += spacing;
            src_p++;
            which_word++;
            num_words_copied++;
        }
    } else { // not word-aligned
        int bit_shift = offset * kBitsPerEdgeChar;
        edge_word_t s = *src_p;
        edge_word_t d = s << bit_shift;
        which_word++;
        while (which_word < words_per_read) {
            s = *(++src_p);
            d |= s >> (kBitsPerEdgeWord - bit_shift);
            *dest_p = d; // write out
            if (++num_words_copied >= words_per_substring) goto here;
            dest_p += spacing;
            d = s << bit_shift;
            which_word++;
        }
        *dest_p = d; // write last word
here:
        ;
    }

    {
        // now mask the extra bits (TODO can be optimized)
        int num_bits_to_copy = num_chars_to_copy * 2;
        int which_word = num_bits_to_copy / kBitsPerEdgeWord;
        edge_word_t *p = dest + which_word * spacing;
        int bits_to_clear = kBitsPerEdgeWord - num_bits_to_copy % kBitsPerEdgeWord;
        if (bits_to_clear < kBitsPerEdgeWord) {
            *p >>= bits_to_clear;
            *p <<= bits_to_clear;
        } else if (which_word < globals.words_per_substring) {
            *p = 0;
        }
        which_word++;
        while (which_word < globals.words_per_substring) { // fill zero
            *(p+=spacing) = 0;
            which_word++;
        }
    }
}

void CopySubstringRC(edge_word_t* dest, edge_word_t* src_read, int offset, int num_chars_to_copy, global_data_t &globals) {
    assert(num_chars_to_copy == globals.kmer_k + 1);
    int spacing = globals.lv2_num_items;
    int which_word = (offset + num_chars_to_copy - 1) / kCharsPerEdgeWord;
    int word_offset = (offset + num_chars_to_copy - 1) % kCharsPerEdgeWord;
    edge_word_t *dest_p = dest;

    if (word_offset == kCharsPerEdgeWord - 1) { // edge_word_t aligned
        for (int i = 0; i < globals.words_per_substring && i <= which_word; ++i) {
            *dest_p = ~ mirror(src_read[which_word - i]);
            dest_p += spacing;
        }
    } else {
        int bit_offset = (kCharsPerEdgeWord - 1 - word_offset) * kBitsPerEdgeChar;
        int i;
        edge_word_t w;
        for (i = 0; i < globals.words_per_substring - 1 && i < which_word; ++i) {
            w = (src_read[which_word - i] >> bit_offset) |
                                      (src_read[which_word - i - 1] << (kBitsPerEdgeWord - bit_offset));
            *dest_p = ~ mirror(w);
            dest_p += spacing;
        }
        // last word
        w = src_read[which_word - i] >> bit_offset;
        if (which_word >= i + 1) {
            w |= (src_read[which_word - i - 1] << (kBitsPerEdgeWord - bit_offset));
        }
        *dest_p = ~ mirror(w);
    }

    {
        // now mask the extra bits (TODO can be optimized)
        int num_bits_to_copy = num_chars_to_copy * 2;
        int which_word = num_bits_to_copy / kBitsPerEdgeWord;
        edge_word_t *p = dest + which_word * spacing;
        int bits_to_clear = kBitsPerEdgeWord - num_bits_to_copy % kBitsPerEdgeWord;
        if (bits_to_clear < kBitsPerEdgeWord) {
            *p >>= bits_to_clear;
            *p <<= bits_to_clear;
        } else if (which_word < globals.words_per_substring) {
            *p = 0;
        }
        which_word++;
        while (which_word < globals.words_per_substring) { // fill zero
            *(p+=spacing) = 0;
            which_word++;
        }
    }
}

/**
 * @brief worker thread for Lv2ExtractSubstrings
 */
void* Lv2ExtractSubstringsThread(void* _data) {
    struct bucketpartition_data_t &bp = *((struct bucketpartition_data_t*) _data);
    struct global_data_t &globals = *(bp.globals);
    int *lv1_p = globals.lv1_items + globals.readpartitions[0].rp_bucket_offsets[ bp.bp_start_bucket ];
    int64_t offset_mask = (1 << globals.offset_num_bits) - 1; // 0000....00011..11
    edge_word_t *substrings_p = globals.lv2_substrings +
                         (globals.readpartitions[0].rp_bucket_offsets[ bp.bp_start_bucket ] - globals.readpartitions[0].rp_bucket_offsets[ globals.lv2_start_bucket ]);
    int64_t *read_info_p = globals.lv2_read_info + 
                       (globals.readpartitions[0].rp_bucket_offsets[ bp.bp_start_bucket ] - globals.readpartitions[0].rp_bucket_offsets[ globals.lv2_start_bucket ]);
    for (int b = bp.bp_start_bucket; b < bp.bp_end_bucket; ++b) {
        for (int t = 0; t < globals.num_cpu_threads; ++t) {
            int64_t full_offset = globals.readpartitions[t].rp_lv1_differential_base;
            int num = globals.readpartitions[t].rp_bucket_sizes[b];
            for (int i = 0; i < num; ++i) {
                if (*lv1_p >= 0) {
                    full_offset += *(lv1_p++);
                } else {
                    full_offset = globals.lv1_items_special[-1 - *(lv1_p++)];
                }
                int64_t read_id = full_offset >> (globals.offset_num_bits + 1);
                int strand = full_offset & 1;
                int offset = (full_offset >> 1) & offset_mask;
                int num_chars_to_copy = globals.kmer_k + 1;
                unsigned char prev, next;
                if (offset > 0) {
                    prev = ExtractNthChar(PACKED_READS(read_id, globals), offset - 1);
                } else {
                    prev = kSentinelValue;
                }

                if (offset + globals.kmer_k + 1 < GetReadLengthByID(read_id, globals)) {
                    next = ExtractNthChar(PACKED_READS(read_id, globals), offset + globals.kmer_k + 1);
                } else {
                    next = kSentinelValue;
                }

                if (strand == 0) {
                    CopySubstring(substrings_p, PACKED_READS(read_id, globals), offset, num_chars_to_copy, globals);
                    *read_info_p = (full_offset << 6) | (prev << 3) | next;
                } else {
                    CopySubstringRC(substrings_p, PACKED_READS(read_id, globals), offset, num_chars_to_copy, globals);
                    *read_info_p = (full_offset << 6) | ((next == kSentinelValue ? kSentinelValue : (3 - next)) << 3)
                                                      | (prev == kSentinelValue ? kSentinelValue : (3 - prev));
                }
#ifdef DBJ_DEBUG
                if ((*substrings_p >> (32 - phase2::kBucketPrefixLength * 2)) != b) {
                    debug("WRONG substring wrong:%d right:%d read_id:%lld offset:%d strand: %d num_chars_to_copy:%d\n", *substrings_p >> (32 - phase2::kBucketPrefixLength * 2), b, read_id, offset, strand, num_chars_to_copy);
                }
#endif
                substrings_p++;
                read_info_p++;
            }
        }
    }
    return NULL;
}

/**
 * @brief extract true substring used for sorting
 */
void Lv2ExtractSubstrings(struct global_data_t &globals) {
    Lv2DistributeBucketPartitions(globals, globals.phase1_num_output_threads);
    // create threads
    for (int t = 0; t < globals.num_cpu_threads-globals.phase1_num_output_threads; ++t) {
        pthread_create(&(globals.bucketpartitions[t].thread), NULL, Lv2ExtractSubstringsThread, &globals.bucketpartitions[t]);
    }
    for (int t = 0; t < globals.num_cpu_threads-globals.phase1_num_output_threads; ++t) {
        pthread_join(globals.bucketpartitions[t].thread, NULL);
    }
#ifdef DBJ_DEBUG
    for (int i = 0; i < globals.lv2_num_items; ++i) {
        edge_word_t *substrings_p = globals.lv2_substrings + i;
        int cur_bucket = *substrings_p >> (32 - phase1::kBucketPrefixLength * 2);
        assert(cur_bucket < globals.lv2_end_bucket && cur_bucket >= globals.lv2_start_bucket);
    }
#endif
}

// helper function for counting
inline bool IsDifferentEdges(edge_word_t *item1, edge_word_t* item2, int num_words, int spacing) {
    for (int i = num_words - 1; i >= 0; --i) {
        if (*(item1 + (int64_t)i * spacing) != *(item2 + (int64_t)i * spacing)) {
            return true;
        }
    }
    return false;
}

/**
 * @brief pack an edge and its multiplicity to word-aligned spaces
 */
inline void PackEdge(edge_word_t *dest, edge_word_t *item, int counting, struct global_data_t &globals) {
    for (int i = 0; i < globals.words_per_edge && i < globals.words_per_substring; ++i) {
        dest[i] = *(item + (int64_t)i * globals.lv2_num_items_to_output);
    }
    int chars_in_last_word = (globals.kmer_k + 1) % kCharsPerEdgeWord;
    int which_word = (globals.kmer_k + 1) / kCharsPerEdgeWord;
    if (chars_in_last_word > 0) {
        dest[which_word] >>= (kCharsPerEdgeWord - chars_in_last_word) * 2;
        dest[which_word] <<= (kCharsPerEdgeWord - chars_in_last_word) * 2;
    } else {
        dest[which_word] = 0;
    }
    while (++which_word < globals.words_per_edge) {
        dest[which_word] = 0;
    }

    dest[globals.words_per_edge - 1] |= std::min(kMaxMulti_t, counting);
}

void* Lv2CountingThread(void *_op) {
    struct outputpartition_data_t *op = (struct outputpartition_data_t*) _op;
    struct global_data_t &globals = *(op->globals);
    int64_t op_start_index = op->op_start_index;
    int64_t op_end_index = op->op_end_index;
    int thread_id = op->op_id;
    xtimer_t local_timer;
    local_timer.start();
    local_timer.reset();
    int start_idx;
    int end_idx;
    edge_word_t packed_edge[32];
    int count_prev[5], count_next[5];
    int64_t offset_mask = (1 << globals.offset_num_bits) - 1;
    int64_t *thread_edge_counting = globals.thread_edge_counting + thread_id * (kMaxMulti_t + 1);

    for (int i = op_start_index; i < op_end_index; i = end_idx) {
        start_idx = i;
        end_idx = i + 1;
        edge_word_t *first_item = globals.lv2_substrings_to_output + (globals.permutation_to_output[i]);
        while (end_idx < op_end_index) {
            if (IsDifferentEdges(first_item,
                               globals.lv2_substrings_to_output + globals.permutation_to_output[end_idx],
                               globals.words_per_substring, globals.lv2_num_items_to_output)) {
                break;
            }
            ++end_idx;
        }
        int count = end_idx - start_idx;

        // update read's first and last

        memset(count_prev, 0, sizeof(int) * 4);
        memset(count_next, 0, sizeof(int) * 4);
        bool has_in = false;
        bool has_out = false;
        for (int j = start_idx; j < end_idx; ++j) {
            int prev_and_next = globals.lv2_read_info_to_output[globals.permutation_to_output[j]] & ((1 << 6) - 1);
            count_prev[prev_and_next >> 3]++;
            count_next[prev_and_next & 7]++;
        }

        for (int j = 0; j < 4; ++j) {
            if (count_prev[j] >= globals.kmer_freq_threshold) { has_in = true; }
            if (count_next[j] >= globals.kmer_freq_threshold) { has_out = true; }
        }

        if (!has_in && count >= globals.kmer_freq_threshold) {
            for (int j = start_idx; j < end_idx; ++j) {
                int64_t read_info = globals.lv2_read_info_to_output[globals.permutation_to_output[j]] >> 6;
                int strand = read_info & 1;
                int offset = (read_info >> 1) & offset_mask;
                int64_t read_id = read_info >> (1 + globals.offset_num_bits);

                if (strand == 0) {
                    // update last
                    while (true) {
                        auto old_value = globals.last_0_in[read_id];
                        if (old_value != kSentinelOffset && old_value >= offset) { break; }
                        if (__sync_bool_compare_and_swap(globals.last_0_in + read_id, old_value, offset)) {
                            break;
                        }
                    }
                } else {
                    // update first
                    offset++;
                    while (true) {
                        auto old_value = globals.first_0_out[read_id];
                        if (old_value <= offset) { break; }
                        if (__sync_bool_compare_and_swap(globals.first_0_out + read_id, old_value, offset)) {
                            break;
                        }
                    }
                }
            }
        }

        if (!has_out && count >= globals.kmer_freq_threshold) {
            for (int j = start_idx; j < end_idx; ++j) {
                int64_t read_info = globals.lv2_read_info_to_output[globals.permutation_to_output[j]] >> 6;
                int strand = read_info & 1;
                int offset = (read_info >> 1) & offset_mask;
                int64_t read_id = read_info >> (1 + globals.offset_num_bits);

                if (strand == 0) {
                    // update first
                    offset++;
                    while (true) {
                        auto old_value = globals.first_0_out[read_id];
                        if (old_value <= offset) { break; }
                        if (__sync_bool_compare_and_swap(globals.first_0_out + read_id, old_value, offset)) {
                            break;
                        }
                    }
                } else {
                    // update last
                    while (true) {
                        auto old_value = globals.last_0_in[read_id];
                        if (old_value != kSentinelOffset && old_value >= offset) { break; }
                        if (__sync_bool_compare_and_swap(globals.last_0_in + read_id, old_value, offset)) {
                            break;
                        }
                    }
                }
            }
        }

        ++thread_edge_counting[std::min(count, kMaxMulti_t)];
        if (count >= globals.kmer_freq_threshold) {
            PackEdge(packed_edge, first_item, count, globals);
            for (int x = 0; x < globals.words_per_edge; ++x) {
                globals.word_writer[thread_id].output(packed_edge[x]);
            }
        }
    }
    local_timer.stop();

    if (sdbg_builder_verbose >= 4) {
        log("[B::%s] Counting time elapsed: %.4lfs\n", __func__, local_timer.elapsed());
    }
    return NULL;
}

/**
 * @brief count and output solid (k+1)-mer
 */
void Lv2Counting(struct global_data_t &globals) {
    int64_t last_end_index = 0;
    int64_t items_per_thread = globals.lv2_num_items_to_output / globals.phase1_num_output_threads;

    for (int thread_id = 0; thread_id < globals.phase1_num_output_threads - 1; ++thread_id) {
        int64_t this_start_index = last_end_index;
        int64_t this_end_index = this_start_index + items_per_thread;
        if (this_end_index > globals.lv2_num_items_to_output) { this_end_index = globals.lv2_num_items_to_output; }
        if (this_end_index > 0) {
            while (this_end_index < globals.lv2_num_items_to_output) {
                edge_word_t *prev_item = globals.lv2_substrings_to_output + (globals.permutation_to_output[this_end_index - 1]);
                edge_word_t *item = globals.lv2_substrings_to_output + (globals.permutation_to_output[this_end_index]);
                if (IsDifferentEdges(prev_item, item, globals.words_per_substring, globals.lv2_num_items_to_output)) {
                    break;
                }
                ++this_end_index;
            }
        }
        globals.outputpartitions[thread_id].op_start_index = this_start_index;
        globals.outputpartitions[thread_id].op_end_index = this_end_index;
        last_end_index = this_end_index;
    }

    // last partition
    globals.outputpartitions[globals.phase1_num_output_threads - 1].op_start_index = last_end_index;
    globals.outputpartitions[globals.phase1_num_output_threads - 1].op_end_index = globals.lv2_num_items_to_output;

    memset(globals.thread_edge_counting, 0, sizeof(int64_t) * (kMaxMulti_t + 1) * globals.phase1_num_output_threads);
    for (int thread_id = 0; thread_id < globals.phase1_num_output_threads; ++thread_id) {
        globals.outputpartitions[thread_id].op_id = thread_id;
        globals.outputpartitions[thread_id].globals = &globals;
        pthread_create(&globals.output_threads[thread_id], NULL, Lv2CountingThread, &globals.outputpartitions[thread_id]);
    }
}

void Lv2CountingJoin(struct global_data_t &globals) {
    for (int thread_id = 0; thread_id < globals.phase1_num_output_threads; ++thread_id) {
        pthread_join(globals.output_threads[thread_id], NULL);
        for (int i = 1; i <= kMaxMulti_t; ++i) {
            globals.edge_counting[i] += globals.thread_edge_counting[thread_id * (kMaxMulti_t + 1) + i];
        }
    }
}

void Phase1Clean(struct global_data_t &globals) {
    pthread_mutex_destroy(&globals.lv1_items_scanning_lock);
    free(globals.packed_reads);
    free(globals.bucket_sizes);
    free(globals.lv1_items);
    free(globals.lv2_substrings);
    free(globals.permutation);
    free(globals.lv2_substrings_to_output);
    free(globals.permutation_to_output);
    free(globals.lv2_read_info);
    free(globals.lv2_read_info_to_output);
    free(globals.first_0_out);
    free(globals.last_0_in);
    free(globals.edge_counting);
    free(globals.thread_edge_counting);
    for (int t = 0; t < globals.num_cpu_threads; ++t) {
        free(globals.readpartitions[t].rp_bucket_sizes);
        free(globals.readpartitions[t].rp_bucket_offsets);
    }
    for (int t = 0; t < globals.phase1_num_output_threads; ++t) {        
        globals.word_writer[t].destroy();
    }

#ifdef DISABLE_GPU
    free(globals.cpu_sort_space);
#endif
}

void Phase1Entry(struct global_data_t &globals) {
    xtimer_t timer;
    // --- read queries ---
    timer.reset();
    timer.start();

    if (sdbg_builder_verbose >= 2) {
        log("[B::%s] Reading input...\n", __func__);
    }
    InitDNAMap();
    ReadInputFile(globals);

    if (sdbg_builder_verbose >= 2) {
        log("[B::%s] Done reading input, %lld reads in total.\n", __func__, globals.num_reads);
    }
    timer.stop();
    if (sdbg_builder_verbose >= 2) {
        log("[B::%s] Time elapsed: %.4lfs\n", __func__, timer.elapsed());
    }

#ifdef DBJ_DEBUG
    debug("The 5th read is [[");
    DumpPackedReads(globals.packed_reads + 5 * globals.words_per_read, GetReadLengthByID(5, globals));
    debug("]]\n");
    debug("The last read is [[");
    DumpPackedReads(globals.packed_reads + (globals.num_reads-1) * globals.words_per_read, GetReadLengthByID(globals.num_reads-1, globals));
    debug("]]\n");
#endif

    // --- init globals ---
    InitGlobalData(globals); // must be called AFTER reading input

    ////////////////////////////////// Start processing... ////////////////////////////

    //======================================== LEVEL 1 loop ============================================//
    int lv1_iteration = 0;
    globals.lv1_start_bucket = 0;
    timer.reset();
    timer.start();
    bool output_thread_created = false;
    while (globals.lv1_start_bucket < phase1::kNumBuckets) {
        xtimer_t local_timer;
        lv1_iteration++;

        // --- finds the bucket range for this iteration ---
        local_timer.reset();
        local_timer.start();
        globals.lv1_end_bucket = FindEndBucket(globals.bucket_sizes, globals.lv1_start_bucket, phase1::kNumBuckets, globals.max_lv1_items, globals.lv1_num_items);
        if (globals.lv1_num_items == 0) { // i.e. can't even hold a single bucket (very unlikely though)
            err("[ERROR B::%s] Bucket %d too large for lv.1: contains %lld items.\n", __func__, globals.lv1_end_bucket, globals.bucket_sizes[globals.lv1_end_bucket]);
            exit(1);
        }

        // --- LV1 scan to fill offsets in host mem ---

        if (sdbg_builder_verbose >= 3) {
            log("[B::%s] Iteration %d, from bucket %d to %d\n", __func__, lv1_iteration, globals.lv1_start_bucket, globals.lv1_end_bucket-1);
            log("[B::%s] Scanning and filling offsets...\n", __func__);
        }
        Lv1ScanToFillOffests(globals);
        local_timer.stop();

        if (sdbg_builder_verbose >= 3) {
            log("[B::%s] Scanning and filling offsets...done. Number of large differentials: %llu\n", __func__, globals.lv1_items_special.size());
            log("[B::%s] time elapsed: %.4lfs\n", __func__, local_timer.elapsed());
        }

        if (globals.lv1_items_special.size() > kDifferentialLimit) {
            err("[ERROR B::%s] Too many large differentials!\n", __func__);
            exit(1);
        }

        //======================================== LEVEL 2 loop ==========================================//
        int lv2_iteration = 0;
        globals.lv2_start_bucket = globals.lv1_start_bucket;
        while (globals.lv2_start_bucket < globals.lv1_end_bucket) {
            lv2_iteration++;
            // --- finds the bucket range for this iteration ---
            local_timer.reset();
            local_timer.start();
            globals.lv2_end_bucket = FindEndBucket(globals.bucket_sizes, globals.lv2_start_bucket, globals.lv1_end_bucket, globals.max_lv2_items, globals.lv2_num_items);
            if (globals.lv2_num_items == 0) { // i.e. can't even hold a single bucket
                err("[ERROR B::%s] Bucket %d too large for lv.2: contains %lld items\n", __func__, globals.lv2_end_bucket, globals.bucket_sizes[globals.lv2_end_bucket]);
                exit(1);
            }

            // --- extract substring to host mem ---

            if (sdbg_builder_verbose >= 4) {
                log("[B::%s] > Iteration [%d,%d], from bucket %d to %d\n", __func__, lv1_iteration, lv2_iteration, globals.lv2_start_bucket, globals.lv2_end_bucket-1);
            }
            Lv2ExtractSubstrings(globals);
            local_timer.stop();

            if (sdbg_builder_verbose >= 4) {
                log("[B::%s] Extracting substrings...done. Time elapsed: %.4lfs\n", __func__, local_timer.elapsed());
            }
            // --- sorting ---
#ifdef DISABLE_GPU
            omp_set_num_threads(globals.num_cpu_threads - globals.phase1_num_output_threads);
            local_timer.reset();
            local_timer.start();
            lv2_cpu_sort(globals.lv2_substrings, globals.permutation, globals.cpu_sort_space, globals.words_per_substring, globals.lv2_num_items);
            omp_set_num_threads(globals.num_cpu_threads);
            local_timer.stop();

            if (sdbg_builder_verbose >= 4) {
                log("[B::%s] Sorting substrings with CPU...done. Time elapsed: %.4lfs\n", __func__, local_timer.elapsed());
            }
#else
            local_timer.reset();
            local_timer.start();
            lv2_gpu_sort(globals.lv2_substrings, globals.permutation, globals.words_per_substring, globals.lv2_num_items);
            local_timer.stop();

            if (sdbg_builder_verbose >= 4) {
                log("[B::%s] Sorting substrings with GPU...done. Time elapsed: %.4lfs\n", __func__, local_timer.elapsed());
            }
#endif
            // --- output is pipelined ---
            if (output_thread_created) {
                Lv2CountingJoin(globals);
            }

            globals.lv2_num_items_to_output = globals.lv2_num_items;
            std::swap(globals.lv2_substrings_to_output, globals.lv2_substrings);
            std::swap(globals.permutation_to_output, globals.permutation);
            std::swap(globals.lv2_read_info_to_output, globals.lv2_read_info);
            Lv2Counting(globals);
            output_thread_created = true;

            globals.lv2_start_bucket = globals.lv2_end_bucket;
        } // end LEVEL 2 loop

        globals.lv1_start_bucket = globals.lv1_end_bucket;
    } // end LEVEL 1 loop

    if (output_thread_created) {
        Lv2CountingJoin(globals);
    }

    timer.stop();

    if (sdbg_builder_verbose >= 3) {
        log("[B::%s] Done all counting. Time elapsed: %.4lf\n", __func__, timer.elapsed());
    }

    // --- output reads for mercy ---
    int64_t num_candidate_reads = 0;
    int64_t num_has_tips = 0;
    FILE *candidate_file = OpenFileAndCheck((std::string(globals.output_prefix) + ".cand").c_str(), "wb");
    for (int64_t i = 0; i < globals.num_reads; ++i) {
        auto first = globals.first_0_out[i];
        auto last = globals.last_0_in[i];
        if (first != kSentinelOffset && last != kSentinelOffset) {
            ++num_has_tips;
            if (last > first) {
                ++num_candidate_reads;
                fwrite(PACKED_READS(i, globals), sizeof(uint32_t), globals.words_per_read, candidate_file);   
            }
        }
    }
    fclose(candidate_file);

    if (sdbg_builder_verbose >= 2) {
        log("[B::%s] Total number of candidate reads: %lld(%lld)\n", __func__, num_candidate_reads, num_has_tips);
    }

    // --- stat ---
    int64_t num_solid_edges = 0;
    for (int i = globals.kmer_freq_threshold; i <= kMaxMulti_t; ++i) {
        num_solid_edges += globals.edge_counting[i];
    }
    globals.num_edges = globals.num_outgoing_zero_nodes + globals.num_incoming_zero_nodes + num_solid_edges;

    if (sdbg_builder_verbose >= 2) {
        log("[B::%s] Total number of solid edges: %llu\n", __func__, num_solid_edges);
    }

    FILE *counting_file = OpenFileAndCheck((std::string(globals.output_prefix)+".counting").c_str(), "w");
    for (int64_t i = 1, acc = 0; i <= kMaxMulti_t; ++i) {
        acc += globals.edge_counting[i];
        fprintf(counting_file, "%lld %lld\n", (long long)i, (long long)acc);
    }
    fclose(counting_file);
    Phase1Clean(globals);
}

}

namespace phase2 {
/**
 * @brief encode edge id and offset using int64_t
 */
inline int64_t EncodeEdgeOffset(int64_t edge_id, int offset, int strand, int k_num_bits) {
    return (edge_id << (k_num_bits + 1)) | (strand << k_num_bits) | offset;
}

/**
 * @brief read edges from disk
 */
int64_t ReadEdges(global_data_t &globals) {
    // --- init reader ---
    EdgeReader edge_reader;
    edge_reader.init((std::string(globals.phase2_input_prefix) + ".edges").c_str(), globals.phase1_num_output_threads);
    // --- calc memory type and max_num_edge ---
    globals.kmer_k = edge_reader.kmer_k;
    globals.words_per_edge = DivCeiling(globals.kmer_k + 1, kCharsPerEdgeWord);
    int free_bits_in_edge = globals.words_per_edge * kBitsPerEdgeWord - (globals.kmer_k + 1) * kBitsPerEdgeChar;
    if (free_bits_in_edge >= kBitsPerMulti_t) {
        globals.mult_mem_type = 0;
    } else if (free_bits_in_edge >= 8) {
        globals.mult_mem_type = 1;
    } else {
        globals.mult_mem_type = 2;
    }
    int64_t bytes_per_edge = globals.words_per_edge * sizeof(edge_word_t) + globals.mult_mem_type;
    int64_t max_num_edges = globals.host_mem * 0.9 / bytes_per_edge; // TODO: more accurate

    if (sdbg_builder_verbose >= 2) {
        log("[B::%s] kmer_k: %d, words_per_edge: %d\n", __func__, globals.kmer_k, globals.words_per_edge);
        log("[B::%s] Max host mem: %ld, max number of edges can be loaded: %lld\n", __func__, globals.host_mem, (long long)max_num_edges);
    }

    // --- alloc memory for edges ---
    globals.capacity = std::min(max_num_edges, int64_t(10485760)); // 10M
    globals.packed_edges = (edge_word_t *) MallocAndCheck(sizeof(edge_word_t) * globals.words_per_edge * globals.capacity, __FILE__, __LINE__);
    if (globals.mult_mem_type == 1) {
        globals.multiplicity8 = (uint8_t*) MallocAndCheck(sizeof(uint8_t) * globals.capacity, __FILE__, __LINE__);
    } else if (globals.mult_mem_type == 2) {
        globals.multiplicity16 = (uint16_t*) MallocAndCheck(sizeof(uint16_t) * globals.capacity, __FILE__, __LINE__);
    }

    // --- read edges ---
    edge_word_t *edge_p = globals.packed_edges;
    int64_t num_edges = 0;
    while (edge_reader.NextEdge(edge_p)) {
        if (num_edges >= globals.capacity) {
            if (num_edges >= max_num_edges) {
                err("[B::%s ERROR] reach max_num_edges: %ld... No enough memory to build the graph.\n", __func__, max_num_edges);
                exit(1);
            }

            globals.capacity = std::min(max_num_edges, globals.capacity * 2);
            edge_word_t *new_packed_edge = (edge_word_t*) realloc(globals.packed_edges, sizeof(edge_word_t) * globals.words_per_edge * globals.capacity);
            if (new_packed_edge == NULL) {
                err("[B::%s ERROR] Fail to alloc %lld bytes; the setting of host_mem = %lld is too big?\n", 
                    __func__, sizeof(edge_word_t) * globals.words_per_edge * globals.capacity, globals.host_mem);
                exit(1);
            }
            globals.packed_edges = new_packed_edge;
            edge_p = new_packed_edge + num_edges * globals.words_per_edge;

            if (globals.mult_mem_type == 1) {
                uint8_t *new_multi8 = (uint8_t*) realloc(globals.multiplicity8, sizeof(uint8_t) * globals.capacity);
                if (new_packed_edge == NULL) {
                    err("[B::%s ERROR] Fail to alloc %lld bytes; the setting of host_mem = %lld is too big?\n", 
                        __func__, sizeof(edge_word_t) * sizeof(uint8_t) * globals.capacity, globals.host_mem);
                    exit(1);
                }
                globals.multiplicity8 = new_multi8;
            } else if (globals.mult_mem_type == 2) {
                uint16_t *new_multi16 = (uint16_t*) realloc(globals.multiplicity16, sizeof(uint16_t) * globals.capacity);
                if (new_packed_edge == NULL) {
                    err("[B::%s ERROR] Fail to alloc %lld bytes; the setting of host_mem = %lld is too big?\n", 
                        __func__, sizeof(edge_word_t) * sizeof(uint16_t) * globals.capacity, globals.host_mem);
                    exit(1);
                }
                globals.multiplicity16 = new_multi16;
            }
        }

        ++num_edges;
        if (globals.mult_mem_type == 1) {
            edge_p[globals.words_per_edge - 1] &= 0xFFFFFF00U;
            edge_p[globals.words_per_edge - 1] |= (edge_p[globals.words_per_edge] >> 8) & 0xFFU;
            globals.multiplicity8[num_edges - 1] = edge_p[globals.words_per_edge] & 0xFFU;
        } else if (globals.mult_mem_type == 2) {
            globals.multiplicity16[num_edges - 1] = edge_p[globals.words_per_edge] & kMaxMulti_t;
        }

        edge_p += globals.words_per_edge;
    }
    edge_reader.destroy();

    // --- realloc if no mercy ---
    if (!globals.need_mercy) {
        globals.packed_edges = (edge_word_t *) ReAllocAndCheck(globals.packed_edges, sizeof(edge_word_t) * globals.words_per_edge * num_edges, __FILE__, __LINE__);
        if (globals.mult_mem_type == 1) {
            globals.multiplicity8 = (uint8_t*) ReAllocAndCheck(globals.multiplicity8, sizeof(uint8_t) * num_edges, __FILE__, __LINE__);
        } else if (globals.mult_mem_type == 2) {
            globals.multiplicity16 = (uint16_t*) ReAllocAndCheck(globals.multiplicity16, sizeof(uint16_t) * num_edges, __FILE__, __LINE__);
        }
        globals.mem_packed_edges = bytes_per_edge * num_edges;
    }
    globals.num_edges = num_edges;

    if (sdbg_builder_verbose >= 3) {
        log("[B::%s] Number of edges: %lld\n", __func__, num_edges);
    }
    return num_edges;
}

/**
 * @brief read mercy edges from disk
 */
int64_t ReadMercyEdges(global_data_t &globals) {
    // --- init reader ---
    EdgeReader edge_reader;
    edge_reader.InitUnsorted((std::string(globals.phase2_input_prefix) + ".mercy").c_str(),
                             globals.num_cpu_threads - 1,
                             globals.kmer_k,
                             globals.words_per_edge + (globals.mult_mem_type > 0));

    // --- read mercy edges ---
    edge_word_t *edge_p = globals.packed_edges + globals.num_edges * globals.words_per_edge;
    int64_t bytes_per_edge = globals.words_per_edge * sizeof(edge_word_t) + globals.mult_mem_type;
    int64_t max_num_edges = globals.host_mem * 0.9 / bytes_per_edge; // TODO: more accurate
    int64_t num_edges = globals.num_edges;
    while (edge_reader.NextEdgeUnsorted(edge_p)) {
        if (num_edges >= globals.capacity) {
            if (num_edges >= max_num_edges) {
                err("[B::%s ERROR] reach max_num_edges: %ld... No enough memory to build the graph.\n", __func__, max_num_edges);
                num_edges = globals.num_edges; // reset
                exit(1);
            }

            globals.capacity = std::min(max_num_edges, globals.capacity * 2);
            edge_word_t *new_packed_edge = (edge_word_t*) realloc(globals.packed_edges, sizeof(edge_word_t) * globals.words_per_edge * globals.capacity);
            if (new_packed_edge == NULL) {
                err("[B::%s ERROR] Fail to alloc %lld bytes; the setting of host_mem = %lld is too big?\n", 
                    __func__, sizeof(edge_word_t) * globals.words_per_edge * globals.capacity, globals.host_mem);
                exit(1);
            }
            globals.packed_edges = new_packed_edge;
            edge_p = new_packed_edge + num_edges * globals.words_per_edge;

            if (globals.mult_mem_type == 1) {
                uint8_t *new_multi8 = (uint8_t*) realloc(globals.multiplicity8, sizeof(uint8_t) * globals.capacity);
                if (new_packed_edge == NULL) {
                    err("[B::%s ERROR] Fail to alloc %lld bytes; the setting of host_mem = %lld is too big?\n", 
                        __func__, sizeof(edge_word_t) * sizeof(uint8_t) * globals.capacity, globals.host_mem);
                    exit(1);
                }
                globals.multiplicity8 = new_multi8;
            } else if (globals.mult_mem_type == 2) {
                uint16_t *new_multi16 = (uint16_t*) realloc(globals.multiplicity16, sizeof(uint16_t) * globals.capacity);
                if (new_packed_edge == NULL) {
                    err("[B::%s ERROR] Fail to alloc %lld bytes; the setting of host_mem = %lld is too big?\n", 
                        __func__, sizeof(edge_word_t) * sizeof(uint16_t) * globals.capacity, globals.host_mem);
                    exit(1);
                }
                globals.multiplicity16 = new_multi16;
            }
        }

        ++num_edges;
        if (globals.mult_mem_type == 1) {
            edge_p[globals.words_per_edge - 1] &= 0xFFFFFF00U;
            edge_p[globals.words_per_edge - 1] |= (edge_p[globals.words_per_edge] >> 8) & 0xFFU;
            globals.multiplicity8[num_edges - 1] = edge_p[globals.words_per_edge] & 0xFFU;
        } else if (globals.mult_mem_type == 2) {
            globals.multiplicity16[num_edges - 1] = edge_p[globals.words_per_edge] & kMaxMulti_t;
        }

        edge_p += globals.words_per_edge;
    }
    edge_reader.destroy();

    // --- realloc ---
    globals.packed_edges = (edge_word_t *) ReAllocAndCheck(globals.packed_edges, sizeof(edge_word_t) * globals.words_per_edge * num_edges, __FILE__, __LINE__);
    if (globals.mult_mem_type == 1) {
        globals.multiplicity8 = (uint8_t*) ReAllocAndCheck(globals.multiplicity8, sizeof(uint8_t) * num_edges, __FILE__, __LINE__);
    } else if (globals.mult_mem_type == 2) {
        globals.multiplicity16 = (uint16_t*) ReAllocAndCheck(globals.multiplicity16, sizeof(uint16_t) * num_edges, __FILE__, __LINE__);
    }


    if (sdbg_builder_verbose >= 2) {
        log("[B::%s] Number of mercy edges: %lld\n", __func__, num_edges - globals.num_edges);
    }
    globals.num_edges = num_edges;
    globals.mem_packed_edges = bytes_per_edge * num_edges;
    return num_edges;
}

/**
 * @brief build lkt for faster binary search for mercy
 */
void InitLookupTable(int64_t *lookup_table, uint32_t *packed_edges, int64_t num_edges, int words_per_edge) {
    memset(lookup_table, 0xFF, sizeof(int64_t) * kLookUpSize * 2);
    uint32_t *edge_p = packed_edges;
    uint32_t cur_prefix = *packed_edges >> kLookUpShift;
    lookup_table[cur_prefix * 2] = 0;
    edge_p += words_per_edge;
    for (int64_t i = 1; i < num_edges; ++i) {
        if ((*edge_p >> kLookUpShift) > cur_prefix) {
            lookup_table[cur_prefix * 2 + 1] = i - 1;
            cur_prefix = (*edge_p >> kLookUpShift);
            lookup_table[cur_prefix * 2] = i;
        } else {
            assert(cur_prefix == (*edge_p >> kLookUpShift));   
        }
        edge_p += words_per_edge;
    }
    lookup_table[cur_prefix * 2 + 1] = num_edges - 1;
}

/**
 * @brief init memory for bucket scan
 */
void PrepareBucketScan(struct global_data_t &globals) {
    // init read partitions
    for (int t = 0; t < globals.num_cpu_threads; ++t) {
        struct readpartition_data_t &rp = globals.readpartitions[t];
        rp.rp_id = t;
        rp.globals = &globals;
        rp.rp_bucket_sizes = (int64_t *) MallocAndCheck(phase2::kNumBuckets * sizeof(int64_t), __FILE__, __LINE__);
        assert(rp.rp_bucket_sizes != NULL);
        rp.rp_bucket_offsets = (int64_t *) MallocAndCheck(phase2::kNumBuckets * sizeof(int64_t), __FILE__, __LINE__);
        assert(rp.rp_bucket_offsets != NULL);
        // distribute reads to partitions
        int64_t average = globals.num_edges / globals.num_cpu_threads;
        rp.rp_start_id = t * average;
        rp.rp_end_id = t < globals.num_cpu_threads - 1 ? (t + 1) * average : globals.num_edges;
        rp.rp_lv1_differential_base = EncodeEdgeOffset(rp.rp_start_id, 0, 0, globals.k_num_bits);
    }
    // init bucket partitions
    for (int t = 0; t < globals.num_cpu_threads - globals.phase2_num_output_threads; ++t) {
        struct bucketpartition_data_t &bp = globals.bucketpartitions[t];
        bp.bp_id = t;
        bp.globals = &globals;
    }
    globals.bucket_sizes = (int64_t *) MallocAndCheck(phase2::kNumBuckets * sizeof(int64_t), __FILE__, __LINE__);
}

void* PreprocessScanToFillBucketSizesThread(void *_data) {
    struct readpartition_data_t &rp = *((struct readpartition_data_t*) _data);
    struct global_data_t &globals = *(rp.globals);
    int64_t *bucket_sizes = rp.rp_bucket_sizes;
    memset(bucket_sizes, 0, phase2::kNumBuckets * sizeof(int64_t));
    edge_word_t *edge_p = globals.packed_edges + rp.rp_start_id * globals.words_per_edge;
    for (int64_t read_id = rp.rp_start_id; read_id < rp.rp_end_id; ++read_id, edge_p += globals.words_per_edge) {
        edge_word_t key = 0; // $$$$$$$$
        edge_word_t *word_p = edge_p;
        edge_word_t w = *(word_p++);
        // build initial partial key
        for (int i = 0; i < phase2::kBucketPrefixLength - 1; ++i) {
            key = key * phase2::kBucketBase + (w >> kTopCharShift) + 1;
            w <<= kBitsPerEdgeChar;
        }
        // 3 edges Sb$, aSb, $aS
        for (int i = phase2::kBucketPrefixLength - 1; i <= phase2::kBucketPrefixLength + 1; ++i) {
            if (i % kCharsPerEdgeWord == 0) {
                w = *(word_p++);
            }
            key = (key * phase2::kBucketBase + (w >> kTopCharShift) + 1) % phase2::kNumBuckets;
            w <<= kBitsPerEdgeChar;
            bucket_sizes[key]++;
        }

        // reverse complement very sucking
        key = 0;
        for (int i = 0; i < phase2::kBucketPrefixLength - 1; ++i) {
            key = key * phase2::kBucketBase + (3 - ExtractNthChar(edge_p, globals.kmer_k - i)) + 1; // complement
        }
        for (int i = phase2::kBucketPrefixLength - 1; i <= phase2::kBucketPrefixLength + 1; ++i) {
            key = key * phase2::kBucketBase + (3 - ExtractNthChar(edge_p, globals.kmer_k - i)) + 1;
            key %= phase2::kNumBuckets;
            bucket_sizes[key]++;
        }
    }
    return NULL;
}

/**
 * @brief Fill bucket sizes
 */
void PreprocessScanToFillBucketSizes(struct global_data_t &globals) {
    // create threads
    for (int t = 1; t < globals.num_cpu_threads; ++t) {
        pthread_create(&(globals.readpartitions[t].thread), NULL, PreprocessScanToFillBucketSizesThread, &globals.readpartitions[t]);
    }
    PreprocessScanToFillBucketSizesThread(&globals.readpartitions[0]);
    for (int t = 1; t < globals.num_cpu_threads; ++t) {
        pthread_join(globals.readpartitions[t].thread, NULL);
    }
    // sum up readpartitions bucketsizes to form global bucketsizes
    int64_t *bucket_sizes = globals.bucket_sizes;
    memset(bucket_sizes, 0, phase2::kNumBuckets * sizeof(int64_t));
    // the array accesses in this loop are optimized by the compiler??
    for (int t = 0; t < globals.num_cpu_threads; ++t) {
        for (int b = 0; b < phase2::kNumBuckets; ++b) {
            bucket_sizes[b] += globals.readpartitions[t].rp_bucket_sizes[b];
        }
    }
}

void InitGlobalData(global_data_t &globals) {
    // --- init output ---
    globals.sdbg_writer.init((std::string(globals.output_prefix)+".w").c_str(),
        (std::string(globals.output_prefix)+".last").c_str(),
        (std::string(globals.output_prefix)+".isd").c_str());
    globals.dummy_nodes_writer.init((std::string(globals.output_prefix)+".dn").c_str());
    globals.output_f_file = OpenFileAndCheck((std::string(globals.output_prefix)+".f").c_str(), "w");
    globals.output_multiplicity_file = OpenFileAndCheck((std::string(globals.output_prefix)+".mul").c_str(), "wb");
    globals.output_multiplicity_file2 = OpenFileAndCheck((std::string(globals.output_prefix)+".mul2").c_str(), "wb");

    // --- compute k_num_bits ---
    {
        globals.k_num_bits = 0;
        int len = 1;
        while (len < globals.kmer_k + 1) {
            globals.k_num_bits++;
            len *= 2;
        }
    }

    // --- init stat ---
    globals.cur_prefix = -1;
    globals.cur_suffix_first_char = -1;
    globals.num_ones_in_last = 0;
    globals.total_number_edges = 0;
    globals.num_dollar_nodes = 0;
    memset(globals.num_chars_in_w, 0, sizeof(globals.num_chars_in_w));

    // --- fill bucket size ---
    xtimer_t timer;
    timer.reset();
    timer.start();

    if (sdbg_builder_verbose >= 3) {
        log("[B::%s] Filling edge partition buckets...\n", __func__);
    }
    pthread_mutex_init(&globals.lv1_items_scanning_lock, NULL); // init lock
    PrepareBucketScan(globals);
    PreprocessScanToFillBucketSizes(globals); // Multithread: fill the read partition buckets, then sum up into the global buckets
    globals.max_bucket_size = *std::max_element(globals.bucket_sizes, globals.bucket_sizes + phase2::kNumBuckets);
    globals.tot_bucket_size = 0;
    for (int i = 0; i < phase2::kNumBuckets; ++i) { globals.tot_bucket_size += globals.bucket_sizes[i]; }
    timer.stop();

    if (sdbg_builder_verbose >= 3) {
        log("[B::%s] Done. Time elapsed: %.4lfs\n", __func__, timer.elapsed());
    }

    // --- calculate lv2 memory ---
#ifdef DISABLE_GPU
    globals.max_lv2_items = std::max(globals.max_bucket_size, kMinLv2BatchSize);
#else
    int64_t lv2_mem = globals.gpu_mem - 1073741824; // should reserver ~1G for GPU sorting
    globals.max_lv2_items = lv2_mem / GPU_BYTES_PER_ITEM;
    globals.max_lv2_items = std::min(globals.max_lv2_items, globals.tot_bucket_size);
    if (globals.max_bucket_size > globals.max_lv2_items) {
        err("[ERROR B::%s] Bucket too large for GPU: contains %lld items. Please try CPU version.\n", __func__, globals.max_bucket_size);
        // TODO: auto switch to CPU version
        exit(1);
    }
#endif
    globals.words_per_substring = DivCeiling(globals.kmer_k * kBitsPerEdgeChar + kBWTCharNumBits + 1 + kBitsPerMulti_t, kBitsPerEdgeWord);
    globals.words_per_dummy_node = DivCeiling(globals.kmer_k * kBitsPerEdgeChar, kBitsPerEdgeWord);
    // lv2 bytes: substring (double buffer), permutation, aux
    int64_t lv2_bytes_per_item = (globals.words_per_substring * sizeof(edge_word_t) + sizeof(uint32_t)) * 2 + sizeof(unsigned char);
#ifdef DISABLE_GPU
    lv2_bytes_per_item += sizeof(uint64_t) * 2; // simulate GPU
#endif

    if (sdbg_builder_verbose >= 2) {
        log("[B::%s] %d words per substring, k_num_bits: %d, words per dummy node ($v): %d\n", __func__, globals.words_per_substring, globals.k_num_bits, globals.words_per_dummy_node);
    }

    // --- memory stuff ---
    int64_t mem_remained = globals.host_mem
                         - globals.mem_packed_edges
                         - phase2::kNumBuckets * sizeof(int64_t) * (globals.num_cpu_threads * 3 + 1);
    if (globals.mem_flag == 1) {
        // auto set memory
        globals.max_lv1_items = std::max(globals.max_lv2_items, int64_t(globals.tot_bucket_size / (kDefaultLv1ScanTime - 0.5)));
        int64_t mem_needed = globals.max_lv1_items * LV1_BYTES_PER_ITEM + globals.max_lv2_items * lv2_bytes_per_item;
        if (mem_needed > mem_remained) {
            GetNumItems(globals, mem_remained, lv2_bytes_per_item);
        }
    } else if (globals.mem_flag == 0) {
        // min memory
        globals.max_lv1_items = std::max(globals.max_lv2_items, int64_t(globals.tot_bucket_size / (kMaxLv1ScanTime - 0.5)));
        int64_t mem_needed = globals.max_lv1_items * LV1_BYTES_PER_ITEM + globals.max_lv2_items * lv2_bytes_per_item;
        if (mem_needed > mem_remained) {
            GetNumItems(globals, mem_remained, lv2_bytes_per_item);
        } else {
            GetNumItems(globals, mem_needed, lv2_bytes_per_item);
        }
    } else {
        // use all
        GetNumItems(globals, mem_remained, lv2_bytes_per_item);
    }

    if (sdbg_builder_verbose >= 2) {
        log("[B::%s] Memory for edges: %lld\n", __func__, globals.mem_packed_edges);
        log("[B::%s] max # lv.1 items = %lld\n", __func__, globals.max_lv1_items);
        log("[B::%s] max # lv.2 items = %lld\n", __func__, globals.max_lv2_items);
    }

    // --- alloc memory ---
    globals.lv1_items = (int*) MallocAndCheck(globals.max_lv1_items * sizeof(int), __FILE__, __LINE__);
    globals.lv2_substrings = (edge_word_t*) MallocAndCheck(globals.max_lv2_items * globals.words_per_substring * sizeof(edge_word_t), __FILE__, __LINE__);
    globals.permutation = (uint32_t *) MallocAndCheck(globals.max_lv2_items * sizeof(uint32_t), __FILE__, __LINE__);
    globals.lv2_substrings_to_output = (edge_word_t*) MallocAndCheck(globals.max_lv2_items * globals.words_per_substring * sizeof(edge_word_t), __FILE__, __LINE__);
    globals.permutation_to_output = (uint32_t *) MallocAndCheck(globals.max_lv2_items * sizeof(uint32_t), __FILE__, __LINE__);
    globals.lv2_aux = (unsigned char*) MallocAndCheck(globals.max_lv2_items * sizeof(unsigned char), __FILE__, __LINE__);
#ifdef DISABLE_GPU
    globals.cpu_sort_space = (uint64_t*) MallocAndCheck(sizeof(uint64_t) * globals.max_lv2_items, __FILE__, __LINE__); // simulate GPU
#endif

    // --- write header ---
    fprintf(globals.output_f_file, "-1\n");
    globals.dummy_nodes_writer.output(globals.words_per_dummy_node);
}

struct ReadReadsThreadData {
    ReadPackage *read_package;
    gzFile *read_file;
};

static void* ReadReadsThread(void* data) {
    ReadPackage &package = *(((ReadReadsThreadData*)data)->read_package);
    gzFile &read_file = *(((ReadReadsThreadData*)data)->read_file);
    package.ReadBinaryReads(read_file);
    return NULL;
}

/**
 * @brief search mercy kmer
 */
inline int64_t BinarySearchKmer(uint32_t *packed_edges, int64_t *lookup_table, int words_per_edge, 
    int words_per_kmer, int last_shift, uint32_t *kmer) {
    // --- first look up ---
    int64_t l = lookup_table[(*kmer >> kLookUpShift) * 2];
    if (l == -1) { return -1; }
    int64_t r = lookup_table[(*kmer >> kLookUpShift) * 2 + 1];
    int64_t mid;

    // --- search the words before the last word ---
    for (int i = 0; i < words_per_kmer - 1; ++i) {
        while (l <= r) {
            mid = (l + r) / 2;
            if (packed_edges[mid * words_per_edge + i] < kmer[i]) {
                l = mid + 1;
            } else if (packed_edges[mid * words_per_edge + i] > kmer[i]) {
                r = mid - 1;
            } else {
                int64_t ll = l, rr = mid, mm;
                while (ll < rr) {
                    mm = (ll + rr) / 2;
                    if (packed_edges[mm * words_per_edge + i] < kmer[i]) {
                        ll = mm + 1;
                    } else {
                        rr = mm;
                    }
                }
                l = ll;

                ll = mid, rr = r;
                while (ll < rr) {
                    mm = (ll + rr + 1) / 2;
                    if (packed_edges[mm * words_per_edge + i] > kmer[i]) {
                        rr = mm - 1;
                    } else {
                        ll = mm;
                    }
                }
                r = rr;
                break;
            }
        }
        if (l > r) { return -1; }
    }

    // --- search the last word ---
    while (l <= r) {
        mid = (l + r) / 2;
        if ((packed_edges[mid * words_per_edge + words_per_kmer - 1] >> last_shift) ==
            (kmer[words_per_kmer - 1] >> last_shift)) {
            return mid;
        } else if ((packed_edges[mid * words_per_edge + words_per_kmer - 1] >> last_shift) > 
            (kmer[words_per_kmer - 1] >> last_shift)) {
            r = mid - 1;
        } else {
            l = mid + 1;
        }
    }
    return -1;
}

//TODO: many hard-code in this function, feel tired to love...
/**
 * @brief read candidate reads and search mercy kmer
 */
void ReadReadsAndGetMercyEdges(global_data_t &globals) {
    assert((globals.edge_lookup = (int64_t *) MallocAndCheck(kLookUpSize * 2 * sizeof(int64_t), __FILE__, __LINE__)) != NULL);
    InitLookupTable(globals.edge_lookup, globals.packed_edges, globals.num_edges, globals.words_per_edge);

    uint32_t *packed_edges = globals.packed_edges;
    int64_t *lookup_table = globals.edge_lookup;
    int kmer_k = globals.kmer_k;
    int words_per_edge = globals.words_per_edge;
    const char *edge_file_prefix = globals.phase2_input_prefix;

    gzFile candidate_file = gzopen((std::string(globals.phase2_input_prefix)+".cand").c_str(), "r");
    ReadPackage read_package[2];
    read_package[0].init(globals.max_read_length);
    read_package[1].init(globals.max_read_length);

    ReadReadsThreadData input_thread_data;
    input_thread_data.read_package = &read_package[0];
    input_thread_data.read_file = &candidate_file;
    int input_thread_idx = 0;
    pthread_t input_thread;

    pthread_create(&input_thread, NULL, ReadReadsThread, &input_thread_data);
    int num_threads = globals.num_cpu_threads - 1;
    omp_set_num_threads(num_threads);

    std::vector<FILE*> out_files;
    for (int i = 0; i < num_threads; ++i) {
        char file_name[10240];
        sprintf(file_name, "%s.mercy.%d", edge_file_prefix, i);
        out_files.push_back(OpenFileAndCheck(file_name, "wb"));
        assert(out_files.back() != NULL);
    }

    // parameters for binary search
    int words_per_kmer = DivCeiling(kmer_k * kBitsPerEdgeChar, kBitsPerEdgeWord);
    int last_shift_k = (kmer_k * kBitsPerEdgeChar) % kBitsPerEdgeWord;
    if (last_shift_k > 0) {
        last_shift_k = kBitsPerEdgeWord - last_shift_k;
    }
    int words_per_k_plus_one = DivCeiling((kmer_k + 1) * kBitsPerEdgeChar, kBitsPerEdgeWord);
    int last_shift_k_plus_one = ((kmer_k + 1) * kBitsPerEdgeChar) % kBitsPerEdgeWord;
    if (last_shift_k_plus_one > 0) {
        last_shift_k_plus_one = kBitsPerEdgeWord - last_shift_k_plus_one;
    }
    // log("%d %d %d %d\n", words_per_kmer, words_per_k_plus_one, last_shift_k, last_shift_k_plus_one);
    uint32_t *kmers = (uint32_t *) MallocAndCheck(sizeof(uint32_t) * (omp_get_max_threads()) * words_per_edge, __FILE__, __LINE__);
    uint32_t *rev_kmers = (uint32_t *) MallocAndCheck(sizeof(uint32_t) * (omp_get_max_threads()) * words_per_edge, __FILE__, __LINE__);
    bool *has_ins = (bool*) MallocAndCheck(sizeof(uint32_t) * (omp_get_max_threads()) * read_package[0].max_read_len, __FILE__, __LINE__);
    bool *has_outs = (bool*) MallocAndCheck(sizeof(uint32_t) * (omp_get_max_threads()) * read_package[0].max_read_len, __FILE__, __LINE__);
    assert(kmers != NULL);
    assert(rev_kmers != NULL);
    assert(has_ins != NULL);
    assert(has_outs != NULL);

    int64_t num_mercy_edges = 0;
    int64_t num_reads = 0;

    while (true) {
        pthread_join(input_thread, NULL);
        ReadPackage &package = read_package[input_thread_idx];
        if (package.num_of_reads == 0) {
            break;
        }

        input_thread_idx ^= 1;
        input_thread_data.read_package = &read_package[input_thread_idx];
        pthread_create(&input_thread, NULL, ReadReadsThread, &input_thread_data);

        num_reads += package.num_of_reads;

#pragma omp parallel for reduction(+:num_mercy_edges)
        for (int read_id = 0; read_id < package.num_of_reads; ++read_id) {
            int read_length = package.length(read_id);
            if (read_length < kmer_k + 2) { continue; }
            bool *has_in = has_ins + omp_get_thread_num() * package.max_read_len;
            bool *has_out = has_outs + omp_get_thread_num() * package.max_read_len;
            memset(has_in, 0, sizeof(bool) * (read_length - kmer_k + 1));
            memset(has_out, 0, sizeof(bool) * (read_length - kmer_k + 1));
            // construct the first kmer
            uint32_t *kmer = kmers + words_per_edge * omp_get_thread_num();
            uint32_t *rev_kmer = rev_kmers + words_per_edge * omp_get_thread_num();
            memcpy(kmer, package.GetReadPtr(read_id), sizeof(uint32_t) * words_per_k_plus_one);
            // construct the rev_kmer
            for (int i = 0; i < words_per_kmer; ++i) {
                rev_kmer[words_per_kmer - 1 - i] = ~ mirror(kmer[i]);
            }
            for (int i = 0; i < words_per_kmer; ++i) {
                rev_kmer[i] <<= last_shift_k;
                rev_kmer[i] |= (i == words_per_kmer - 1) ? 0 : (rev_kmer[i + 1] >> (kBitsPerEdgeWord - last_shift_k));
            }

            int last_index = std::min(read_length - 1, kCharsPerEdgeWord * words_per_k_plus_one - 1);

            // first determine which kmer has in or out
            for (int first_index = 0; first_index + kmer_k <= read_length; ++first_index) {
                if (!has_in[first_index]) {
                    // search the reverse complement
                    if (BinarySearchKmer(packed_edges, lookup_table, words_per_edge, 
                            words_per_kmer, last_shift_k, rev_kmer) != -1) {
                        has_in[first_index] = true;
                    } else {
                        // check whether it has incomings
                        int last_char = kmer[words_per_k_plus_one - 1] & 3;
                        for (int i = words_per_k_plus_one - 1; i > 0; --i) {
                            kmer[i] = (kmer[i] >> 2) | (kmer[i - 1] << 30);
                        }
                        kmer[0] >>= 2;
                        // set the highest char to c
                        for (int c = 0; c < 4; ++c) {
                            kmer[0] &= 0x3FFFFFFF;
                            kmer[0] |= c << 30;
                            if (kmer[0] > rev_kmer[0]) {
                                break;
                            }
                            if (BinarySearchKmer(packed_edges, lookup_table, words_per_edge, 
                                     words_per_k_plus_one, last_shift_k_plus_one, kmer) != -1) {
                                has_in[first_index] = true;
                                break;
                            }
                        }
                        for (int i = 0; i < words_per_k_plus_one - 1; ++i) {
                            kmer[i] = (kmer[i] << 2) | (kmer[i + 1] >> 30);
                        }
                        kmer[words_per_k_plus_one - 1] = (kmer[words_per_k_plus_one - 1] << 2) | last_char;
                    }
                }

                if (true) {
                    // check whether it has outgoing
                    int64_t search_idx = BinarySearchKmer(packed_edges, lookup_table, words_per_edge, 
                                                          words_per_kmer, last_shift_k, kmer);
                    if (search_idx != -1) {
                        has_out[first_index] = true;
                        // a quick check whether next has in
                        if (first_index + kmer_k < read_length && 
                            (packed_edges[search_idx * words_per_edge + words_per_k_plus_one - 1] >> last_shift_k_plus_one) ==
                            (kmer[words_per_k_plus_one - 1] >> last_shift_k_plus_one)) {
                            has_in[first_index + 1] = true;
                        }
                    } else {
                        // search the rc
                        int rc_last_char = rev_kmer[words_per_k_plus_one - 1] & 3;
                        for (int i = words_per_k_plus_one - 1; i > 0; --i) {
                            rev_kmer[i] = (rev_kmer[i] >> 2) | (rev_kmer[i - 1] << 30);
                        }
                        rev_kmer[0] >>= 2;
                        int next_c = first_index + kmer_k < read_length ?
                                     (3 - package.CharAt(read_id, first_index + kmer_k)) :
                                     3;
                        rev_kmer[0] &= 0x3FFFFFFF;
                        rev_kmer[0] |= next_c << 30;
                        if (rev_kmer[0] <= kmer[0] &&
                            BinarySearchKmer(packed_edges, lookup_table, words_per_edge, 
                                words_per_k_plus_one, last_shift_k_plus_one, rev_kmer) != -1) {
                            has_out[first_index] = true;
                            has_in[first_index + 1] = true;
                        }

                        for (int c = 0; !has_out[first_index] && c < 4; ++c) {
                            if (c == next_c) { continue; }
                            rev_kmer[0] &= 0x3FFFFFFF;
                            rev_kmer[0] |= c << 30;
                            if (rev_kmer[0] > kmer[0]) {
                                break;
                            }
                            if (BinarySearchKmer(packed_edges, lookup_table, words_per_edge, 
                                    words_per_k_plus_one, last_shift_k_plus_one, rev_kmer) != -1) {
                                has_out[first_index] = true;
                                break;
                            }
                        }
                        for (int i = 0; i < words_per_k_plus_one - 1; ++i) {
                            rev_kmer[i] = (rev_kmer[i] << 2) | (rev_kmer[i + 1] >> 30);
                        }
                        rev_kmer[words_per_k_plus_one - 1] = (rev_kmer[words_per_k_plus_one - 1] << 2) | rc_last_char;
                    }
                }

                // shift kmer and rev_kmer
                for (int i = 0; i < words_per_k_plus_one - 1; ++i) {
                    kmer[i] = (kmer[i] << 2) | (kmer[i + 1] >> 30);
                }
                kmer[words_per_k_plus_one - 1] <<= 2;
                if (++last_index < read_length) {
                    kmer[words_per_k_plus_one - 1] |= package.CharAt(read_id, last_index);
                }

                for (int i = words_per_k_plus_one - 1; i > 0; --i) {
                    rev_kmer[i] = (rev_kmer[i] >> 2) | (rev_kmer[i - 1] << 30);
                }
                rev_kmer[0] = (rev_kmer[0] >> 2) | ((3 - package.CharAt(read_id, first_index + kmer_k)) << 30);
            }

            // adding mercy edges
            int last_no_out = -1;
            std::vector<bool> is_mercy_edges(read_length - kmer_k, false);
            for (int i = 0; i + kmer_k <= read_length; ++i) {
                switch (has_in[i] | (int(has_out[i]) << 1)) {
                    case 1: { // has incoming only
                        last_no_out = i;
                        break;
                    }
                    case 2: { // has outgoing only
                        if (last_no_out >= 0) {
                            for (int j = last_no_out; j < i; ++j) {
                                is_mercy_edges[j] = true;
                            }
                            num_mercy_edges += i - last_no_out;
                        }
                        last_no_out = -1;
                        break;
                    }
                    case 3: { // has in and out
                        last_no_out = -1;
                        break;
                    }
                    default: {
                        // do nothing
                        break;
                    }
                }
            }
            
            memcpy(kmer, package.GetReadPtr(read_id), sizeof(uint32_t) * words_per_k_plus_one);
            last_index = std::min(read_length - 1, kCharsPerEdgeWord * words_per_k_plus_one - 1);
            for (int i = 0; i + kmer_k < read_length; ++i) {
                if (is_mercy_edges[i]) {
                    uint32_t last_word = kmer[words_per_k_plus_one - 1];
                    kmer[words_per_k_plus_one - 1] >>= last_shift_k_plus_one;
                    kmer[words_per_k_plus_one - 1] <<= last_shift_k_plus_one;
                    for (int j = words_per_k_plus_one; j < words_per_edge; ++j) {
                        kmer[j] = 0;
                    }
                    if (globals.mult_mem_type == 0) {
                        kmer[words_per_edge - 1] |= 1; // WARNING: only accurate when m=2, but I think doesn't matter a lot
                    }
                    fwrite(kmer, sizeof(uint32_t), words_per_edge, out_files[omp_get_thread_num()]);
                    if (globals.mult_mem_type > 0) {
                        uint32_t kMercyMult = 1;
                        fwrite(&kMercyMult, sizeof(uint32_t), 1, out_files[omp_get_thread_num()]);
                    }
                    kmer[words_per_k_plus_one - 1] = last_word;
                }

                for (int i = 0; i < words_per_k_plus_one - 1; ++i) {
                    kmer[i] = (kmer[i] << 2) | (kmer[i + 1] >> 30);
                }
                kmer[words_per_k_plus_one - 1] <<= 2;
                if (++last_index < read_length) {
                    kmer[words_per_k_plus_one - 1] |= package.CharAt(read_id, last_index);
                }
            }
        }
        if (num_reads % (16 * package.kMaxNumReads) == 0) {
            if (sdbg_builder_verbose >= 4) {
                log("[B::%s] Number of reads: %ld, Number of mercy edges: %ld\n", __func__, num_reads, num_mercy_edges);
            }
        }
    }


    if (sdbg_builder_verbose >= 2) {
        log("[B::%s] Number of reads: %ld, Number of mercy edges: %ld\n", __func__, num_reads, num_mercy_edges);
    }

    free(kmers);
    free(rev_kmers);
    free(has_ins);
    free(has_outs);
    free(globals.edge_lookup);
    for (unsigned i = 0; i < out_files.size(); ++i) {
        fclose(out_files[i]);
    }
}

/**
 * @brief worker thread for Lv1ScanToFillOffests
 */
void* Lv1ScanToFillOffsetsThread(void *_data) {
    struct readpartition_data_t &rp = *((struct readpartition_data_t*) _data);
    struct global_data_t &globals = *(rp.globals);
    int64_t *prev_full_offsets = (int64_t *)MallocAndCheck(phase2::kNumBuckets * sizeof(int64_t), __FILE__, __LINE__); // temporary array for computing differentials
    assert(prev_full_offsets != NULL);
    for (int b = globals.lv1_start_bucket; b < globals.lv1_end_bucket; ++b)
        prev_full_offsets[b] = rp.rp_lv1_differential_base;
    // this loop is VERY similar to that in PreprocessScanToFillBucketSizesThread
    edge_word_t *edge_p = globals.packed_edges + rp.rp_start_id * globals.words_per_edge;

    // ===== this is a macro to save some copy&paste ================
#define CHECK_AND_SAVE_OFFSET_PHASE2(offset, strand)                                   \
    do {                                                                \
      if (((key - globals.lv1_start_bucket) ^ (key - globals.lv1_end_bucket)) & kSignBitMask) { \
        int64_t full_offset = EncodeEdgeOffset(read_id, offset, strand, globals.k_num_bits); \
        int64_t differential = full_offset - prev_full_offsets[key];      \
        if (differential > kDifferentialLimit) {                      \
          pthread_mutex_lock(&globals.lv1_items_scanning_lock); \
          globals.lv1_items[ rp.rp_bucket_offsets[key]++ ] = -globals.lv1_items_special.size() - 1; \
          globals.lv1_items_special.push_back(full_offset);                  \
          pthread_mutex_unlock(&globals.lv1_items_scanning_lock); \
        } else {                                                              \
          assert(differential >= 0); \
          globals.lv1_items[ rp.rp_bucket_offsets[key]++ ] = (int) differential; \
        } \
        prev_full_offsets[key] = full_offset;                           \
      }                                                                 \
    } while (0)
    // ^^^^^ why is the macro surrounded by a do-while? please ask Google
    // =========== end macro ==========================

    for (int64_t read_id = rp.rp_start_id; read_id < rp.rp_end_id; ++read_id, edge_p += globals.words_per_edge) {
        edge_word_t key = 0; // $$$$$$$$
        edge_word_t *word_p = edge_p;
        edge_word_t w = *(word_p++);
        // build initial partial key
        for (int i = 0; i < phase2::kBucketPrefixLength - 1; ++i) {
            key = key * phase2::kBucketBase + (w >> kTopCharShift) + 1;
            w <<= kBitsPerEdgeChar;
        }
        // 3 edges Sb$, aSb, $aS
        for (int i = kBucketPrefixLength - 1; i <= kBucketPrefixLength + 1; ++i) {
            if (i % kCharsPerEdgeWord == 0) {
                w = *(word_p++);
            }
            key = (key * phase2::kBucketBase + (w >> kTopCharShift) + 1) % phase2::kNumBuckets;
            w <<= kBitsPerEdgeChar;
            CHECK_AND_SAVE_OFFSET_PHASE2(i - kBucketPrefixLength + 1, 0);
        }

        // reverse complement very sucking
        key = 0;
        for (int i = 0; i < phase2::kBucketPrefixLength - 1; ++i) {
            key = key * phase2::kBucketBase + (3 - ExtractNthChar(edge_p, globals.kmer_k - i)) + 1; // complement
        }
        for (int i = phase2::kBucketPrefixLength - 1; i <= phase2::kBucketPrefixLength + 1; ++i) {
            key = (key * phase2::kBucketBase) + (3 - ExtractNthChar(edge_p, globals.kmer_k - i)) + 1;
            key %= phase2::kNumBuckets;
            CHECK_AND_SAVE_OFFSET_PHASE2(i - kBucketPrefixLength + 1, 1);
        }
    }

    free(prev_full_offsets);
    return NULL;
}

/**
 * @brief Fill Lv1 items to host mem
 */
void Lv1ScanToFillOffests(struct global_data_t &globals) {
    globals.lv1_items_special.clear();
    Lv1ComputeBucketOffset(globals);
    // create threads
    for (int t = 1; t < globals.num_cpu_threads; ++t) {
        pthread_create(&(globals.readpartitions[t].thread), NULL, Lv1ScanToFillOffsetsThread, &globals.readpartitions[t]);
    }
    Lv1ScanToFillOffsetsThread(&globals.readpartitions[0]);
    for (int t = 1; t < globals.num_cpu_threads; ++t) {
        pthread_join(globals.readpartitions[t].thread, NULL);
    }
    // revert rp_bucket_offsets
    Lv1ComputeBucketOffset(globals);
}

void CopySubstring(edge_word_t* dest, edge_word_t* src_read, int offset, int num_chars_to_copy, int counting, global_data_t &globals) {
    int64_t spacing = globals.lv2_num_items;
    int words_per_edge = globals.words_per_edge;
    int64_t words_per_substring = globals.words_per_substring;
    int kmer_k = globals.kmer_k;

    // copy words of the suffix to the suffix pool
    int which_word = offset / kCharsPerEdgeWord;
    int word_offset = offset % kCharsPerEdgeWord;
    edge_word_t *src_p = src_read + which_word;
    edge_word_t *dest_p = dest;
    int num_words_copied = 0;
    if (!word_offset) { // special case (word aligned), easy
        while (which_word < words_per_edge && num_words_copied < words_per_substring) {
            *dest_p = *src_p; // write out
            dest_p += spacing;
            src_p++;
            which_word++;
            num_words_copied++;
        }
    } else { // not word-aligned
        int bit_shift = offset * kBitsPerEdgeChar;
        edge_word_t s = *src_p;
        edge_word_t d = s << bit_shift;
        which_word++;
        while (which_word < words_per_edge) {
            s = *(++src_p);
            d |= s >> (kBitsPerEdgeWord - bit_shift);
            *dest_p = d; // write out
            if (++num_words_copied >= words_per_substring) goto here2;
            dest_p += spacing;
            d = s << bit_shift;
            which_word++;
        }
        *dest_p = d; // write last word
here2:
        ;
    }

    {
        // now mask the extra bits (TODO can be optimized)
        int num_bits_to_copy = num_chars_to_copy * 2;
        int which_word = num_bits_to_copy / kBitsPerEdgeWord;
        edge_word_t *p = dest + which_word * spacing;
        int bits_to_clear = kBitsPerEdgeWord - num_bits_to_copy % kBitsPerEdgeWord;
        if (bits_to_clear < kBitsPerEdgeWord) {
            *p >>= bits_to_clear;
            *p <<= bits_to_clear;
        } else if (which_word < globals.words_per_substring) {
            *p = 0;
        }
        ++which_word;
        while (which_word < words_per_substring) { // fill zero
            *(p+=spacing) = 0;
            which_word++;
        }
    }

    int prev_char;
    if (offset == 0) {
        assert(num_chars_to_copy == globals.kmer_k);
        prev_char = kSentinelValue;
    } else {
        prev_char = ExtractNthChar(src_read, offset - 1);
    }

    edge_word_t *last_word = dest + (words_per_substring - 1) * spacing;
    *last_word |= int(num_chars_to_copy == kmer_k) << (kBWTCharNumBits + kBitsPerMulti_t);
    *last_word |= prev_char << kBitsPerMulti_t;
    *last_word |= std::min(counting, kMaxMulti_t);
}

void CopySubstringRC(edge_word_t* dest, edge_word_t* src_read, int offset, int num_chars_to_copy, int counting, global_data_t &globals) {
    int64_t spacing = globals.lv2_num_items;
    int which_word = (globals.kmer_k - offset) / kCharsPerEdgeWord;
    int word_offset = (globals.kmer_k - offset) % kCharsPerEdgeWord;
    edge_word_t *dest_p = dest;

    if (word_offset == kCharsPerEdgeWord - 1) { // edge_word_t aligned
        for (int i = 0; i < globals.words_per_substring && i <= which_word; ++i) {
            *dest_p = ~ mirror(src_read[which_word - i]);
            dest_p += spacing;
        }
    } else {
        int bit_offset = (kCharsPerEdgeWord - 1 - word_offset) * kBitsPerEdgeChar;
        int i;
        edge_word_t w;
        for (i = 0; i < globals.words_per_substring - 1 && i < which_word; ++i) {
            w = (src_read[which_word - i] >> bit_offset) |
                                      (src_read[which_word - i - 1] << (kBitsPerEdgeWord - bit_offset));
            *dest_p = ~ mirror(w);
            dest_p += spacing;
        }
        // last word
        w = src_read[which_word - i] >> bit_offset;
        if (which_word >= i + 1) {
            w |= (src_read[which_word - i - 1] << (kBitsPerEdgeWord - bit_offset));
        }
        *dest_p = ~ mirror(w);
    }

    {
        // now mask the extra bits (TODO can be optimized)
        int num_bits_to_copy = num_chars_to_copy * 2;
        int which_word = num_bits_to_copy / kBitsPerEdgeWord;
        edge_word_t *p = dest + which_word * spacing;
        int bits_to_clear = kBitsPerEdgeWord - num_bits_to_copy % kBitsPerEdgeWord;
        if (bits_to_clear < kBitsPerEdgeWord) {
            *p >>= bits_to_clear;
            *p <<= bits_to_clear;
        } else if (which_word < globals.words_per_substring) {
            *p = 0;
        }
        ++which_word;
        while (which_word < globals.words_per_substring) { // fill zero
            *(p+=spacing) = 0;
            which_word++;
        }
    }

    int prev_char;
    if (offset == 0) {
        assert(num_chars_to_copy == globals.kmer_k);
        prev_char = kSentinelValue;
    } else {
        prev_char = 3 - ExtractNthChar(src_read, globals.kmer_k - (offset - 1));
    }

    edge_word_t *last_word = dest + (globals.words_per_substring - 1) * spacing;
    *last_word |= int(num_chars_to_copy == globals.kmer_k) << (kBWTCharNumBits + kBitsPerMulti_t);
    *last_word |= prev_char << kBitsPerMulti_t;
    *last_word |= std::min(counting, kMaxMulti_t);
}

inline int BucketToPrefix(int x) {
    int y = 0;
    for (int i=0; i < phase2::kBucketPrefixLength; ++i) {
        int z = x % phase2::kBucketBase;
        if (z > 0) { --z; }
        y |= (z << (i * kBitsPerEdgeChar));
        x /= phase2::kBucketBase;
    }
    return y;
}

/**
 * @brief worker thread for Lv2ExtractSubstrings
 */
void* Lv2ExtractSubstringsThread(void* _data) {
    struct bucketpartition_data_t &bp = *((struct bucketpartition_data_t*) _data);
    struct global_data_t &globals = *(bp.globals);
    int *lv1_p = globals.lv1_items + globals.readpartitions[0].rp_bucket_offsets[ bp.bp_start_bucket ];
    int64_t offset_mask = (1 << globals.k_num_bits) - 1; // 0000....00011..11
    edge_word_t *substrings_p = globals.lv2_substrings +
                         (globals.readpartitions[0].rp_bucket_offsets[ bp.bp_start_bucket ] - globals.readpartitions[0].rp_bucket_offsets[ globals.lv2_start_bucket ]);
    for (int bucket = bp.bp_start_bucket; bucket < bp.bp_end_bucket; ++bucket) {
        for (int t = 0; t < globals.num_cpu_threads; ++t) {
            int64_t full_offset = globals.readpartitions[t].rp_lv1_differential_base;
            int num = globals.readpartitions[t].rp_bucket_sizes[bucket];
            for (int i = 0; i < num; ++i) {
                if (*lv1_p >= 0) {
                    full_offset += *(lv1_p++);
                } else {
                    full_offset = globals.lv1_items_special[-1 - *(lv1_p++)];
                }
                int64_t read_id = full_offset >> (1 + globals.k_num_bits);
                int offset = full_offset & offset_mask;
                int strand = (full_offset >> globals.k_num_bits) & 1;
                edge_word_t *edge_p = globals.packed_edges + read_id * globals.words_per_edge;
                int num_chars_to_copy = globals.kmer_k - (offset >= 2);
                int counting = 0;
                if (offset == 1) {
                    switch (globals.mult_mem_type) {
                      case 0:
                        counting = *(edge_p + globals.words_per_edge - 1) & kMaxMulti_t;
                        break;
                      case 1:
                        counting = ((*(edge_p + globals.words_per_edge - 1) & 0xFF) << 8) | globals.multiplicity8[read_id];
                        break;
                      case 2:
                        counting = globals.multiplicity16[read_id];
                        break;
                      default: assert(false);
                    }
                }
                if (strand == 0) {
                    CopySubstring(substrings_p, edge_p, offset, num_chars_to_copy, counting, globals);
                } else {
                    CopySubstringRC(substrings_p, edge_p, offset, num_chars_to_copy, counting, globals);
                }
#ifdef DBJ_DEBUG
                if ((*substrings_p >> (32 - phase2::kBucketPrefixLength * 2)) != BucketToPrefix(bucket)) {
                    debug("WRONG substring wrong:%d right:%d(%d) read_id:%lld offset:%d strand: %d num_chars_to_copy:%d\n", *substrings_p >> 16, BucketToPrefix(bucket), bucket, read_id, offset, strand, num_chars_to_copy);
                }
#endif
                substrings_p++;
            }
        }
    }
    return NULL;
}

/**
 * @brief extract true substrings for Lv2 sorting
 */
void Lv2ExtractSubstrings(struct global_data_t &globals) {
    Lv2DistributeBucketPartitions(globals, globals.phase2_num_output_threads);
    // create threads
    for (int t = 0; t < globals.num_cpu_threads-globals.phase2_num_output_threads; ++t) {
        pthread_create(&(globals.bucketpartitions[t].thread), NULL, Lv2ExtractSubstringsThread, &globals.bucketpartitions[t]);
    }
    for (int t = 0; t < globals.num_cpu_threads-globals.phase2_num_output_threads; ++t) {
        pthread_join(globals.bucketpartitions[t].thread, NULL);
    }
#ifdef DBJ_DEBUG
    for (int i = 0; i < globals.lv2_num_items; ++i) {
        edge_word_t *substrings_p = globals.lv2_substrings + i;
        int cur_bucket = *substrings_p >> (32 - phase2::kBucketPrefixLength * 2);
        if (cur_bucket > BucketToPrefix(globals.lv2_end_bucket - 1) || cur_bucket < BucketToPrefix(globals.lv2_start_bucket)) {
            debug("Start: %d, end: %d\n", globals.lv2_start_bucket, globals.lv2_end_bucket);
            debug("%d %d %d\n", cur_bucket, BucketToPrefix(globals.lv2_start_bucket), BucketToPrefix(globals.lv2_end_bucket - 1));
            exit(1);
        }
    }
#endif
}

inline bool IsDiffKMinusOneMer(edge_word_t *item1, edge_word_t *item2, int64_t spacing, int kmer_k) {
    // mask extra bits
    int chars_in_last_word = (kmer_k - 1) % kCharsPerEdgeWord;
    int num_full_words = (kmer_k - 1) / kCharsPerEdgeWord;
    if (chars_in_last_word > 0) {
        edge_word_t w1 = item1[num_full_words * spacing];
        edge_word_t w2 = item2[num_full_words * spacing];
        if ((w1 >> (kCharsPerEdgeWord - chars_in_last_word) * kBitsPerEdgeChar) != (w2 >> (kCharsPerEdgeWord - chars_in_last_word) * kBitsPerEdgeChar)) {
            return true;
        }
    } 

    for (int i = num_full_words - 1; i >= 0; --i) {
        if (item1[i * spacing] != item2[i * spacing]) {
            return true;
        }
    }
    return false;
}

// helper
inline int ExtractFirstChar(edge_word_t *item) {
    return *item >> kTopCharShift;
}

// bS'a
inline int Extract_a(edge_word_t *item, int num_words, int64_t spacing, int kmer_k) {
    int non_dollar = (item[(num_words - 1) * spacing] >> (kBWTCharNumBits + kBitsPerMulti_t)) & 1;
    if (non_dollar) {
        int which_word = (kmer_k - 1) / kCharsPerEdgeWord;
        int word_index = (kmer_k - 1) % kCharsPerEdgeWord;
        return (item[which_word * spacing] >> (kCharsPerEdgeWord - 1 - word_index) * kBitsPerEdgeChar) & kEdgeCharMask;
    } else {
        return kSentinelValue;
    }
}

inline int Extract_b(edge_word_t *item, int num_words, int64_t spacing) {
    return (item[(num_words - 1) * spacing] >> kBitsPerMulti_t) & ((1 << kBWTCharNumBits) - 1);
}

inline int ExtractCounting(edge_word_t *item, int num_words, int64_t spacing) {
    return item[(num_words - 1) * spacing] & kMaxMulti_t; 
}

inline int Extract_a_aux(unsigned char aux) {
    return (aux >> kBWTCharNumBits) & ((1 << kBWTCharNumBits) - 1);
}

inline int Extract_b_aux(unsigned char aux) {
    return aux & ((1 << kBWTCharNumBits) - 1);
}

void *Lv2OutputThread(void *_op) {
    struct outputpartition_data_t *op = (struct outputpartition_data_t*) _op;
    struct global_data_t &globals = *(op->globals);
    int64_t op_start_index = op->op_start_index;
    int64_t op_end_index = op->op_end_index;
    int start_idx, end_idx;
    int has_solid_a = 0; // has solid (k+1)-mer aSb
    int has_solid_b = 0; // has solid aSb
    int last_a[4], outputed_b;

    for (start_idx = op_start_index; start_idx < op_end_index; start_idx = end_idx) {
        end_idx = start_idx + 1;
        edge_word_t *item = globals.lv2_substrings_to_output + globals.permutation_to_output[start_idx];
        while (end_idx < op_end_index && 
               !IsDiffKMinusOneMer(
                    item, 
                    globals.lv2_substrings_to_output + globals.permutation_to_output[end_idx],
                    globals.lv2_num_items_to_output,
                    globals.kmer_k)) {
            ++end_idx;
        }

        // clean marking
        has_solid_a = has_solid_b = 0;
        outputed_b = 0;
        for (int i = start_idx; i < end_idx; ++i) {
            edge_word_t *cur_item = globals.lv2_substrings_to_output + globals.permutation_to_output[i];
            int a = Extract_a(cur_item, globals.words_per_substring, globals.lv2_num_items_to_output, globals.kmer_k);
            int b = Extract_b(cur_item, globals.words_per_substring, globals.lv2_num_items_to_output);

            if (a != kSentinelValue && b != kSentinelValue) {
                has_solid_a |= 1 << a;
                has_solid_b |= 1 << b;
            }
            if (a != kSentinelValue && 
                (b != kSentinelValue || !(has_solid_a & (1 << a)))) {
                last_a[a] = i;
            }
        }

        for (int i = start_idx, j; i < end_idx; i = j) {
            edge_word_t *cur_item = globals.lv2_substrings_to_output + globals.permutation_to_output[i];
            int a = Extract_a(cur_item, globals.words_per_substring, globals.lv2_num_items_to_output, globals.kmer_k);
            int b = Extract_b(cur_item, globals.words_per_substring, globals.lv2_num_items_to_output);

            j = i + 1;
            while (j < end_idx) {
                edge_word_t *next_item = globals.lv2_substrings_to_output + globals.permutation_to_output[j];
                if (Extract_a(next_item, globals.words_per_substring, globals.lv2_num_items_to_output, globals.kmer_k) != a ||
                    Extract_b(next_item, globals.words_per_substring, globals.lv2_num_items_to_output) != b) {
                    break;
                } else {
                    ++j;
                }
            }

            int w, last, is_dollar = 0;

            if (a == kSentinelValue) {
                assert(b != kSentinelValue);
                if (has_solid_b & (1 << b)) {
                    continue;
                }
                is_dollar = 1;
            }

            if (b == kSentinelValue) {
                assert(a != kSentinelValue);
                if (has_solid_a & (1 << a)) {
                    continue;
                }
            }

            w = (b == kSentinelValue) ? 0 : ((outputed_b & (1 << b)) ? b + 5 : b + 1);
            outputed_b |= 1 << b;
            last = (a == kSentinelValue) ? 0 : ((last_a[a] == j - 1) ? 1 : 0);

            assert(!(globals.lv2_aux[i] & (1 << 7)));
            globals.lv2_aux[i] = w | (last << 4) | (is_dollar << 5) | (1 << 7);
        }
    }

    pthread_barrier_wait(&globals.output_barrier);

    if (op_start_index == 0) {
        xtimer_t local_timer;
        local_timer.reset();
        local_timer.start();
        for (int i = 0; i < globals.lv2_num_items_to_output; ++i) {
            if (globals.lv2_aux[i] & (1 << 7)) {
                edge_word_t *item = globals.lv2_substrings_to_output + globals.permutation_to_output[i];
                while (ExtractFirstChar(item) > globals.cur_suffix_first_char) {
                    ++globals.cur_suffix_first_char;
                    fprintf(globals.output_f_file, "%lld\n", (long long)globals.total_number_edges);
                }

                multi_t counting_to_output = std::min(kMaxMulti_t, 
                    ExtractCounting(item, globals.words_per_substring, globals.lv2_num_items_to_output));
                // output
                globals.sdbg_writer.outputW(globals.lv2_aux[i] & 0xF);
                globals.sdbg_writer.outputLast((globals.lv2_aux[i] >> 4) & 1);
                globals.sdbg_writer.outputIsDollar((globals.lv2_aux[i] >> 5) & 1);
                if (counting_to_output <= kMaxMulti2_t) {
                    multi2_t c = counting_to_output;
                    fwrite(&c, sizeof(multi2_t), 1, globals.output_multiplicity_file);   
                } else {
                    int64_t c = counting_to_output | (globals.total_number_edges << 16);
                    fwrite(&c, sizeof(int64_t), 1, globals.output_multiplicity_file2);
                    fwrite(&kMulti2Sp, sizeof(multi2_t), 1, globals.output_multiplicity_file);
                }

                globals.total_number_edges++;
                globals.num_chars_in_w[globals.lv2_aux[i] & 0xF]++;
                globals.num_ones_in_last += (globals.lv2_aux[i] >> 4) & 1;

                if ((globals.lv2_aux[i] >> 5) & 1) {
                    globals.num_dollar_nodes++;
                    if (globals.num_dollar_nodes >= phase2::kMaxDummyEdges) {
                        err("[ERROR B::%s] Too many dummy nodes (>= %lld)! The graph contains too many tips!\n", __func__, (long long)phase2::kMaxDummyEdges);
                        exit(1);
                    }
                    for (int64_t i = 0; i < globals.words_per_dummy_node; ++i) {
                        globals.dummy_nodes_writer.output(item[i * globals.lv2_num_items_to_output]);
                    }
                }
                if ((globals.lv2_aux[i] & 0xF) == 0) {
                    globals.num_dummy_edges++;
                }
            }
        }
        local_timer.stop();

        if (sdbg_builder_verbose >= 4) {
            log("[B::%s] Linear part: %lf\n", __func__, local_timer.elapsed());
        }
    }
    return NULL;
}

void Lv2Output(global_data_t &globals) {
    globals.phase2_output_timer.reset();
    globals.phase2_output_timer.start();
    // distribute threads
    int64_t last_end_index = 0;
    int64_t items_per_thread = globals.lv2_num_items_to_output / globals.phase2_num_output_threads;

    for (int thread_id = 0; thread_id < globals.phase2_num_output_threads - 1; ++thread_id) {
        int64_t this_start_index = last_end_index;
        int64_t this_end_index = this_start_index + items_per_thread;
        if (this_end_index > globals.lv2_num_items_to_output) { this_end_index = globals.lv2_num_items_to_output; }
        if (this_end_index > 0) {
            while (this_end_index < globals.lv2_num_items_to_output) {
                edge_word_t *prev_item = globals.lv2_substrings_to_output + globals.permutation_to_output[this_end_index - 1];
                edge_word_t *item = globals.lv2_substrings_to_output + globals.permutation_to_output[this_end_index];
                if (IsDiffKMinusOneMer(item, prev_item, globals.lv2_num_items_to_output, globals.kmer_k)) {
                    break;
                }
                ++this_end_index;
            }
        }
        globals.outputpartitions[thread_id].op_start_index = this_start_index;
        globals.outputpartitions[thread_id].op_end_index = this_end_index;
        last_end_index = this_end_index;
    }

    // last partition
    globals.outputpartitions[globals.phase2_num_output_threads - 1].op_start_index = last_end_index;
    globals.outputpartitions[globals.phase2_num_output_threads - 1].op_end_index = globals.lv2_num_items_to_output;

    memset(globals.lv2_aux, 0, sizeof(globals.lv2_aux[0]) * globals.lv2_num_items_to_output);
    pthread_barrier_init(&globals.output_barrier, NULL, globals.phase2_num_output_threads);
    for (int thread_id = 0; thread_id < globals.phase2_num_output_threads; ++thread_id) {
        globals.outputpartitions[thread_id].op_id = thread_id;
        globals.outputpartitions[thread_id].globals = &globals;
        pthread_create(&globals.output_threads[thread_id], NULL, Lv2OutputThread, &globals.outputpartitions[thread_id]);
    }
}

void Lv2OutputJoin(global_data_t &globals) {
    for (int thread_id = 0; thread_id < globals.phase2_num_output_threads; ++thread_id) {
        pthread_join(globals.output_threads[thread_id], NULL);
    }

    pthread_barrier_destroy(&globals.output_barrier);
    globals.phase2_output_timer.stop();
}

void Phase2Clean(struct global_data_t &globals) {
    pthread_mutex_destroy(&globals.lv1_items_scanning_lock);
    free(globals.packed_edges);
    free(globals.bucket_sizes);
    free(globals.lv1_items);
    free(globals.lv2_substrings);
    free(globals.permutation);
    free(globals.lv2_substrings_to_output);
    free(globals.permutation_to_output);
    free(globals.lv2_aux);
    fclose(globals.output_f_file);
    fclose(globals.output_multiplicity_file);
    if (globals.mult_mem_type == 1) {
        free(globals.multiplicity8);
    } else if (globals.mult_mem_type == 2) {
        free(globals.multiplicity16);
    }
    globals.dummy_nodes_writer.destroy();
    for (int t = 0; t < globals.num_cpu_threads; ++t) {
        free(globals.readpartitions[t].rp_bucket_sizes);
        free(globals.readpartitions[t].rp_bucket_offsets);
    }
#ifdef DISABLE_GPU
    free(globals.cpu_sort_space);
#endif
}

void Phase2Entry(struct global_data_t &globals) {
    xtimer_t timer;
    // --- read edges ---
    timer.reset();
    timer.start();

    if (sdbg_builder_verbose >= 3) {
        log("[B::%s] Reading edges from temporary files...\n", __func__);
    }
    InitDNAMap();
    ReadEdges(globals);
    timer.stop();

    if (sdbg_builder_verbose >= 3) {
        log("[B::%s] Done. Time elapsed: %.4lfs\n", __func__, timer.elapsed());
    }

    if (globals.need_mercy) {
        timer.reset();
        timer.start();
        if (sdbg_builder_verbose >= 3) {
            log("[B::%s] Adding mercy edges...\n", __func__);
        }
        ReadReadsAndGetMercyEdges(globals);
        ReadMercyEdges(globals);
        timer.stop();

        if (sdbg_builder_verbose >= 3) {
            log("[B::%s] Done. Time elapsed: %.4lfs\n", __func__, timer.elapsed());
        }
    }

    // --- init global data ---
    InitGlobalData(globals);
#ifdef DBJ_DEBUG
    debug("The first 10 edges:\n");
    for (int i = 0; i < 10; ++i) {
        DumpPackedEdge(globals.packed_edges + i * globals.words_per_edge, globals.words_per_edge, globals.kmer_k);
        debug("\n");
    }
    debug("The last 10 edge: ");
    for (int i = 9; i >= 0; --i) {
        DumpPackedEdge(globals.packed_edges + (globals.num_edges - 1 - i) * globals.words_per_edge, globals.words_per_edge, globals.kmer_k);
        debug("\n");
    }
#endif

    ////////////////////////////////// Start processing... ////////////////////////////
    int lv1_iteration = 0;
    //======================================== LEVEL 1 loop ============================================//
    globals.lv1_start_bucket = 0;
    timer.reset();
    timer.start();
    bool output_thread_created = false;
    while (globals.lv1_start_bucket < phase2::kNumBuckets) {
        xtimer_t local_timer;
        lv1_iteration++;
        // --- finds the bucket range for this iteration ---

        if (sdbg_builder_verbose >= 3) {
            log("[B::%s] Finding end bucket...\n", __func__);
        }
        local_timer.reset();
        local_timer.start();
        globals.lv1_end_bucket = FindEndBucket(globals.bucket_sizes, globals.lv1_start_bucket, phase2::kNumBuckets, globals.max_lv1_items, globals.lv1_num_items);
        if (globals.lv1_num_items == 0) { // i.e. can't even hold a single bucket (very unlikely though)
            err("[ERROR] Bucket %d too large for lv.1: contains %lld items\n", globals.lv1_end_bucket, globals.bucket_sizes[globals.lv1_end_bucket]);
            exit(1);
        }

        if (sdbg_builder_verbose >= 3) {
            log("[B::%s] Iteration %d, from bucket %d to %d\n", __func__, lv1_iteration, globals.lv1_start_bucket, globals.lv1_end_bucket-1);
            log("[B::%s] Scanning and filling offsets...\n", __func__);
        }
        Lv1ScanToFillOffests(globals);
        local_timer.stop();

        if (sdbg_builder_verbose >= 3) {
            log("[B::%s] Number of large differentials: %llu\n", __func__, globals.lv1_items_special.size());
            log("[B::%s] Lv1 scanning time elapsed: %.4lfs\n", __func__, local_timer.elapsed());
        }

        if (globals.lv1_items_special.size() > kDifferentialLimit) {
            err("[ERROER B::%s]Too many large differentials!\n", __func__);
            exit(1);
        }

        int lv2_iteration = 0;
        //======================================== LEVEL 2 loop ==========================================//
        globals.lv2_start_bucket = globals.lv1_start_bucket;
        while (globals.lv2_start_bucket < globals.lv1_end_bucket) {
            lv2_iteration++;
            // --- finds the bucket range for this iteration ---
            local_timer.reset();
            local_timer.start();
            globals.lv2_end_bucket = FindEndBucket(globals.bucket_sizes, globals.lv2_start_bucket, globals.lv1_end_bucket, globals.max_lv2_items, globals.lv2_num_items);
            if (globals.lv2_num_items == 0) { // i.e. can't even hold a single bucket
                err("[ERROR B::%s] Bucket %d too large for lv.2: contains %lld items\n", __func__, globals.lv2_end_bucket, globals.bucket_sizes[globals.lv2_end_bucket]);
                exit(1);
            }

            if (sdbg_builder_verbose >= 4) {
                log("[B::%s] > Iteration [%d,%d], from bucket %d to %d\n", __func__, lv1_iteration, lv2_iteration, globals.lv2_start_bucket, globals.lv2_end_bucket-1);
            }

            Lv2ExtractSubstrings(globals);
            local_timer.stop();

            if (sdbg_builder_verbose >= 4) {
                log("[B::%s] Extracting substrings... done. Time elapsed: %.4lfs\n", __func__, local_timer.elapsed());
            }

            // --- sorting ---
#ifdef DISABLE_GPU
            omp_set_num_threads(globals.num_cpu_threads - globals.phase2_num_output_threads);
            local_timer.reset();
            local_timer.start();
            lv2_cpu_sort(globals.lv2_substrings, globals.permutation, globals.cpu_sort_space, globals.words_per_substring, globals.lv2_num_items);
            omp_set_num_threads(globals.num_cpu_threads);
            local_timer.stop();

            if (sdbg_builder_verbose >= 4) {
                log("[B::%s] Sorting substrings with CPU...Sorting time elapsed: %.4lfs\n", __func__, local_timer.elapsed());
            }
#else
            local_timer.reset();
            local_timer.start();
            lv2_gpu_sort(globals.lv2_substrings, globals.permutation, globals.words_per_substring, globals.lv2_num_items);
            local_timer.stop();

            if (sdbg_builder_verbose >= 4) {
                log("[B::%s] Sorting substrings with GPU...Sorting time elapsed: %.4lfs\n", __func__, local_timer.elapsed());
            }
#endif
            // --- output ---
            if (output_thread_created) {
                Lv2OutputJoin(globals);
            }

            globals.lv2_num_items_to_output = globals.lv2_num_items;
            std::swap(globals.lv2_substrings_to_output, globals.lv2_substrings);
            std::swap(globals.permutation_to_output, globals.permutation);

            Lv2Output(globals);
            output_thread_created = true;

            globals.lv2_start_bucket = globals.lv2_end_bucket;
        } // end LEVEL 2 loop

        globals.lv1_start_bucket = globals.lv1_end_bucket;
    } // end LEVEL 1 loop

    if (output_thread_created) {
        Lv2OutputJoin(globals);
    }
    timer.stop();

    if (sdbg_builder_verbose >= 3) {
        log("[B::%s] Done sorting! Time elapsed: %.4lf\n", __func__, timer.elapsed());
        log("[B::%s] Number of $ A C G T A- C- G- T-:\n", __func__);
    }

    for (int i = 0; i < 9; ++i) {
        log("%lld ", globals.num_chars_in_w[i]);
    }
    log("\n");

    // --- write tails ---
    fprintf(globals.output_f_file, "%lld\n", (long long)globals.total_number_edges);
    fprintf(globals.output_f_file, "%d\n", globals.kmer_k);
    fprintf(globals.output_f_file, "%lld\n", (long long)globals.num_dollar_nodes);


    if (sdbg_builder_verbose >= 2) {
        log("[B::%s] Total number of edges: %llu\n", __func__, globals.total_number_edges);
        log("[B::%s] Total number of ONEs: %llu\n", __func__, globals.num_ones_in_last);
        log("[B::%s] Total number of v$ edges: %llu\n", __func__, globals.num_dummy_edges);
        log("[B::%s] Total number of $v edges: %llu\n", __func__, globals.num_dollar_nodes);
    }

    ////////////////////////////////// Cleaning up... /////////////////////////////////
    Phase2Clean(globals);
}

}