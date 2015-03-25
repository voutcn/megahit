/*
 *  cx1_functions_1pass.cpp
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
#include <parallel/algorithm>

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

inline int GetReadLength(edge_word_t* read_p, int words_per_read, int mask) {
    return *(read_p + words_per_read - 1) & mask;
}

inline int GetReadLengthByID(int64_t id, global_data_t &globals) {
    return *(globals.packed_reads + (id + 1) * globals.words_per_read - 1) & globals.read_length_mask;
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

// helper: see whether two lv2 items have the same (k-1)-mer
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


namespace phase1 {
/**
 * @brief encode read_id and its offset in one int64_t
 */
inline int64_t EncodeOffset(int64_t read_id, int offset, int strand, int length_num_bits) {
    return (read_id << (length_num_bits + 1)) | (offset << 1) | strand;
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

    // --- allocate memory for is_solid bit_vector
    globals.num_k1_per_read = globals.max_read_length - globals.kmer_k;
    globals.is_solid.reset(globals.num_k1_per_read * globals.num_reads);
    globals.mem_packed_reads += DivCeiling(globals.num_k1_per_read * globals.num_reads, 8);

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
    KmerUint32 k_minus1_mer, rev_k_minus1_mer; // (k-1)-mer and its rc
    for (int64_t read_id = rp.rp_start_id; read_id < rp.rp_end_id; ++read_id, read_p += globals.words_per_read) {
        int read_length = GetReadLength(read_p, globals.words_per_read, globals.read_length_mask);
        if (read_length < globals.kmer_k + 1) { continue; }
        k_minus1_mer.init(read_p, globals.kmer_k - 1);
        rev_k_minus1_mer.clean();
        for (int i = 0; i < globals.kmer_k - 1; ++i) {
            rev_k_minus1_mer.Append(3 - ExtractNthChar(read_p, globals.kmer_k - 2 - i));
        }
        
        // the first one special handling
        bucket_sizes[k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar]++;
        bucket_sizes[rev_k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar]++;

        int last_char_offset = globals.kmer_k - 1;
        int c = ExtractNthChar(read_p, last_char_offset);
        k_minus1_mer.ShiftLeftAppend(c);
        rev_k_minus1_mer.ShiftRightAppend(3 - c);

        while (last_char_offset < read_length - 1) {
            int cmp = k_minus1_mer.cmp(rev_k_minus1_mer);
            if (cmp > 0) {
                bucket_sizes[rev_k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar]++;
            } else {
                bucket_sizes[k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar]++;
            }

            int c = ExtractNthChar(read_p, ++last_char_offset);
            k_minus1_mer.ShiftLeftAppend(c);
            rev_k_minus1_mer.ShiftRightAppend(3 - c);
        }

        // last one special handling
        bucket_sizes[k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar]++;
        bucket_sizes[rev_k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar]++;
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
    // --- compute offset bits ---
    {
        globals.offset_num_bits = 0;
        int len = 1;
        while (len - 1 < globals.max_read_length) {
            globals.offset_num_bits++;
            len *= 2;
        }
    }

    // --- initialize output mercy files ---
    globals.num_mercy_files = 1;
    while (globals.num_mercy_files * 10485760LL < globals.num_reads && globals.num_mercy_files < 64) {
        globals.num_mercy_files <<= 1;
    }
    if (sdbg_builder_verbose >= 3) {
        log("[B::%s] Number of files for mercy candidate reads: %d\n", __func__, globals.num_mercy_files);
    }
    globals.mercy_output_locks.reset(globals.num_mercy_files);
    for (int i = 0; i < globals.num_mercy_files; ++i) {
        char file_name[10240];
        sprintf(file_name, "%s.mercy_cand.%d", globals.output_prefix, i);
        globals.mercy_files.push_back(OpenFileAndCheck(file_name, "wb"));
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
    globals.max_lv2_items = std::max(globals.max_bucket_size, kMinLv2BatchSizeGPU);
 #endif
    // to count (k+1)-mers, sort by the internal (k-1)-mer
    // (k+1)-mer = abS[0..k-2]cd
    // is solid: number of bSc >= threshold
    // bS has in coming: for some a, num of abS >= threshold
    // Sc has outgoing: for some a, num of Scd >= threshold
    globals.words_per_substring = DivCeiling((globals.kmer_k - 1) * kBitsPerEdgeChar + 2 * kBWTCharNumBits, kBitsPerEdgeWord);
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
 #ifdef DISABLE_GPU
    globals.cpu_sort_space = (uint64_t*) MallocAndCheck(sizeof(uint64_t) * globals.max_lv2_items, __FILE__, __LINE__);
 #endif
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
    KmerUint32 k_minus1_mer, rev_k_minus1_mer; // (k+1)-mer and its rc
    int key;
    for (int64_t read_id = rp.rp_start_id; read_id < rp.rp_end_id; ++read_id, read_p += globals.words_per_read) {
        int read_length = GetReadLength(read_p, globals.words_per_read, globals.read_length_mask);
        if (read_length < globals.kmer_k + 1) { continue; }
        k_minus1_mer.init(read_p, globals.kmer_k - 1);
        rev_k_minus1_mer.clean();
        for (int i = 0; i < globals.kmer_k - 1; ++i) {
            rev_k_minus1_mer.Append(3 - ExtractNthChar(read_p, globals.kmer_k - 2 - i));
        }

        // ===== this is a macro to save some copy&paste ================
 #define CHECK_AND_SAVE_OFFSET(offset, strand)                                   \
    do {                                                                \
      assert(offset + globals.kmer_k - 1 <= read_length); \
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

        // the first one special handling
        key = k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar;
        CHECK_AND_SAVE_OFFSET(0, 0);
        key = rev_k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar;
        CHECK_AND_SAVE_OFFSET(0, 1);

        int last_char_offset = globals.kmer_k - 1;
        int c = ExtractNthChar(read_p, last_char_offset);
        k_minus1_mer.ShiftLeftAppend(c);
        rev_k_minus1_mer.ShiftRightAppend(3 - c);

        // shift the key char by char
        while (last_char_offset < read_length - 1) {
            int cmp = k_minus1_mer.cmp(rev_k_minus1_mer);
            if (cmp > 0) {
                key = rev_k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar;
                CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k + 2, 1);
            } else if (cmp < 0) {
                key = k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar;
                CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k + 2, 0);
            } else {
                int prev = ExtractNthChar(read_p, last_char_offset - (globals.kmer_k - 1));
                int next = ExtractNthChar(read_p, last_char_offset + 1);
                if (prev <= 3 - next) {
                    key = rev_k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar;
                    CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k + 2, 0);
                } else {
                    key = rev_k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar;
                    CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k + 2, 1);
                }
            }

            int c = ExtractNthChar(read_p, ++last_char_offset);
            k_minus1_mer.ShiftLeftAppend(c);
            rev_k_minus1_mer.ShiftRightAppend(3 - c);
        }

        // the last one special handling
        key = k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar;
        CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k + 2, 0);
        key = rev_k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar;
        CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k + 2, 1);
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
    }\
    // revert rp_bucket_offsets
    Lv1ComputeBucketOffset(globals);
}

// single thread helper function
// 'spacing' is the strip length for read-word "coalescing"
void CopySubstring(edge_word_t* dest, edge_word_t* src_read, int offset, int num_chars_to_copy, global_data_t &globals, uint8_t head, uint8_t tail) {
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

    edge_word_t *last_word = dest + (globals.words_per_substring - 1) * spacing;
    *last_word |= (head << kBWTCharNumBits) | tail;
}

void CopySubstringRC(edge_word_t* dest, edge_word_t* src_read, int offset, int num_chars_to_copy, global_data_t &globals, uint8_t head, uint8_t tail) {
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

    edge_word_t *last_word = dest + (globals.words_per_substring - 1) * spacing;
    *last_word |= (head << kBWTCharNumBits) | tail;
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

                int num_chars_to_copy = globals.kmer_k - 1;
                unsigned char prev, next, head, tail; // (k+1)=abScd, prev=a, head=b, tail=c, next=d
                if (offset > 1) {
                    head = ExtractNthChar(PACKED_READS(read_id, globals), offset - 1);
                    prev = ExtractNthChar(PACKED_READS(read_id, globals), offset - 2);
                } else {
                    prev = kSentinelValue;
                    if (offset > 0) {
                        head = ExtractNthChar(PACKED_READS(read_id, globals), offset - 1);
                    } else {
                        head = kSentinelValue;
                    }
                }

                int read_length = GetReadLengthByID(read_id, globals);
                if (offset + globals.kmer_k < read_length) {
                    tail = ExtractNthChar(PACKED_READS(read_id, globals), offset + globals.kmer_k - 1);
                    next = ExtractNthChar(PACKED_READS(read_id, globals), offset + globals.kmer_k);
                } else {
                    next = kSentinelValue;
                    if (offset + globals.kmer_k - 1 < read_length) {
                        tail = ExtractNthChar(PACKED_READS(read_id, globals), offset + globals.kmer_k - 1);
                    } else {
                        tail = kSentinelValue;
                    }
                }

                if (strand == 0) {
                    CopySubstring(substrings_p, PACKED_READS(read_id, globals), offset, num_chars_to_copy, globals, head, tail);
                    *read_info_p = (full_offset << 6) | (prev << 3) | next;
                } else {
                    CopySubstringRC(substrings_p, PACKED_READS(read_id, globals), offset, num_chars_to_copy, globals, 
                                    tail == kSentinelValue ? kSentinelValue : (3 - tail), head == kSentinelValue ? kSentinelValue : (3 - head));
                    *read_info_p = (full_offset << 6) | ((next == kSentinelValue ? kSentinelValue : (3 - next)) << 3)
                                                      | (prev == kSentinelValue ? kSentinelValue : (3 - prev));
                }

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
}

// helper
inline uint8_t ExtractHeadTail(edge_word_t *item, int64_t spacing, int words_per_substring) {
    return *(item + spacing * (words_per_substring - 1)) & ((1 << 2 * kBWTCharNumBits) - 1);
}

inline uint8_t ExtractPrevNext(int i, global_data_t &globals) {
    return globals.lv2_read_info_to_output[i] & ((1 << 2 * kBWTCharNumBits) - 1);
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
    int count_prev_head[5][5];
    int count_tail_next[5][5];
    int64_t offset_mask = (1 << globals.offset_num_bits) - 1;
    int64_t *thread_edge_counting = globals.thread_edge_counting + thread_id * (kMaxMulti_t + 1);

    for (int i = op_start_index; i < op_end_index; i = end_idx) {
        start_idx = i;
        end_idx = i + 1;
        edge_word_t *first_item = globals.lv2_substrings_to_output + (globals.permutation_to_output[i]);

        while (end_idx < op_end_index) {
            if (IsDiffKMinusOneMer(first_item,
                                 globals.lv2_substrings_to_output + globals.permutation_to_output[end_idx],
                                 globals.lv2_num_items_to_output,
                                 globals.kmer_k))
            {
                break;
            }
            ++end_idx;
        }

        memset(count_prev_head, 0, sizeof(count_prev_head));
        memset(count_tail_next, 0, sizeof(count_tail_next));
        for (int j = start_idx; j < end_idx; ++j) {
            uint8_t prev_and_next = ExtractPrevNext(globals.permutation_to_output[j], globals);
            uint8_t head_and_tail = ExtractHeadTail(globals.lv2_substrings_to_output + globals.permutation_to_output[j], globals.lv2_num_items_to_output, globals.words_per_substring);
            count_prev_head[prev_and_next >> 3][head_and_tail >> 3]++;
            count_tail_next[head_and_tail & 7][prev_and_next & 7]++;
        }

        int has_in = 0, has_out = 0;
        for (int j = 0; j < 4; ++j) {
            for (int x = 0; x < 4; ++x) {
                if (count_prev_head[x][j] >= globals.kmer_freq_threshold) {
                    has_in |= 1 << j;
                    break;
                }
            }

            for (int x = 0; x < 4; ++x) {
                if (count_tail_next[j][x] >= globals.kmer_freq_threshold) {
                    has_out |= 1 << j;
                    break;
                }
            }
        }

        while (i < end_idx) {
            uint8_t head_and_tail = ExtractHeadTail(globals.lv2_substrings_to_output + globals.permutation_to_output[i], globals.lv2_num_items_to_output, globals.words_per_substring);
            uint8_t head = head_and_tail >> 3;
            uint8_t tail = head_and_tail & 7;
            int count = 1;
            ++i;
            if (head == kSentinelValue || tail == kSentinelValue) {
                continue;
            }

            while (i < end_idx && ExtractHeadTail(globals.lv2_substrings_to_output + globals.permutation_to_output[i], globals.lv2_num_items_to_output, globals.words_per_substring) == head_and_tail) {
                ++i;
                ++count;
            }

            ++thread_edge_counting[std::min(count, kMaxMulti_t)];
            if (count < globals.kmer_freq_threshold) { continue; }

            for (int j = i - count; j < i; ++j) {
                int64_t read_info = globals.lv2_read_info_to_output[globals.permutation_to_output[j]] >> 6;
                int strand = read_info & 1;
                int offset = ((read_info >> 1) & offset_mask) - 1;
                int64_t read_id = read_info >> (1 + globals.offset_num_bits);

                // mark this is a solid edge
                globals.is_solid.set(globals.num_k1_per_read * read_id + offset); 

                // check if this is a mercy candidate
                if ((!(has_in & (1 << head)) && strand == 0) || (!(has_out & (1 << tail)) && strand == 1)) {
                    assert(offset < globals.num_k1_per_read);
                    // no in
                    int64_t packed_mercy_cand = (read_id << (globals.offset_num_bits + 1)) | (offset << 1);
                    while (globals.mercy_output_locks.lock(read_id & (globals.num_mercy_files - 1))) {
                        continue;
                    }
                    fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1, globals.mercy_files[read_id & (globals.num_mercy_files - 1)]);
                    globals.mercy_output_locks.unset(read_id & (globals.num_mercy_files - 1));
                }

                if ((!(has_in & (1 << head)) && strand == 1) || (!(has_out & (1 << tail)) && strand == 0)) {
                    assert(offset < globals.num_k1_per_read);
                    // no out
                    int64_t packed_mercy_cand = (read_id << (globals.offset_num_bits + 1)) | (offset << 1) | 1;
                    while (globals.mercy_output_locks.lock(read_id & (globals.num_mercy_files - 1))) {
                        continue;
                    }
                    fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1, globals.mercy_files[read_id & (globals.num_mercy_files - 1)]);
                    globals.mercy_output_locks.unset(read_id & (globals.num_mercy_files - 1));
                }
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
                if (IsDiffKMinusOneMer(prev_item, item, globals.lv2_num_items_to_output, globals.kmer_k)) {
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
    // free(globals.packed_reads);
    free(globals.bucket_sizes);
    free(globals.lv1_items);
    free(globals.lv2_substrings);
    free(globals.permutation);
    free(globals.lv2_substrings_to_output);
    free(globals.permutation_to_output);
    free(globals.lv2_read_info);
    free(globals.lv2_read_info_to_output);
    free(globals.edge_counting);
    free(globals.thread_edge_counting);
    for (int t = 0; t < globals.num_cpu_threads; ++t) {
       free(globals.readpartitions[t].rp_bucket_sizes);
       free(globals.readpartitions[t].rp_bucket_offsets);
    }

    for (int i = 0; i < globals.num_mercy_files; ++i) {
        fclose(globals.mercy_files[i]);
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

inline int64_t EncodeOffset(int64_t read_id, int offset, int strand, int length_num_bits, int edge_type) {
    // edge_type: 0 left $; 1 solid; 2 right $
    return (read_id << (length_num_bits + 3)) | (offset << 3) | (edge_type << 1) | strand;
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
        rp.rp_lv1_differential_base = EncodeOffset(rp.rp_start_id, 0, 0, globals.offset_num_bits, 0);
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
        int64_t full_offset = globals.num_k1_per_read * read_id;
        while (true) {
            if (globals.is_solid.get(full_offset)) {
                bool is_palindrome = (rev_edge == edge);
                bucket_sizes[(edge.data_[0] << 2) >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar]++;
                if (!is_palindrome)
                    bucket_sizes[(rev_edge.data_[0] << 2) >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar]++;

                if (last_char_offset == globals.kmer_k || !globals.is_solid.get(full_offset - 1)) {
                    bucket_sizes[edge.data_[0] >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar]++;
                    if (!is_palindrome)
                        bucket_sizes[(rev_edge.data_[0] << 4) >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar]++;
                }

                if (last_char_offset == read_length - 1 || !globals.is_solid.get(full_offset + 1)) {
                    bucket_sizes[(edge.data_[0] << 4) >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar]++;
                    if (!is_palindrome)
                        bucket_sizes[rev_edge.data_[0] >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar]++;
                }
            }

            ++full_offset;

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

void AddMercyEdge(global_data_t &globals) {
    std::vector<uint64_t> mercy_cand;
    uint64_t offset_mask = (1 << globals.offset_num_bits) - 1; // 0000....00011..11
    uint64_t num_mercy = 0;

    for (int fid = 0; fid < globals.num_mercy_files; ++fid) {
        char file_name[10240];
        sprintf(file_name, "%s.mercy_cand.%d", globals.output_prefix, fid);
        FILE *fp = OpenFileAndCheck(file_name, "rb");
        mercy_cand.clear();

        int num_read = 0;
        uint64_t buf[4096];
        while ((num_read = fread(buf, sizeof(uint64_t), 4096, fp)) > 0) {
            mercy_cand.insert(mercy_cand.end(), buf, buf + num_read);
        }

        omp_set_num_threads(globals.num_cpu_threads);
        __gnu_parallel::sort(mercy_cand.begin(), mercy_cand.end());

        // multi threading
        uint64_t avg = DivCeiling(mercy_cand.size(), globals.num_cpu_threads);
        std::vector<uint64_t> start_idx(globals.num_cpu_threads), end_idx(globals.num_cpu_threads);

        // manually distribute threads
        for (int tid = 0; tid < globals.num_cpu_threads; ++tid) {
            if (tid == 0) { start_idx[tid] = 0; }
            else { start_idx[tid] = end_idx[tid - 1]; }

            uint64_t this_end = avg * (tid + 1);
            uint64_t read_id = mercy_cand[this_end] >> (globals.offset_num_bits + 1);
            while (this_end < mercy_cand.size() && (mercy_cand[this_end] >> (globals.offset_num_bits + 1)) == read_id) {
                ++this_end;
            }
            end_idx[tid] = this_end;
        }

 #pragma omp parallel for reduction(+:num_mercy)
        for (int tid = 0; tid < globals.num_cpu_threads; ++tid) {
            std::vector<bool> no_in(globals.max_read_length);
            std::vector<bool> no_out(globals.max_read_length);

            uint64_t i = start_idx[tid];
            // go read by read
            while (i != end_idx[tid]) {
                uint64_t read_id = mercy_cand[i] >> (globals.offset_num_bits + 1);
                int first_0_out = globals.max_read_length + 1;
                int last_0_in = -1;

                std::fill(no_in.begin(), no_in.end(), false);
                std::fill(no_out.begin(), no_out.end(), false);

                while (i != end_idx[tid] && (mercy_cand[i] >> (globals.offset_num_bits + 1)) == read_id) {
                    if (mercy_cand[i] & 1) {
                        no_out[(mercy_cand[i] >> 1) & offset_mask] = true;
                        first_0_out = std::min(first_0_out, int((mercy_cand[i] >> 1) & offset_mask));
                    } else {
                        no_in[(mercy_cand[i] >> 1) & offset_mask] = true;
                        last_0_in = std::max(last_0_in, int((mercy_cand[i] >> 1) & offset_mask));
                    }
                    ++i;
                }
                if (last_0_in < first_0_out) { continue; }

                int read_length = GetReadLengthByID(read_id, globals);
                int last_no_out = -1;

                for (int i = 0; i + globals.kmer_k < read_length; ++i) {
                    if (no_in[i] && last_no_out != -1) {
                        assert(globals.is_solid.get(read_id * globals.num_k1_per_read + i));
                        for (int j = last_no_out + 1; j < i; ++j) {
                            globals.is_solid.set(read_id * globals.num_k1_per_read + j);
                        }
                        num_mercy += i - last_no_out - 1;
                    }
                    if (globals.is_solid.get(read_id * globals.num_k1_per_read + i)) {
                        last_no_out = -1;
                    }
                    if (no_out[i]) {
                        assert(globals.is_solid.get(read_id * globals.num_k1_per_read + i));
                        last_no_out = i;
                    }
                }

            }
        }

        fclose(fp);
    }

    log("[B::%s] Number of mercy edges (reads): %ld\n", __func__, num_mercy);
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
    globals.max_lv2_items = std::max(globals.max_bucket_size, kMinLv2BatchSizeGPU);
 #endif
    globals.words_per_substring = DivCeiling(globals.kmer_k * kBitsPerEdgeChar + kBWTCharNumBits + 1, kBitsPerEdgeWord);
    globals.words_per_dummy_node = DivCeiling(globals.kmer_k * kBitsPerEdgeChar, kBitsPerEdgeWord);
    // lv2 bytes: substring (double buffer), permutation, aux
    int64_t lv2_bytes_per_item = (globals.words_per_substring * sizeof(edge_word_t) + sizeof(uint32_t)) * 2 + sizeof(int64_t);
 #ifdef DISABLE_GPU
    lv2_bytes_per_item += sizeof(uint64_t) * 2; // simulate GPU
 #endif

    if (sdbg_builder_verbose >= 2) {
        log("[B::%s] %d words per substring, words per dummy node ($v): %d\n", __func__, globals.words_per_substring, globals.words_per_dummy_node);
    }

    // --- memory stuff ---
    int64_t mem_remained = globals.host_mem
                         - globals.mem_packed_reads
                         - phase1::kNumBuckets * sizeof(int64_t) * (globals.num_cpu_threads * 3 + 1);
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
        log("[B::%s] Memory for edges: %lld\n", __func__, globals.mem_packed_reads);
        log("[B::%s] max # lv.1 items = %lld\n", __func__, globals.max_lv1_items);
        log("[B::%s] max # lv.2 items = %lld\n", __func__, globals.max_lv2_items);
    }

    // --- alloc memory ---
    globals.lv1_items = (int*) MallocAndCheck(globals.max_lv1_items * sizeof(int), __FILE__, __LINE__);
    globals.lv2_substrings = (edge_word_t*) MallocAndCheck(globals.max_lv2_items * globals.words_per_substring * sizeof(edge_word_t), __FILE__, __LINE__);
    globals.permutation = (uint32_t *) MallocAndCheck(globals.max_lv2_items * sizeof(uint32_t), __FILE__, __LINE__);
    globals.lv2_substrings_to_output = (edge_word_t*) MallocAndCheck(globals.max_lv2_items * globals.words_per_substring * sizeof(edge_word_t), __FILE__, __LINE__);
    globals.permutation_to_output = (uint32_t *) MallocAndCheck(globals.max_lv2_items * sizeof(uint32_t), __FILE__, __LINE__);
 #ifdef DISABLE_GPU
    globals.cpu_sort_space = (uint64_t*) MallocAndCheck(sizeof(uint64_t) * globals.max_lv2_items, __FILE__, __LINE__); // simulate GPU
 #endif
    globals.lv2_output_items.resize(globals.phase2_num_output_threads);

    // --- write header ---
    fprintf(globals.output_f_file, "-1\n");
    globals.dummy_nodes_writer.output(globals.words_per_dummy_node);
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
 #define CHECK_AND_SAVE_OFFSET(offset, strand, edge_type)                                   \
    do {                                                                \
      assert(edge_type < 3 && edge_type >= 0 && strand >= 0 && strand <= 1 && offset + globals.kmer_k < read_length); \
      if (((key - globals.lv1_start_bucket) ^ (key - globals.lv1_end_bucket)) & kSignBitMask) { \
        int64_t full_offset = EncodeOffset(read_id, offset, strand, globals.offset_num_bits, edge_type); \
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
        int64_t full_offset = globals.num_k1_per_read * read_id;
        while (true) {
            if (globals.is_solid.get(full_offset)) {
                bool is_palindrome = (rev_edge == edge);

                // left $
                if (last_char_offset == globals.kmer_k || !globals.is_solid.get(full_offset - 1)) {
                    key = edge.data_[0] >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar;
                    CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 0, 0);
                    if (!is_palindrome) {
                        key = (rev_edge.data_[0] << 4) >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar;
                        CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 1, 0);
                    }
                }

                // solid
                key = (edge.data_[0] << 2) >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar;
                CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 0, 1);

                if (!is_palindrome) {
                    key = (rev_edge.data_[0] << 2) >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar;
                    CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 1, 1);
                }

                // right $
                if (last_char_offset == read_length - 1 || !globals.is_solid.get(full_offset + 1)) {
                    key = (edge.data_[0] << 4) >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar;
                    CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 0, 2);
                    if (!is_palindrome) {
                        key = rev_edge.data_[0] >> (kCharsPerEdgeWord - phase1::kBucketPrefixLength) * kBitsPerEdgeChar;
                        CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 1, 2);
                    }
                }
            }

            ++full_offset;

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
void CopySubstring(edge_word_t* dest, edge_word_t* src_read, int offset, int num_chars_to_copy, global_data_t &globals, uint8_t prev_char = 0) {
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

    edge_word_t *last_word = dest + (words_per_substring - 1) * spacing;
    *last_word |= int(num_chars_to_copy == globals.kmer_k) << kBWTCharNumBits;
    *last_word |= prev_char;
}

void CopySubstringRC(edge_word_t* dest, edge_word_t* src_read, int offset, int num_chars_to_copy, global_data_t &globals, uint8_t prev_char = 0) {
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

    edge_word_t *last_word = dest + (globals.words_per_substring - 1) * spacing;
    *last_word |= int(num_chars_to_copy == globals.kmer_k) << kBWTCharNumBits;
    *last_word |= prev_char;
}

void* Lv2ExtractSubstringsThread(void* _data) {
    struct bucketpartition_data_t &bp = *((struct bucketpartition_data_t*) _data);
    struct global_data_t &globals = *(bp.globals);
    int *lv1_p = globals.lv1_items + globals.readpartitions[0].rp_bucket_offsets[ bp.bp_start_bucket ];
    int64_t offset_mask = (1 << globals.offset_num_bits) - 1; // 0000....00011..11
    edge_word_t *substrings_p = globals.lv2_substrings +
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
                int64_t read_id = full_offset >> (globals.offset_num_bits + 3);
                int offset = (full_offset >> 3) & offset_mask;
                int strand = full_offset & 1;
                int edge_type = (full_offset >> 1) & 3;

                if (strand == 0) {
                    int num_chars_to_copy = globals.kmer_k;
                    uint8_t prev = kSentinelValue;

                    switch (edge_type) {
                      case 0:
                        break;
                      case 1:
                        prev = ExtractNthChar(PACKED_READS(read_id, globals), offset);
                        offset++;
                        break;
                      case 2:
                        prev = ExtractNthChar(PACKED_READS(read_id, globals), offset + 1);
                        offset += 2;
                        num_chars_to_copy--;
                        break;
                      default:
                        assert(false);
                    }

                    CopySubstring(substrings_p, PACKED_READS(read_id, globals), offset, num_chars_to_copy, globals, prev);
                } else {
                    int num_chars_to_copy = globals.kmer_k;
                    uint8_t prev = kSentinelValue;

                    switch (edge_type) {
                      case 0:
                        num_chars_to_copy--;
                        prev = 3 - ExtractNthChar(PACKED_READS(read_id, globals), offset + globals.kmer_k - 1);
                        break;
                      case 1:
                        prev = 3 - ExtractNthChar(PACKED_READS(read_id, globals), offset + globals.kmer_k);
                        break;
                      case 2:
                        offset++;
                        break;
                      default:
                        assert(false);
                    }

                    CopySubstringRC(substrings_p, PACKED_READS(read_id, globals), offset, num_chars_to_copy, globals, prev);
                }

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
}

// helper
inline int ExtractFirstChar(edge_word_t *item) {
    return *item >> kTopCharShift;
}

// bS'a
inline int Extract_a(edge_word_t *item, int num_words, int64_t spacing, int kmer_k) {
    int non_dollar = (item[(num_words - 1) * spacing] >> kBWTCharNumBits) & 1;
    if (non_dollar) {
        int which_word = (kmer_k - 1) / kCharsPerEdgeWord;
        int word_index = (kmer_k - 1) % kCharsPerEdgeWord;
        return (item[which_word * spacing] >> (kCharsPerEdgeWord - 1 - word_index) * kBitsPerEdgeChar) & kEdgeCharMask;
    } else {
        return kSentinelValue;
    }
}

inline int Extract_b(edge_word_t *item, int num_words, int64_t spacing) {
    return item[(num_words - 1) * spacing] & ((1 << kBWTCharNumBits) - 1);
}

void *Lv2OutputThread(void *_op) {
    struct outputpartition_data_t *op = (struct outputpartition_data_t*) _op;
    struct global_data_t &globals = *(op->globals);
    int64_t op_start_index = op->op_start_index;
    int64_t op_end_index = op->op_end_index;
    int thread_id = op->op_id;
    int start_idx, end_idx;
    int has_solid_a = 0; // has solid (k+1)-mer aSb
    int has_solid_b = 0; // has solid aSb
    int last_a[4], outputed_b;

    globals.lv2_output_items[thread_id].clear();

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
            int count = std::min(j - i, kMaxMulti_t);

            if (a == kSentinelValue) {
                assert(b != kSentinelValue);
                if (has_solid_b & (1 << b)) {
                    continue;
                }
                is_dollar = 1;
                count = 0;
            }

            if (b == kSentinelValue) {
                assert(a != kSentinelValue);
                if (has_solid_a & (1 << a)) {
                    continue;
                }
                count = 0;
            }

            w = (b == kSentinelValue) ? 0 : ((outputed_b & (1 << b)) ? b + 5 : b + 1);
            outputed_b |= 1 << b;
            last = (a == kSentinelValue) ? 0 : ((last_a[a] == j - 1) ? 1 : 0);

            // save this item to the out_item array
            uint64_t out_item = (uint64_t(i) << 32) | (count << 16) | (is_dollar << 5) | (last << 4) | w;
            globals.lv2_output_items[thread_id].push_back(out_item);
        }
    }

    pthread_barrier_wait(&globals.output_barrier);

    if (op_start_index == 0) {
        xtimer_t local_timer;
        local_timer.reset();
        local_timer.start();
        for (int tid = 0; tid < globals.phase2_num_output_threads; ++tid) {
            for (auto it = globals.lv2_output_items[tid].begin(); it != globals.lv2_output_items[tid].end(); ++it) {
                int i = (*it) >> 32;
                edge_word_t *item = globals.lv2_substrings_to_output + globals.permutation_to_output[i];
                while (ExtractFirstChar(item) > globals.cur_suffix_first_char) {
                    ++globals.cur_suffix_first_char;
                    fprintf(globals.output_f_file, "%lld\n", (long long)globals.total_number_edges);
                }

                multi_t counting_to_output = (*it) >> 16 & kMaxMulti_t;
                // output
                globals.sdbg_writer.outputW(*it & 0xF);
                globals.sdbg_writer.outputLast((*it >> 4) & 1);
                globals.sdbg_writer.outputIsDollar((*it >> 5) & 1);
                if (counting_to_output <= kMaxMulti2_t) {
                    multi2_t c = counting_to_output;
                    fwrite(&c, sizeof(multi2_t), 1, globals.output_multiplicity_file);   
                } else {
                    int64_t c = counting_to_output | (globals.total_number_edges << 16);
                    fwrite(&c, sizeof(int64_t), 1, globals.output_multiplicity_file2);
                    fwrite(&kMulti2Sp, sizeof(multi2_t), 1, globals.output_multiplicity_file);
                }

                // stat
                globals.total_number_edges++;
                globals.num_chars_in_w[*it & 0xF]++;
                globals.num_ones_in_last += (*it >> 4) & 1;

                if ((*it >> 5) & 1) {
                    globals.num_dollar_nodes++;
                    if (globals.num_dollar_nodes >= phase2::kMaxDummyEdges) {
                        err("[ERROR B::%s] Too many dummy nodes (>= %lld)! The graph contains too many tips!\n", __func__, (long long)phase2::kMaxDummyEdges);
                        exit(1);
                    }
                    for (int64_t i = 0; i < globals.words_per_dummy_node; ++i) {
                        globals.dummy_nodes_writer.output(item[i * globals.lv2_num_items_to_output]);
                    }
                }
                if ((*it & 0xF) == 0) {
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
    free(globals.packed_reads);
    free(globals.bucket_sizes);
    free(globals.lv1_items);
    free(globals.lv2_substrings);
    free(globals.permutation);
    free(globals.lv2_substrings_to_output);
    free(globals.permutation_to_output);
    fclose(globals.output_f_file);
    fclose(globals.output_multiplicity_file);
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

    InitDNAMap();

    if (globals.need_mercy) {
        timer.reset();
        timer.start();
        if (sdbg_builder_verbose >= 3) {
            log("[B::%s] Adding mercy edges...\n", __func__);
        }
        AddMercyEdge(globals);
        timer.stop();

        if (sdbg_builder_verbose >= 3) {
            log("[B::%s] Done. Time elapsed: %.4lfs\n", __func__, timer.elapsed());
        }
    }

    // --- init global data ---
    InitGlobalData(globals);

    ////////////////////////////////// Start processing... ////////////////////////////
    int lv1_iteration = 0;
    //======================================== LEVEL 1 loop ============================================//
    globals.lv1_start_bucket = 0;
    timer.reset();
    timer.start();
    bool output_thread_created = false;
    while (globals.lv1_start_bucket < phase1::kNumBuckets) {
        xtimer_t local_timer;
        lv1_iteration++;
        // --- finds the bucket range for this iteration ---

        if (sdbg_builder_verbose >= 3) {
            log("[B::%s] Finding end bucket...\n", __func__);
        }
        local_timer.reset();
        local_timer.start();
        globals.lv1_end_bucket = FindEndBucket(globals.bucket_sizes, globals.lv1_start_bucket, phase1::kNumBuckets, globals.max_lv1_items, globals.lv1_num_items);
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