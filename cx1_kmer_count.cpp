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


#include "cx1_kmer_count.h"

#include <string.h>
#include <algorithm>
#include <zlib.h>
#include <omp.h>

#include "mem_file_checker-inl.h"
#include "kseq.h"
#include "utils.h"
#include "kmer.h"
#include "packed_reads.h"
#include "read_lib_functions-inl.h"

#include "lv2_cpu_sort.h"
#include "lv2_gpu_functions.h"

extern void kt_dfor(int n_threads, void (*func)(void *, long, int), void *data, long n);

namespace cx1_kmer_count {

// helpers
typedef CX1<count_global_t, kNumBuckets> cx1_t;
typedef CX1<count_global_t, kNumBuckets>::readpartition_data_t readpartition_data_t;
typedef CX1<count_global_t, kNumBuckets>::bucketpartition_data_t bucketpartition_data_t;
typedef CX1<count_global_t, kNumBuckets>::outputpartition_data_t outputpartition_data_t;

/**
 * @brief encode read_id and its offset in one int64_t
 */
inline int64_t EncodeOffset(int64_t read_id, int offset, int strand, SequencePackage &p) {
    return ((p.get_start_index(read_id) + offset) << 1) | strand;
}

inline bool IsDifferentEdges(uint32_t *item1, uint32_t *item2, int num_words, int spacing) {
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
inline void PackEdge(uint32_t *dest, uint32_t *item, int counting, struct count_global_t &globals, int num_items) {
    for (int i = 0; i < globals.words_per_edge && i < globals.words_per_substring; ++i) {
        dest[i] = *(item + (int64_t)i * num_items);
    }

    int chars_in_last_word = (globals.kmer_k + 1) % kCharsPerEdgeWord;
    int which_word = (globals.kmer_k + 1) / kCharsPerEdgeWord;

    if (chars_in_last_word > 0) {
        dest[which_word] >>= (kCharsPerEdgeWord - chars_in_last_word) * 2;
        dest[which_word] <<= (kCharsPerEdgeWord - chars_in_last_word) * 2;
    }
    else {
        dest[which_word] = 0;
    }

    while (++which_word < globals.words_per_edge) {
        dest[which_word] = 0;
    }

    dest[globals.words_per_edge - 1] |= std::min(kMaxMulti_t, counting);
}

// function pass to CX1

int64_t encode_lv1_diff_base(int64_t read_id, count_global_t &g) {
    return EncodeOffset(read_id, 0, 0, g.package);
}

void read_input_prepare(count_global_t &globals) { // num_items_, num_cpu_threads_ and num_output_threads_ must be set here
    bool is_reverse = true;

    int64_t num_bases, num_reads;
    GetBinaryLibSize(globals.read_lib_file, num_bases, num_reads);

    if (globals.assist_seq_file != "") {
        FILE *assist_seq_info = OpenFileAndCheck((globals.assist_seq_file + ".info").c_str(), "r");
        long long num_ass_bases, num_ass_seq;
        assert(fscanf(assist_seq_info, "%lld%lld", &num_ass_seq, &num_ass_bases) == 2);
        fclose(assist_seq_info);

        num_bases += num_ass_bases;
        num_reads += num_ass_seq;
    }
    
    globals.package.reserve_num_seq(num_reads);
    globals.package.reserve_bases(num_bases);

    ReadBinaryLibs(globals.read_lib_file, globals.package, globals.lib_info, is_reverse);

    if (globals.assist_seq_file != "") {
        SequenceManager seq_manager;
        seq_manager.set_file_type(SequenceManager::kFastxReads);
        seq_manager.set_file(globals.assist_seq_file);

        bool reverse_read = true;
        bool append_to_package = true;
        bool trimN = false;

        seq_manager.ReadShortReads(1LL << 60, 1LL << 60, append_to_package, reverse_read, trimN);
        seq_manager.clear();
    }
    
    globals.package.BuildLookup();
    globals.max_read_length = globals.package.max_read_len();
    globals.num_reads = globals.package.size();

    xlog("%ld reads, %d max read length\n", globals.num_reads, globals.max_read_length);

    // calc words_per_xxx
    int bits_read_length = 1; // bit needed to store read_length

    while ((1 << bits_read_length) - 1 < globals.max_read_length) {
        ++bits_read_length;
    }

    globals.offset_num_bits = bits_read_length;
    globals.read_length_mask = (1 << bits_read_length) - 1;

    globals.mem_packed_reads = globals.package.size_in_byte();

    int64_t mem_low_bound = globals.mem_packed_reads
                            + globals.num_reads * sizeof(unsigned char) * 2 // first_in0 & last_out0
                            + kNumBuckets * sizeof(int64_t) * (globals.num_cpu_threads * 3 + 1)
                            + (kMaxMulti_t + 1) * (globals.num_output_threads + 1) * sizeof(int64_t);
    mem_low_bound *= 1.05;

    if (mem_low_bound > globals.host_mem) {
        xerr_and_exit("%lld bytes is not enough for CX1 sorting, please set -m parameter to at least %lld\n", globals.host_mem, mem_low_bound);
    }

    // set cx1 param
    globals.cx1.num_cpu_threads_ = globals.num_cpu_threads;
    globals.cx1.num_output_threads_ = globals.num_output_threads;
    globals.cx1.num_items_ = globals.num_reads;
}

void *lv0_calc_bucket_size(void *_data) {
    readpartition_data_t &rp = *((readpartition_data_t *) _data);
    count_global_t &globals = *(rp.globals);
    int64_t *bucket_sizes = rp.rp_bucket_sizes;
    memset(bucket_sizes, 0, sizeof(bucket_sizes[0]) * kNumBuckets);
    GenericKmer edge, rev_edge; // (k+1)-mer and its rc

    for (int64_t read_id = rp.rp_start_id; read_id < rp.rp_end_id; ++read_id) {
        int read_length = globals.package.length(read_id);

        if (read_length < globals.kmer_k + 1) {
            continue;
        }

        uint32_t *read_p = &globals.package.packed_seq[globals.package.get_start_index(read_id) / 16];
        edge.init(read_p, globals.package.get_start_index(read_id) % 16, globals.kmer_k + 1);
        rev_edge = edge;
        rev_edge.ReverseComplement(globals.kmer_k + 1);

        int last_char_offset = globals.kmer_k;

        while (true) {
            if (rev_edge.cmp(edge, globals.kmer_k + 1) < 0) {
                bucket_sizes[rev_edge.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
            }
            else {
                bucket_sizes[edge.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
            }

            if (++last_char_offset >= read_length) {
                break;
            }
            else {
                int c = globals.package.get_base(read_id, last_char_offset);
                edge.ShiftAppend(c, globals.kmer_k + 1);
                rev_edge.ShiftPreappend(3 - c, globals.kmer_k + 1);
            }
        }
    }

    return NULL;
}

void init_global_and_set_cx1(count_global_t &globals) {
    // --- calculate lv2 memory ---
    globals.max_bucket_size = *std::max_element(globals.cx1.bucket_sizes_, globals.cx1.bucket_sizes_ + kNumBuckets);
    globals.tot_bucket_size = 0;
    int num_non_empty = 0;

    for (int i = 0; i < kNumBuckets; ++i) {
        globals.tot_bucket_size += globals.cx1.bucket_sizes_[i];

        if (globals.cx1.bucket_sizes_[i] > 0) {
            num_non_empty++;
        }
    }

    globals.words_per_substring = DivCeiling((globals.kmer_k + 1) * kBitsPerEdgeChar, kBitsPerEdgeWord);
    globals.words_per_edge = DivCeiling((globals.kmer_k + 1) * kBitsPerEdgeChar + kBitsPerMulti_t, kBitsPerEdgeWord);

    if (cx1_t::kCX1Verbose >= 2) {
        xlog("%d words per substring, %d words per edge\n", globals.words_per_substring, globals.words_per_edge);
    }

    // FILE *buc = fopen("bucket.txt", "w");
    // fprintf(buc, "avg %lld\n", globals.tot_bucket_size/ num_non_empty);
    // for (int i = 0; i < kNumBuckets; ++i) {
    //     fprintf(buc, "%d %lld\n", i, globals.cx1.bucket_sizes_[i]);
    // }
    // fclose(buc);

#ifdef USE_GPU
    globals.cx1.lv1_just_go_ = false;
    int64_t lv2_mem = globals.gpu_mem - 1073741824; // should reserver ~1G for GPU sorting
    globals.cx1.max_lv2_items_ = std::min(lv2_mem / cx1_t::kGPUBytePerItem, std::max(globals.max_bucket_size, kMinLv2BatchSizeGPU));

    if (globals.max_bucket_size > globals.cx1.max_lv2_items_) {
        xerr_and_exit("Bucket too large for GPU: contains %lld items. Please try CPU version.\n", globals.max_bucket_size);
        // TODO: auto switch to CPU version
    }

    // lv2 bytes: substring, permutation, readinfo
    int64_t lv2_bytes_per_item = globals.words_per_substring * sizeof(uint32_t) + sizeof(uint32_t) + sizeof(int64_t);
    lv2_bytes_per_item = lv2_bytes_per_item * 2; // double buffering
    // --- memory stuff ---
    int64_t mem_remained = globals.host_mem
                           - globals.mem_packed_reads
                           - globals.num_reads * sizeof(unsigned char) * 2 // first_in0 & last_out0
                           - kNumBuckets * sizeof(int64_t) * (globals.num_cpu_threads * 3 + 1)
                           - (kMaxMulti_t + 1) * (globals.num_output_threads + 1) * sizeof(int64_t);
    int64_t min_lv1_items = globals.tot_bucket_size / (kMaxLv1ScanTime - 0.5);
    int64_t min_lv2_items = std::max(globals.max_bucket_size, kMinLv2BatchSize);

    if (globals.mem_flag == 1) {
        // auto set memory
        globals.cx1.max_lv1_items_ = std::max(globals.cx1.max_lv2_items_, int64_t(globals.tot_bucket_size / (kDefaultLv1ScanTime - 0.5)));
        int64_t mem_needed = globals.cx1.max_lv1_items_ * cx1_t::kLv1BytePerItem + globals.cx1.max_lv2_items_ * lv2_bytes_per_item;

        if (mem_needed > mem_remained) {
            globals.cx1.adjust_mem(mem_remained, lv2_bytes_per_item, min_lv1_items, min_lv2_items);
        }
    }
    else if (globals.mem_flag == 0) {
        // min memory
        globals.cx1.max_lv1_items_ = std::max(globals.cx1.max_lv2_items_, int64_t(globals.tot_bucket_size / (kMaxLv1ScanTime - 0.5)));
        int64_t mem_needed = globals.cx1.max_lv1_items_ * cx1_t::kLv1BytePerItem + globals.cx1.max_lv2_items_ * lv2_bytes_per_item;

        if (mem_needed > mem_remained) {
            globals.cx1.adjust_mem(mem_remained, lv2_bytes_per_item, min_lv1_items, min_lv2_items);
        }
        else {
            globals.cx1.adjust_mem(mem_needed, lv2_bytes_per_item, min_lv1_items, min_lv2_items);
        }
    }
    else {
        // use all
        globals.cx1.adjust_mem(mem_remained, lv2_bytes_per_item, min_lv1_items, min_lv2_items);
    }

    // --- alloc memory ---
    globals.lv1_items = (int *) MallocAndCheck(globals.cx1.max_lv1_items_ * sizeof(int), __FILE__, __LINE__);
    globals.lv2_substrings = (uint32_t *) MallocAndCheck(globals.cx1.max_lv2_items_ * globals.words_per_substring * sizeof(uint32_t), __FILE__, __LINE__);
    globals.lv2_substrings_db = (uint32_t *) MallocAndCheck(globals.cx1.max_lv2_items_ * globals.words_per_substring * sizeof(uint32_t), __FILE__, __LINE__);
    globals.permutation = (uint32_t *) MallocAndCheck(globals.cx1.max_lv2_items_ * sizeof(uint32_t), __FILE__, __LINE__);
    globals.permutation_db = (uint32_t *) MallocAndCheck(globals.cx1.max_lv2_items_ * sizeof(uint32_t), __FILE__, __LINE__);
    globals.lv2_read_info = (int64_t *) MallocAndCheck(globals.cx1.max_lv2_items_ * sizeof(int64_t), __FILE__, __LINE__);
    globals.lv2_read_info_db = (int64_t *) MallocAndCheck(globals.cx1.max_lv2_items_ * sizeof(int64_t), __FILE__, __LINE__);
    alloc_gpu_buffers(globals.gpu_key_buffer1, globals.gpu_key_buffer2, globals.gpu_value_buffer1, globals.gpu_value_buffer2, (size_t)globals.cx1.max_lv2_items_);
#else

    globals.cx1.lv1_just_go_ = true;
    globals.num_output_threads = globals.num_cpu_threads;

    num_non_empty = std::max(1, num_non_empty);

    for (int i = 0; i < kNumBuckets; ++i) {
        if (globals.cx1.bucket_sizes_[i] > 2 * globals.tot_bucket_size / num_non_empty) {
            // xlog("Bucket %d size = %lld > %lld = 2 * avg\n", i, (long long)globals.cx1.bucket_sizes_[i], (long long)2 * globals.tot_bucket_size / num_non_empty);
        }
    }

    globals.max_sorting_items = std::max(3 * globals.tot_bucket_size / num_non_empty * globals.num_cpu_threads, globals.max_bucket_size * 2);

    // lv2 bytes: substring, permutation, readinfo
    int64_t lv2_bytes_per_item = globals.words_per_substring * sizeof(uint32_t) + sizeof(uint32_t) + sizeof(int64_t);
    lv2_bytes_per_item += sizeof(uint32_t); // cpu sorting space

    int64_t mem_remained = globals.host_mem
                           - globals.mem_packed_reads
                           - globals.num_cpu_threads * 65536 * sizeof(uint64_t) // radix sort buckets
                           - globals.num_reads * sizeof(unsigned char) * 2 // first_in0 & last_out0
                           - kNumBuckets * sizeof(int64_t) * (globals.num_cpu_threads * 3 + 1)
                           - (kMaxMulti_t + 1) * (globals.num_output_threads + 1) * sizeof(int64_t);
    int64_t min_lv1_items = globals.tot_bucket_size / (kMaxLv1ScanTime - 0.5);

    if (globals.mem_flag == 1) {
        // auto set memory
        globals.cx1.max_lv1_items_ = int64_t(globals.tot_bucket_size / (kDefaultLv1ScanTime - 0.5));
        int64_t mem_needed = globals.cx1.max_lv1_items_ * cx1_t::kLv1BytePerItem;

        if (mem_needed > mem_remained) {
            globals.cx1.adjust_mem_just_go(mem_remained, lv2_bytes_per_item, min_lv1_items, globals.max_bucket_size,
                                           globals.max_sorting_items, globals.cx1.max_lv1_items_, globals.max_sorting_items);
        }

    }
    else if (globals.mem_flag == 0) {
        // min memory
        globals.cx1.max_lv1_items_ = int64_t(globals.tot_bucket_size / (kMaxLv1ScanTime - 0.5));
        int64_t mem_needed = globals.cx1.max_lv1_items_ * cx1_t::kLv1BytePerItem;

        if (mem_needed > mem_remained) {
            globals.cx1.adjust_mem_just_go(mem_remained, lv2_bytes_per_item, min_lv1_items, globals.max_bucket_size,
                                           globals.max_sorting_items, globals.cx1.max_lv1_items_, globals.max_sorting_items);
        }
        else {
            globals.cx1.adjust_mem_just_go(mem_needed, lv2_bytes_per_item, min_lv1_items, globals.max_bucket_size,
                                           globals.max_sorting_items, globals.cx1.max_lv1_items_, globals.max_sorting_items);
        }

    }
    else {
        // use all
        globals.cx1.adjust_mem_just_go(mem_remained, lv2_bytes_per_item, min_lv1_items, globals.max_bucket_size,
                                       globals.max_sorting_items, globals.cx1.max_lv1_items_, globals.max_sorting_items);
    }

    if (globals.cx1.max_lv1_items_ < min_lv1_items) {
        xerr_and_exit("No enough memory to process.");
    }

    globals.cx1.max_mem_remain_ = globals.cx1.max_lv1_items_ * cx1_t::kLv1BytePerItem + globals.max_sorting_items * lv2_bytes_per_item;
    globals.cx1.bytes_per_sorting_item_ = lv2_bytes_per_item;
    globals.lv1_items = (int32_t *)MallocAndCheck(globals.cx1.max_mem_remain_ + globals.num_cpu_threads * sizeof(uint64_t) * 65536, __FILE__, __LINE__);

#endif

    if (cx1_t::kCX1Verbose >= 2) {
        xlog("Memory for reads: %lld\n", globals.mem_packed_reads);
        xlog("max # lv.1 items = %lld\n", globals.cx1.max_lv1_items_);
#ifdef USE_GPU
        xlog("max # lv.2 items = %lld\n", globals.cx1.max_lv2_items_);
#endif
    }

    // --- malloc read first_in / last_out ---
#ifdef LONG_READS
    globals.first_0_out = (uint32_t *) MallocAndCheck(globals.num_reads * sizeof(uint32_t), __FILE__, __LINE__);
    globals.last_0_in = (uint32_t *) MallocAndCheck(globals.num_reads * sizeof(uint32_t), __FILE__, __LINE__);
#else
    globals.first_0_out = (uint8_t *) MallocAndCheck(globals.num_reads * sizeof(uint8_t), __FILE__, __LINE__);
    globals.last_0_in = (uint8_t *) MallocAndCheck(globals.num_reads * sizeof(uint8_t), __FILE__, __LINE__);
#endif
    memset(globals.first_0_out, 0xFF, globals.num_reads * sizeof(globals.first_0_out[0]));
    memset(globals.last_0_in, 0xFF, globals.num_reads * sizeof(globals.last_0_in[0]));

    // --- initialize stat ---
    globals.edge_counting = (int64_t *) MallocAndCheck((kMaxMulti_t + 1) * sizeof(int64_t), __FILE__, __LINE__);
    globals.thread_edge_counting = (int64_t *) MallocAndCheck((kMaxMulti_t + 1) * globals.num_output_threads * sizeof(int64_t), __FILE__, __LINE__);
    memset(globals.edge_counting, 0, (kMaxMulti_t + 1) * sizeof(int64_t));

    // --- initialize lock ---
    pthread_mutex_init(&globals.lv1_items_scanning_lock, NULL); // init lock

    // --- initialize writer ---
    globals.edge_writer.set_file_prefix(globals.output_prefix);
    globals.edge_writer.set_num_threads(globals.num_output_threads);
    globals.edge_writer.set_kmer_size(globals.kmer_k);
    globals.edge_writer.set_num_buckets(kNumBuckets);
    globals.edge_writer.init_files();
}

void *lv1_fill_offset(void *_data) {
    readpartition_data_t &rp = *((readpartition_data_t *) _data);
    count_global_t &globals = *(rp.globals);
    int64_t *prev_full_offsets = (int64_t *)MallocAndCheck(kNumBuckets * sizeof(int64_t), __FILE__, __LINE__); // temporary array for computing differentials
    assert(prev_full_offsets != NULL);

    for (int b = globals.cx1.lv1_start_bucket_; b < globals.cx1.lv1_end_bucket_; ++b)
        prev_full_offsets[b] = rp.rp_lv1_differential_base;

    // this loop is VERY similar to that in PreprocessScanToFillBucketSizesThread
    GenericKmer edge, rev_edge; // (k+1)-mer and its rc
    int key;

    for (int64_t read_id = rp.rp_start_id; read_id < rp.rp_end_id; ++read_id) {
        int read_length = globals.package.length(read_id);

        if (read_length < globals.kmer_k + 1) {
            continue;
        }

        uint32_t *read_p = &globals.package.packed_seq[globals.package.get_start_index(read_id) / 16];
        edge.init(read_p, globals.package.get_start_index(read_id) % 16, globals.kmer_k + 1);
        rev_edge = edge;
        rev_edge.ReverseComplement(globals.kmer_k + 1);

        // ===== this is a macro to save some copy&paste ================
#define CHECK_AND_SAVE_OFFSET(offset, strand)                                                                   \
    do {                                                                                                        \
        if (globals.cx1.cur_lv1_buckets_[key]) {                                                                \
            int key_ = globals.cx1.bucket_rank_[key];                                                           \
            int64_t full_offset = EncodeOffset(read_id, offset, strand, globals.package);                       \
            int64_t differential = full_offset - prev_full_offsets[key_];                                       \
            if (differential > cx1_t::kDifferentialLimit) {                                                     \
                pthread_mutex_lock(&globals.lv1_items_scanning_lock);                                           \
                globals.lv1_items[rp.rp_bucket_offsets[key_]++] = -globals.cx1.lv1_items_special_.size() - 1;   \
                globals.cx1.lv1_items_special_.push_back(full_offset);                                          \
                pthread_mutex_unlock(&globals.lv1_items_scanning_lock);                                         \
            } else {                                                                                            \
                assert((int) differential >= 0);                                                                \
                globals.lv1_items[rp.rp_bucket_offsets[key_]++] = (int) differential;                           \
            }                                                                                                   \
            prev_full_offsets[key_] = full_offset;                                                              \
        }                                                                                                       \
    } while (0)
        // ^^^^^ why is the macro surrounded by a do-while? please ask Google
        // =========== end macro ==========================

        // shift the key char by char
        int last_char_offset = globals.kmer_k;

        while (true) {
            if (rev_edge.cmp(edge, globals.kmer_k + 1) < 0) {
                key = rev_edge.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
                CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 1);
            }
            else {
                key = edge.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
                CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 0);
            }

            if (++last_char_offset >= read_length) {
                break;
            }
            else {
                int c = globals.package.get_base(read_id, last_char_offset);
                edge.ShiftAppend(c, globals.kmer_k + 1);
                rev_edge.ShiftPreappend(3 - c, globals.kmer_k + 1);
            }
        }
    }

#undef CHECK_AND_SAVE_OFFSET

    free(prev_full_offsets);
    return NULL;
}

inline void lv2_extract_substr_(uint32_t *substrings_p, int64_t *read_info_p, count_global_t &globals, int start_bucket, int end_bucket, int num_items) {
    int *lv1_p = globals.lv1_items + globals.cx1.rp_[0].rp_bucket_offsets[start_bucket];

    for (int b = start_bucket; b < end_bucket; ++b) {
        for (int t = 0; t < globals.num_cpu_threads; ++t) {
            int64_t full_offset = globals.cx1.rp_[t].rp_lv1_differential_base;
            int num = globals.cx1.rp_[t].rp_bucket_sizes[b];

            for (int i = 0; i < num; ++i) {
                if (*lv1_p >= 0) {
                    full_offset += *(lv1_p++);
                }
                else {
                    full_offset = globals.cx1.lv1_items_special_[-1 - * (lv1_p++)];
                }

                int64_t read_id = globals.package.get_id(full_offset >> 1);
                int strand = full_offset & 1;
                int offset = (full_offset >> 1) - globals.package.get_start_index(read_id);
                int num_chars_to_copy = globals.kmer_k + 1;

                int read_length = globals.package.length(read_id);
                int64_t which_word = globals.package.get_start_index(read_id) / 16;
                int start_offset = globals.package.get_start_index(read_id) % 16;
                int words_this_seq = DivCeiling(start_offset + read_length, 16);
                uint32_t *read_p = &globals.package.packed_seq[which_word];

                unsigned char prev, next;

                if (offset > 0) {
                    prev = globals.package.get_base(read_id, offset - 1);
                }
                else {
                    prev = kSentinelValue;
                }

                if (offset + globals.kmer_k + 1 < read_length) {
                    next = globals.package.get_base(read_id, offset + globals.kmer_k + 1);
                }
                else {
                    next = kSentinelValue;
                }

                if (strand == 0) {
                    CopySubstring(substrings_p, read_p, offset + start_offset, num_chars_to_copy,
                                  num_items, words_this_seq, globals.words_per_substring);
                    *read_info_p = (full_offset << 6) | (prev << 3) | next;
                }
                else {
                    CopySubstringRC(substrings_p, read_p, offset + start_offset, num_chars_to_copy,
                                    num_items, words_this_seq, globals.words_per_substring);
                    *read_info_p = (full_offset << 6) | ((next == kSentinelValue ? kSentinelValue : (3 - next)) << 3)
                                   | (prev == kSentinelValue ? kSentinelValue : (3 - prev));
                }

                substrings_p++;
                read_info_p++;
            }
        }
    }
}

void *lv2_extract_substr(void *_data) {
    bucketpartition_data_t &bp = *((bucketpartition_data_t *) _data);
    count_global_t &globals = *(bp.globals);
    uint32_t *substrings_p = globals.lv2_substrings +
                             (globals.cx1.rp_[0].rp_bucket_offsets[bp.bp_start_bucket] - globals.cx1.rp_[0].rp_bucket_offsets[globals.cx1.lv2_start_bucket_]);
    int64_t *read_info_p = globals.lv2_read_info +
                           (globals.cx1.rp_[0].rp_bucket_offsets[bp.bp_start_bucket] - globals.cx1.rp_[0].rp_bucket_offsets[globals.cx1.lv2_start_bucket_]);

    lv2_extract_substr_(substrings_p, read_info_p, globals, bp.bp_start_bucket, bp.bp_end_bucket, globals.cx1.lv2_num_items_);
    return NULL;
}

void lv2_sort(count_global_t &globals) {
    xtimer_t local_timer;
#ifndef USE_GPU

    if (cx1_t::kCX1Verbose >= 4) {
        local_timer.reset();
        local_timer.start();
    }

    omp_set_num_threads(globals.num_cpu_threads - globals.num_output_threads);
    lv2_cpu_sort(globals.lv2_substrings, globals.permutation, globals.cpu_sort_space, globals.words_per_substring, globals.cx1.lv2_num_items_);
    omp_set_num_threads(globals.num_cpu_threads);
    local_timer.stop();

    if (cx1_t::kCX1Verbose >= 4) {
        xlog("Sorting substrings with CPU...done. Time elapsed: %.4lf\n", local_timer.elapsed());
    }

#else

    if (cx1_t::kCX1Verbose >= 4) {
        local_timer.reset();
        local_timer.start();
    }

    lv2_gpu_sort(globals.lv2_substrings, globals.permutation, globals.words_per_substring, globals.cx1.lv2_num_items_,
                 globals.gpu_key_buffer1, globals.gpu_key_buffer2, globals.gpu_value_buffer1, globals.gpu_value_buffer2);

    if (cx1_t::kCX1Verbose >= 4) {
        local_timer.stop();
        xlog("Sorting substrings with GPU...done. Time elapsed: %.4lf\n", local_timer.elapsed());
    }

#endif
}

void lv2_pre_output_partition(count_global_t &globals) {
    // swap double buffers
    globals.lv2_num_items_db = globals.cx1.lv2_num_items_;
    std::swap(globals.lv2_substrings_db, globals.lv2_substrings);
    std::swap(globals.permutation_db, globals.permutation);
    std::swap(globals.lv2_read_info_db, globals.lv2_read_info);

    // distribute partition
    int64_t items_per_thread = globals.lv2_num_items_db / globals.num_output_threads;
    int64_t acc = 0, t = 0;
    globals.cx1.op_[t].op_start_index = 0;

    for (int b = globals.cx1.lv2_start_bucket_; b < globals.cx1.lv2_end_bucket_; ++b) {
        if (globals.cx1.bucket_sizes_[b] + acc > (t + 1) * items_per_thread) {
            globals.cx1.op_[t].op_end_index = acc;
            ++t;

            globals.cx1.op_[t].op_start_index = acc;

            if (t == globals.num_output_threads - 1) {
                break;
            }
        }

        acc += globals.cx1.bucket_sizes_[b];
    }

    globals.cx1.op_[t].op_end_index = globals.lv2_num_items_db;

    for (++t; t < globals.num_output_threads; ++t) {
        globals.cx1.op_[t].op_end_index = globals.lv2_num_items_db;
        globals.cx1.op_[t].op_start_index = globals.lv2_num_items_db;
    }
}

void lv2_output_(int64_t start_index, int64_t end_index, int thread_id, count_global_t &globals,
                 uint32_t *substrings, uint32_t *permutation, int64_t *read_infos, int num_items) {
    uint32_t packed_edge[32];
    int count_prev[5], count_next[5];
    int64_t *thread_edge_counting = globals.thread_edge_counting + thread_id * (kMaxMulti_t + 1);

    int from_;
    int to_;

    for (int i = start_index; i < end_index; i = to_) {
        from_ = i;
        to_ = i + 1;
        uint32_t *first_item = substrings + permutation[i];

        while (to_ < end_index) {
            if (IsDifferentEdges(first_item,
                                 substrings + permutation[to_],
                                 globals.words_per_substring, num_items)) {
                break;
            }

            ++to_;
        }

        int count = to_ - from_;

        // update read's first and last

        memset(count_prev, 0, sizeof(int) * 4);
        memset(count_next, 0, sizeof(int) * 4);
        bool has_in = false;
        bool has_out = false;

        for (int j = from_; j < to_; ++j) {
            int prev_and_next = read_infos[permutation[j]] & ((1 << 6) - 1);
            count_prev[prev_and_next >> 3]++;
            count_next[prev_and_next & 7]++;
        }

        for (int j = 0; j < 4; ++j) {
            if (count_prev[j] >= globals.kmer_freq_threshold) {
                has_in = true;
            }

            if (count_next[j] >= globals.kmer_freq_threshold) {
                has_out = true;
            }
        }

        if (!has_in && count >= globals.kmer_freq_threshold) {
            for (int j = from_; j < to_; ++j) {
                int64_t read_info = read_infos[permutation[j]] >> 6;
                int64_t read_id = globals.package.get_id(read_info >> 1);
                int strand = read_info & 1;
                uint32_t offset = (read_info >> 1) - globals.package.get_start_index(read_id);

                if (strand == 0) {
                    // update last
                    while (true) {
                        auto old_value = globals.last_0_in[read_id];

                        if (old_value != kSentinelOffset && old_value >= offset) {
                            break;
                        }

                        if (__sync_bool_compare_and_swap(globals.last_0_in + read_id, old_value, offset)) {
                            break;
                        }
                    }
                }
                else {
                    // update first
                    offset++;

                    while (true) {
                        auto old_value = globals.first_0_out[read_id];

                        if (old_value <= offset) {
                            break;
                        }

                        if (__sync_bool_compare_and_swap(globals.first_0_out + read_id, old_value, offset)) {
                            break;
                        }
                    }
                }
            }
        }

        if (!has_out && count >= globals.kmer_freq_threshold) {
            for (int j = from_; j < to_; ++j) {
                int64_t read_info = read_infos[permutation[j]] >> 6;
                int64_t read_id = globals.package.get_id(read_info >> 1);
                int strand = read_info & 1;
                uint32_t offset = (read_info >> 1) - globals.package.get_start_index(read_id);

                if (strand == 0) {
                    // update first
                    offset++;

                    while (true) {
                        auto old_value = globals.first_0_out[read_id];

                        if (old_value <= offset) {
                            break;
                        }

                        if (__sync_bool_compare_and_swap(globals.first_0_out + read_id, old_value, offset)) {
                            break;
                        }
                    }
                }
                else {
                    // update last
                    while (true) {
                        auto old_value = globals.last_0_in[read_id];

                        if (old_value != kSentinelOffset && old_value >= offset) {
                            break;
                        }

                        if (__sync_bool_compare_and_swap(globals.last_0_in + read_id, old_value, offset)) {
                            break;
                        }
                    }
                }
            }
        }

        ++thread_edge_counting[std::min(count, kMaxMulti_t)];

        if (count >= globals.kmer_freq_threshold) {
            PackEdge(packed_edge, first_item, count, globals, num_items);
            globals.edge_writer.write(packed_edge, packed_edge[0] >> (32 - 2 * kBucketPrefixLength), thread_id);
        }
    }
}

void *lv2_output(void *_op) {
    xtimer_t local_timer;

    if (cx1_t::kCX1Verbose >= 4) {
        local_timer.start();
        local_timer.reset();
    }

    outputpartition_data_t *op = (outputpartition_data_t *) _op;
    count_global_t &globals = *(op->globals);
    int64_t op_start_index = op->op_start_index;
    int64_t op_end_index = op->op_end_index;
    int thread_id = op->op_id;

    lv2_output_(op_start_index, op_end_index, thread_id, globals, globals.lv2_substrings_db, globals.permutation_db, globals.lv2_read_info_db, globals.lv2_num_items_db);

    if (cx1_t::kCX1Verbose >= 4) {
        local_timer.stop();
        xlog("Counting time elapsed: %.4lfs\n", local_timer.elapsed());
    }

    return NULL;
}

struct kt_sort_t {
    count_global_t *globals;
    std::vector<int64_t> thread_offset;
};

void kt_sort(void *_g, long i, int tid) {
    kt_sort_t *kg = (kt_sort_t *)_g;
    int b = kg->globals->cx1.lv1_start_bucket_ + i;

    if (kg->globals->cx1.bucket_sizes_[b] == 0) {
        return;
    }

    if (tid + 1 < kg->globals->num_cpu_threads) {
        kg->thread_offset[tid + 1] = kg->thread_offset[tid] + kg->globals->cx1.bucket_sizes_[b];
    }

    size_t offset = kg->globals->cx1.lv1_num_items_ * sizeof(int32_t) +
                    kg->thread_offset[tid] * kg->globals->cx1.bytes_per_sorting_item_ +
                    tid * sizeof(uint64_t) * 65536;

    uint32_t *substr_ptr = (uint32_t *) ((char *)kg->globals->lv1_items + offset);
    uint64_t *bucket = (uint64_t *)(substr_ptr + kg->globals->cx1.bucket_sizes_[b] * kg->globals->words_per_substring);
    uint32_t *permutation_ptr = (uint32_t *)(bucket + 65536);
    uint32_t *cpu_sort_space_ptr = permutation_ptr + kg->globals->cx1.bucket_sizes_[b];
    int64_t *readinfo_ptr = (int64_t *) (cpu_sort_space_ptr + kg->globals->cx1.bucket_sizes_[b]);

    lv2_extract_substr_(substr_ptr, readinfo_ptr, *(kg->globals), b, b + 1, kg->globals->cx1.bucket_sizes_[b]);
    lv2_cpu_radix_sort_st(substr_ptr, permutation_ptr, cpu_sort_space_ptr, bucket, kg->globals->words_per_substring, kg->globals->cx1.bucket_sizes_[b]);
    lv2_output_(0, kg->globals->cx1.bucket_sizes_[b], tid, *(kg->globals), substr_ptr, permutation_ptr, readinfo_ptr, kg->globals->cx1.bucket_sizes_[b]);
}

void lv1_direct_sort_and_count(count_global_t &globals) {
    kt_sort_t kg;
    kg.globals = &globals;
    kg.thread_offset.clear();
    int64_t acc_size = 0;

    for (int i = 0, b = globals.cx1.lv1_start_bucket_; i < globals.num_cpu_threads && b < globals.cx1.lv1_end_bucket_; ++i, ++b) {
        kg.thread_offset.push_back(acc_size);
        acc_size += globals.cx1.bucket_sizes_[b];
    }

    kt_dfor(globals.num_cpu_threads, kt_sort, &kg, globals.cx1.lv1_end_bucket_ - globals.cx1.lv1_start_bucket_);
}

void lv2_post_output(count_global_t &globals) {
}

void post_proc(count_global_t &globals) {
    for (int t = 0; t < globals.num_output_threads; ++t) {
        for (int i = 1; i <= kMaxMulti_t; ++i) {
            globals.edge_counting[i] += globals.thread_edge_counting[t * (kMaxMulti_t + 1) + i];
        }
    }

    // --- output reads for mercy ---
    int64_t num_candidate_reads = 0;
    int64_t num_has_tips = 0;
    FILE *candidate_file = OpenFileAndCheck((globals.output_prefix + ".cand").c_str(), "wb");
    SequenceManager seq_manager(&globals.package);

    for (int64_t i = 0; i < globals.num_reads; ++i) {
        auto first = globals.first_0_out[i];
        auto last = globals.last_0_in[i];

        if (first != kSentinelOffset && last != kSentinelOffset) {
            ++num_has_tips;

            if (last > first) {
                ++num_candidate_reads;
                seq_manager.WriteBinarySequences(candidate_file, false, i, i);
            }
        }
    }

    fclose(candidate_file);

    if (cx1_t::kCX1Verbose >= 2) {
        xlog("Total number of candidate reads: %lld(%lld)\n", num_candidate_reads, num_has_tips);
    }

    // --- stat ---
    int64_t num_solid_edges = 0;

    for (int i = globals.kmer_freq_threshold; i <= kMaxMulti_t; ++i) {
        num_solid_edges += globals.edge_counting[i];
    }

    if (cx1_t::kCX1Verbose >= 2) {
        xlog("Total number of solid edges: %llu\n", num_solid_edges);
    }

    FILE *counting_file = OpenFileAndCheck((globals.output_prefix + ".counting").c_str(), "w");

    for (int64_t i = 1, acc = 0; i <= kMaxMulti_t; ++i) {
        acc += globals.edge_counting[i];
        fprintf(counting_file, "%lld %lld\n", (long long)i, (long long)acc);
    }

    fclose(counting_file);

    // --- cleaning ---
    pthread_mutex_destroy(&globals.lv1_items_scanning_lock);
    free(globals.lv1_items);
#ifdef USE_GPU
    free(globals.lv2_substrings);
    free(globals.permutation);
    free(globals.permutation_db);
    free(globals.lv2_substrings_db);
    free(globals.lv2_read_info);
    free(globals.lv2_read_info_db);
    free_gpu_buffers(globals.gpu_key_buffer1, globals.gpu_key_buffer2, globals.gpu_value_buffer1, globals.gpu_value_buffer2);
#endif
    free(globals.first_0_out);
    free(globals.last_0_in);
    free(globals.edge_counting);
    free(globals.thread_edge_counting);
    globals.edge_writer.destroy();
}

} // namespace::cx1_kmer_count