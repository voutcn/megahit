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


#include "cx1_kmer_count.h"

#include <string.h>
#include <algorithm>
#include <zlib.h>
#include <omp.h>

#include "sdbg_builder_writers.h"
#include "mem_file_checker-inl.h"
#include "kseq.h"
#include "utils.h"
#include "kmer.h"
#include "packed_reads.h"
#include "read_lib_functions-inl.h"

#ifndef USE_GPU
#include "lv2_cpu_sort.h"
#else
#include "lv2_gpu_functions.h"
#endif

namespace cx1_kmer_count {

// helpers
typedef CX1<count_global_t, kNumBuckets> cx1_t;
typedef CX1<count_global_t, kNumBuckets>::readpartition_data_t readpartition_data_t;
typedef CX1<count_global_t, kNumBuckets>::bucketpartition_data_t bucketpartition_data_t;
typedef CX1<count_global_t, kNumBuckets>::outputpartition_data_t outputpartition_data_t;

/**
 * @brief encode read_id and its offset in one int64_t
 */
inline int64_t EncodeOffset(int64_t read_id, int offset, int strand, int length_num_bits) {
    return (read_id << (length_num_bits + 1)) | (offset << 1) | strand;
}

inline bool IsDifferentEdges(uint32_t *item1, uint32_t* item2, int num_words, int spacing) {
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
inline void PackEdge(uint32_t *dest, uint32_t *item, int counting, struct count_global_t &globals) {
    for (int i = 0; i < globals.words_per_edge && i < globals.words_per_substring; ++i) {
        dest[i] = *(item + (int64_t)i * globals.lv2_num_items_db);
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

// function pass to CX1

int64_t encode_lv1_diff_base(int64_t read_id, count_global_t &g) {
    return EncodeOffset(read_id, 0, 0, g.offset_num_bits);
}

void read_input_prepare(count_global_t &globals) { // num_items_, num_cpu_threads_ and num_output_threads_ must be set here
    bool is_reverse = true;
    ReadMultipleLibs(globals.read_lib_file, globals.package, globals.lib_info, is_reverse);
    globals.max_read_length = globals.package.max_read_len();
    globals.num_reads = globals.package.size();

    xlog("%ld reads, %d max read length\n", globals.num_reads, globals.max_read_length);
    xlog("Read libs:\n")
    for (unsigned i = 0; i < globals.lib_info.size(); ++i) {
        xlog("");
        xlog_ext("%s: %s, %lld reads\n", globals.lib_info[i].is_pe ? "Paired-end" : "Single-end",
                                         globals.lib_info[i].metadata.c_str(), globals.lib_info[i].to - globals.lib_info[i].from + 1);
    }

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

void* lv0_calc_bucket_size(void* _data) {
    readpartition_data_t &rp = *((readpartition_data_t*) _data);
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
            } else {
                bucket_sizes[edge.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
            }

            if (++last_char_offset >= read_length) {
                break;
            } else {
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
    for (int i = 0; i < kNumBuckets; ++i) {
        globals.tot_bucket_size += globals.cx1.bucket_sizes_[i];
    }

#ifndef USE_GPU
    globals.cx1.max_lv2_items_ = std::max(globals.max_bucket_size, kMinLv2BatchSize);
#else
    int64_t lv2_mem = globals.gpu_mem - 1073741824; // should reserver ~1G for GPU sorting
    globals.cx1.max_lv2_items_ = std::min(lv2_mem / cx1_t::kGPUBytePerItem, std::max(globals.max_bucket_size, kMinLv2BatchSizeGPU));
    if (globals.max_bucket_size > globals.cx1.max_lv2_items_) {
        xerr_and_exit("Bucket too large for GPU: contains %lld items. Please try CPU version.\n", globals.max_bucket_size);
        // TODO: auto switch to CPU version
    }
#endif
    globals.words_per_substring = DivCeiling((globals.kmer_k + 1) * kBitsPerEdgeChar, kBitsPerEdgeWord);
    globals.words_per_edge = DivCeiling((globals.kmer_k + 1) * kBitsPerEdgeChar + kBitsPerMulti_t, kBitsPerEdgeWord);
    // lv2 bytes: substring, permutation, readinfo
    int64_t lv2_bytes_per_item = (globals.words_per_substring) * sizeof(uint32_t) + sizeof(uint32_t) + sizeof(int64_t);
    lv2_bytes_per_item = lv2_bytes_per_item * 2; // double buffering
#ifndef USE_GPU
    lv2_bytes_per_item += sizeof(uint64_t) * 2; // CPU memory is used to simulate GPU
#endif

    if (cx1_t::kCX1Verbose >= 2) {
        xlog("%d words per substring, %d words per edge\n", globals.words_per_substring, globals.words_per_edge);
    }
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
    } else if (globals.mem_flag == 0) {
        // min memory
        globals.cx1.max_lv1_items_ = std::max(globals.cx1.max_lv2_items_, int64_t(globals.tot_bucket_size / (kMaxLv1ScanTime - 0.5)));
        int64_t mem_needed = globals.cx1.max_lv1_items_ * cx1_t::kLv1BytePerItem + globals.cx1.max_lv2_items_ * lv2_bytes_per_item;
        if (mem_needed > mem_remained) {
            globals.cx1.adjust_mem(mem_remained, lv2_bytes_per_item, min_lv1_items, min_lv2_items);
        } else {
            globals.cx1.adjust_mem(mem_needed, lv2_bytes_per_item, min_lv1_items, min_lv2_items);
        }
    } else {
        // use all
        globals.cx1.adjust_mem(mem_remained, lv2_bytes_per_item, min_lv1_items, min_lv2_items);
    }

    if (cx1_t::kCX1Verbose >= 2) {
        xlog("Memory for reads: %lld\n", globals.mem_packed_reads);
        xlog("max # lv.1 items = %lld\n", globals.cx1.max_lv1_items_);
        xlog("max # lv.2 items = %lld\n", globals.cx1.max_lv2_items_);
    }

    // --- alloc memory ---
    globals.lv1_items = (int*) MallocAndCheck(globals.cx1.max_lv1_items_ * sizeof(int), __FILE__, __LINE__);
    globals.lv2_substrings = (uint32_t*) MallocAndCheck(globals.cx1.max_lv2_items_ * globals.words_per_substring * sizeof(uint32_t), __FILE__, __LINE__);
    globals.lv2_substrings_db = (uint32_t*) MallocAndCheck(globals.cx1.max_lv2_items_ * globals.words_per_substring * sizeof(uint32_t), __FILE__, __LINE__);
    globals.permutation = (uint32_t *) MallocAndCheck(globals.cx1.max_lv2_items_ * sizeof(uint32_t), __FILE__, __LINE__);
    globals.permutation_db = (uint32_t *) MallocAndCheck(globals.cx1.max_lv2_items_ * sizeof(uint32_t), __FILE__, __LINE__);
    globals.lv2_read_info = (int64_t *) MallocAndCheck(globals.cx1.max_lv2_items_ * sizeof(int64_t), __FILE__, __LINE__);
    globals.lv2_read_info_db = (int64_t *) MallocAndCheck(globals.cx1.max_lv2_items_ * sizeof(int64_t), __FILE__, __LINE__);
#ifdef LONG_READS
    globals.first_0_out = (uint16_t*) MallocAndCheck(globals.num_reads * sizeof(uint16_t), __FILE__, __LINE__);
    globals.last_0_in = (uint16_t*) MallocAndCheck(globals.num_reads * sizeof(uint16_t), __FILE__, __LINE__);
#else
    globals.first_0_out = (unsigned char*) MallocAndCheck(globals.num_reads * sizeof(unsigned char), __FILE__, __LINE__);
    globals.last_0_in = (unsigned char*) MallocAndCheck(globals.num_reads * sizeof(unsigned char), __FILE__, __LINE__);
#endif
#ifndef USE_GPU
    globals.cpu_sort_space = (uint64_t*) MallocAndCheck(sizeof(uint64_t) * globals.cx1.max_lv2_items_, __FILE__, __LINE__);
#else
    alloc_gpu_buffers(globals.gpu_key_buffer1, globals.gpu_key_buffer2, globals.gpu_value_buffer1, globals.gpu_value_buffer2, (size_t)globals.cx1.max_lv2_items_);
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
    globals.word_writer = (WordWriter*) MallocAndCheck(globals.num_output_threads * sizeof(WordWriter), __FILE__, __LINE__);
    for (int t = 0; t < globals.num_output_threads; ++t) {
        globals.word_writer[t].init(FormatString("%s.edges.%d", globals.output_prefix.c_str(), t));
    }

    // --- write the edge file header ---
    globals.word_writer[0].output(globals.kmer_k);
    globals.word_writer[0].output(globals.words_per_edge);
}

void* lv1_fill_offset(void* _data) {
    readpartition_data_t &rp = *((readpartition_data_t*) _data);
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
        assert(offset + globals.kmer_k < read_length);                                                          \
        if (((key - globals.cx1.lv1_start_bucket_) ^ (key - globals.cx1.lv1_end_bucket_)) & kSignBitMask) {     \
            int64_t full_offset = EncodeOffset(read_id, offset, strand, globals.offset_num_bits);               \
            int64_t differential = full_offset - prev_full_offsets[key];                                        \
            if (differential > cx1_t::kDifferentialLimit) {                                                     \
                pthread_mutex_lock(&globals.lv1_items_scanning_lock);                                           \
                globals.lv1_items[rp.rp_bucket_offsets[key]++] = -globals.cx1.lv1_items_special_.size() - 1;    \
                globals.cx1.lv1_items_special_.push_back(full_offset);                                          \
                pthread_mutex_unlock(&globals.lv1_items_scanning_lock);                                         \
            } else {                                                                                            \
                assert((int) differential >= 0);                                                                \
                globals.lv1_items[rp.rp_bucket_offsets[key]++] = (int) differential;                            \
            }                                                                                                   \
            prev_full_offsets[key] = full_offset;                                                               \
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
            } else {
                key = edge.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
                CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 0);
            }

            if (++last_char_offset >= read_length) {
                break;
            } else {
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

void* lv2_extract_substr(void* _data) {
    bucketpartition_data_t &bp = *((bucketpartition_data_t*) _data);
    count_global_t &globals = *(bp.globals);
    int *lv1_p = globals.lv1_items + globals.cx1.rp_[0].rp_bucket_offsets[bp.bp_start_bucket];
    int64_t offset_mask = (1 << globals.offset_num_bits) - 1; // 0000....00011..11
    uint32_t *substrings_p = globals.lv2_substrings +
                             (globals.cx1.rp_[0].rp_bucket_offsets[bp.bp_start_bucket] - globals.cx1.rp_[0].rp_bucket_offsets[globals.cx1.lv2_start_bucket_]);
    int64_t *read_info_p = globals.lv2_read_info +
                           (globals.cx1.rp_[0].rp_bucket_offsets[bp.bp_start_bucket] - globals.cx1.rp_[0].rp_bucket_offsets[globals.cx1.lv2_start_bucket_]);

    for (int b = bp.bp_start_bucket; b < bp.bp_end_bucket; ++b) {
        for (int t = 0; t < globals.num_cpu_threads; ++t) {
            int64_t full_offset = globals.cx1.rp_[t].rp_lv1_differential_base;
            int num = globals.cx1.rp_[t].rp_bucket_sizes[b];
            for (int i = 0; i < num; ++i) {
                if (*lv1_p >= 0) {
                    full_offset += *(lv1_p++);
                } else {
                    full_offset = globals.cx1.lv1_items_special_[-1 - *(lv1_p++)];
                }

                int64_t read_id = full_offset >> (globals.offset_num_bits + 1);
                int strand = full_offset & 1;
                int offset = (full_offset >> 1) & offset_mask;
                int num_chars_to_copy = globals.kmer_k + 1;

                int read_length = globals.package.length(read_id);
                int64_t which_word = globals.package.get_start_index(read_id) / 16;
                int start_offset = globals.package.get_start_index(read_id) % 16;
                int words_this_seq = DivCeiling(start_offset + read_length, 16);
                uint32_t *read_p = &globals.package.packed_seq[which_word];

                unsigned char prev, next;

                if (offset > 0) {
                    prev = globals.package.get_base(read_id, offset - 1);
                } else {
                    prev = kSentinelValue;
                }

                if (offset + globals.kmer_k + 1 < read_length) {
                    next = globals.package.get_base(read_id, offset + globals.kmer_k + 1);
                } else {
                    next = kSentinelValue;
                }

                if (strand == 0) {
                    CopySubstring(substrings_p, read_p, offset + start_offset, num_chars_to_copy,
                                  globals.cx1.lv2_num_items_, words_this_seq, globals.words_per_substring);
                    *read_info_p = (full_offset << 6) | (prev << 3) | next;
                } else {
                    CopySubstringRC(substrings_p, read_p, offset + start_offset, num_chars_to_copy,
                                    globals.cx1.lv2_num_items_, words_this_seq, globals.words_per_substring);
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
    int64_t last_end_index = 0;
    int64_t items_per_thread = globals.lv2_num_items_db / globals.num_output_threads;

    for (int t = 0; t < globals.num_output_threads - 1; ++t) {
        int64_t this_start_index = last_end_index;
        int64_t this_end_index = this_start_index + items_per_thread;
        if (this_end_index > globals.lv2_num_items_db) {
            this_end_index = globals.lv2_num_items_db;
        }
        if (this_end_index > 0) {
            while (this_end_index < globals.lv2_num_items_db) {
                uint32_t *prev_item = globals.lv2_substrings_db + (globals.permutation_db[this_end_index - 1]);
                uint32_t *item = globals.lv2_substrings_db + (globals.permutation_db[this_end_index]);
                if (IsDifferentEdges(prev_item, item, globals.words_per_substring, globals.lv2_num_items_db)) {
                    break;
                }
                ++this_end_index;
            }
        }
        globals.cx1.op_[t].op_start_index = this_start_index;
        globals.cx1.op_[t].op_end_index = this_end_index;
        last_end_index = this_end_index;
    }

    // last partition
    globals.cx1.op_[globals.num_output_threads - 1].op_start_index = last_end_index;
    globals.cx1.op_[globals.num_output_threads - 1].op_end_index = globals.lv2_num_items_db;

    memset(globals.thread_edge_counting, 0, sizeof(int64_t) * (kMaxMulti_t + 1) * globals.num_output_threads);
}

void* lv2_output(void* _op) {
    xtimer_t local_timer;
    if (cx1_t::kCX1Verbose >= 4) {
        local_timer.start();
        local_timer.reset();
    }

    outputpartition_data_t *op = (outputpartition_data_t*) _op;
    count_global_t &globals = *(op->globals);
    int64_t op_start_index = op->op_start_index;
    int64_t op_end_index = op->op_end_index;
    int thread_id = op->op_id;
    int start_idx;
    int end_idx;
    uint32_t packed_edge[32];
    int count_prev[5], count_next[5];
    int64_t offset_mask = (1 << globals.offset_num_bits) - 1;
    int64_t *thread_edge_counting = globals.thread_edge_counting + thread_id * (kMaxMulti_t + 1);

    for (int i = op_start_index; i < op_end_index; i = end_idx) {
        start_idx = i;
        end_idx = i + 1;
        uint32_t *first_item = globals.lv2_substrings_db + (globals.permutation_db[i]);
        while (end_idx < op_end_index) {
            if (IsDifferentEdges(first_item,
                                 globals.lv2_substrings_db + globals.permutation_db[end_idx],
                                 globals.words_per_substring, globals.lv2_num_items_db)) {
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
            int prev_and_next = globals.lv2_read_info_db[globals.permutation_db[j]] & ((1 << 6) - 1);
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
            for (int j = start_idx; j < end_idx; ++j) {
                int64_t read_info = globals.lv2_read_info_db[globals.permutation_db[j]] >> 6;
                int strand = read_info & 1;
                int offset = (read_info >> 1) & offset_mask;
                int64_t read_id = read_info >> (1 + globals.offset_num_bits);

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
                } else {
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
            for (int j = start_idx; j < end_idx; ++j) {
                int64_t read_info = globals.lv2_read_info_db[globals.permutation_db[j]] >> 6;
                int strand = read_info & 1;
                int offset = (read_info >> 1) & offset_mask;
                int64_t read_id = read_info >> (1 + globals.offset_num_bits);

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
                } else {
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
            PackEdge(packed_edge, first_item, count, globals);
            for (int x = 0; x < globals.words_per_edge; ++x) {
                globals.word_writer[thread_id].output(packed_edge[x]);
            }
        }
    }

    if (cx1_t::kCX1Verbose >= 4) {
        local_timer.stop();
        xlog("Counting time elapsed: %.4lfs\n", local_timer.elapsed());
    }
    return NULL;
}

void lv2_post_output(count_global_t &globals) {
    for (int t = 0; t < globals.num_output_threads; ++t) {
        for (int i = 1; i <= kMaxMulti_t; ++i) {
            globals.edge_counting[i] += globals.thread_edge_counting[t * (kMaxMulti_t + 1) + i];
        }
    }
}

void post_proc(count_global_t &globals) {
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

    FILE *counting_file = OpenFileAndCheck((std::string(globals.output_prefix)+".counting").c_str(), "w");
    for (int64_t i = 1, acc = 0; i <= kMaxMulti_t; ++i) {
        acc += globals.edge_counting[i];
        fprintf(counting_file, "%lld %lld\n", (long long)i, (long long)acc);
    }
    fclose(counting_file);

    if (cx1_t::kCX1Verbose >= 2) {
        xlog("Output reads to binary file...\n");
        bool is_reverse = true;
        WriteMultipleLibs(globals.package, globals.lib_info, globals.output_prefix + ".all_reads", is_reverse);
    }

    // --- cleaning ---
    pthread_mutex_destroy(&globals.lv1_items_scanning_lock);
    free(globals.lv1_items);
    free(globals.lv2_substrings);
    free(globals.permutation);
    free(globals.permutation_db);
    free(globals.lv2_substrings_db);
    free(globals.lv2_read_info);
    free(globals.lv2_read_info_db);
    free(globals.first_0_out);
    free(globals.last_0_in);
    free(globals.edge_counting);
    free(globals.thread_edge_counting);
    for (int t = 0; t < globals.num_output_threads; ++t) {
        globals.word_writer[t].destroy();
    }
    free(globals.word_writer);

#ifndef USE_GPU
    free(globals.cpu_sort_space);
#else
    free_gpu_buffers(globals.gpu_key_buffer1, globals.gpu_key_buffer2, globals.gpu_value_buffer1, globals.gpu_value_buffer2);
#endif
}

} // namespace::cx1_kmer_count