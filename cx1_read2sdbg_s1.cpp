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

#include "cx1_read2sdbg.h"

#include <string.h>
#include <algorithm>
#include <zlib.h>
#include <omp.h>

#include "mem_file_checker-inl.h"
#include "kseq.h"
#include "utils.h"
#include "kmer.h"
#include "packed_reads.h"
#include "sequence_package.h"
#include "read_lib_functions-inl.h"

#ifndef USE_GPU
#include "lv2_cpu_sort.h"
#else
#include "lv2_gpu_functions.h"
#endif
// helping functions

namespace cx1_read2sdbg {

namespace s1 {

// helpers
typedef CX1<read2sdbg_global_t, kNumBuckets> cx1_t;
typedef CX1<read2sdbg_global_t, kNumBuckets>::readpartition_data_t readpartition_data_t;
typedef CX1<read2sdbg_global_t, kNumBuckets>::bucketpartition_data_t bucketpartition_data_t;
typedef CX1<read2sdbg_global_t, kNumBuckets>::outputpartition_data_t outputpartition_data_t;

/**
 * @brief encode read_id and its offset in one int64_t
 */
inline int64_t EncodeOffset(int64_t read_id, int offset, int strand, int length_num_bits) {
    return (read_id << (length_num_bits + 1)) | (offset << 1) | strand;
}

// helper: see whether two lv2 items have the same (k-1)-mer
inline bool IsDiffKMinusOneMer(uint32_t *item1, uint32_t *item2, int64_t spacing, int kmer_k) {
    // mask extra bits
    int chars_in_last_word = (kmer_k - 1) % kCharsPerEdgeWord;
    int num_full_words = (kmer_k - 1) / kCharsPerEdgeWord;
    if (chars_in_last_word > 0) {
        uint32_t w1 = item1[num_full_words * spacing];
        uint32_t w2 = item2[num_full_words * spacing];
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

inline uint8_t ExtractHeadTail(uint32_t *item, int64_t spacing, int words_per_substring) {
    return *(item + spacing * (words_per_substring - 1)) & ((1 << 2 * kBWTCharNumBits) - 1);
}

inline uint8_t ExtractPrevNext(int i, read2sdbg_global_t &globals) {
    return globals.lv2_read_info_db[i] & ((1 << 2 * kBWTCharNumBits) - 1);
}

// cx1 core functions

int64_t s1_encode_lv1_diff_base(int64_t read_id, read2sdbg_global_t &g) {
    return EncodeOffset(read_id, 0, 0, g.offset_num_bits);
}

void s1_read_input_prepare(read2sdbg_global_t &globals) {
    bool is_reverse = true;
    ReadMultipleLibs(globals.read_lib_file, globals.package, globals.lib_info, is_reverse);
    globals.max_read_length = globals.package.max_read_len();
    globals.num_reads = globals.package.size();

    int bits_read_length = 1; // bit needed to store read_length
    while ((1 << bits_read_length) - 1 < globals.max_read_length) {
        ++bits_read_length;
    }
    globals.read_length_mask = (1 << bits_read_length) - 1;
    globals.offset_num_bits = bits_read_length;

    // --- allocate memory for is_solid bit_vector
    globals.num_k1_per_read = globals.max_read_length - globals.kmer_k;
    globals.is_solid.reset(globals.num_k1_per_read * globals.num_reads);
    globals.mem_packed_reads = DivCeiling(globals.num_k1_per_read * globals.num_reads, 8) + globals.package.size_in_byte();

    xlog("Total number of reads: %lld, bases: %lld\n", (long long)globals.num_reads, (long long)globals.package.base_size());

    // set cx1 param
    globals.cx1.num_cpu_threads_ = globals.num_cpu_threads;
    globals.cx1.num_output_threads_ = globals.num_output_threads;
    globals.cx1.num_items_ = globals.num_reads;
}

void* s1_lv0_calc_bucket_size(void* _data) {
    readpartition_data_t &rp = *((readpartition_data_t*) _data);
    read2sdbg_global_t &globals = *(rp.globals);
    int64_t *bucket_sizes = rp.rp_bucket_sizes;
    memset(bucket_sizes, 0, kNumBuckets * sizeof(int64_t));
    GenericKmer k_minus1_mer, rev_k_minus1_mer; // (k-1)-mer and its rc

    for (int64_t read_id = rp.rp_start_id; read_id < rp.rp_end_id; ++read_id) {
        int read_length = globals.package.length(read_id);
        if (read_length < globals.kmer_k + 1) {
            continue;
        }

        int64_t which_word = globals.package.get_start_index(read_id) / 16;
        int64_t offset = globals.package.get_start_index(read_id) % 16;
        uint32_t *read_p = &globals.package.packed_seq[which_word];

        k_minus1_mer.init(read_p, offset, globals.kmer_k - 1);
        rev_k_minus1_mer = k_minus1_mer;
        rev_k_minus1_mer.ReverseComplement(globals.kmer_k - 1);

        // the first one special handling
        bucket_sizes[k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
        bucket_sizes[rev_k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;

        int last_char_offset = globals.kmer_k - 1;
        int c = globals.package.get_base(read_id, last_char_offset);
        k_minus1_mer.ShiftAppend(c, globals.kmer_k - 1);
        rev_k_minus1_mer.ShiftPreappend(3 - c, globals.kmer_k - 1);

        while (last_char_offset < read_length - 1) {
            int cmp = k_minus1_mer.cmp(rev_k_minus1_mer, globals.kmer_k - 1);
            if (cmp > 0) {
                bucket_sizes[rev_k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
            } else {
                bucket_sizes[k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
            }

            int c = globals.package.get_base(read_id, ++last_char_offset);
            k_minus1_mer.ShiftAppend(c, globals.kmer_k - 1);
            rev_k_minus1_mer.ShiftPreappend(3 - c, globals.kmer_k - 1);
        }

        // last one special handling
        bucket_sizes[k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
        bucket_sizes[rev_k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
    }
    return NULL;
}

void s1_init_global_and_set_cx1(read2sdbg_global_t &globals) {
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
    // to count (k+1)-mers, sort by the internal (k-1)-mer
    // (k+1)-mer = abS[0..k-2]cd
    // is solid: number of bSc >= threshold
    // bS has in coming: for some a, num of abS >= threshold
    // Sc has outgoing: for some a, num of Scd >= threshold
    globals.words_per_substring = DivCeiling((globals.kmer_k - 1) * kBitsPerEdgeChar + 2 * kBWTCharNumBits, kBitsPerEdgeWord);
    // lv2 bytes: substring, permutation, readinfo
    int64_t lv2_bytes_per_item = (globals.words_per_substring) * sizeof(uint32_t) + sizeof(uint32_t) + sizeof(int64_t);
    lv2_bytes_per_item = lv2_bytes_per_item * 2; // double buffering
#ifndef USE_GPU
    lv2_bytes_per_item += sizeof(uint64_t) * 2; // CPU memory is used to simulate GPU
#endif

    if (cx1_t::kCX1Verbose >= 2) {
        xlog("%d words per substring\n", globals.words_per_substring);
    }
    // --- memory stuff ---
    int64_t mem_remained = globals.host_mem
                           - globals.mem_packed_reads
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
    globals.permutation = (uint32_t *) MallocAndCheck(globals.cx1.max_lv2_items_ * sizeof(uint32_t), __FILE__, __LINE__);
    globals.lv2_substrings_db = (uint32_t*) MallocAndCheck(globals.cx1.max_lv2_items_ * globals.words_per_substring * sizeof(uint32_t), __FILE__, __LINE__);
    globals.permutation_db = (uint32_t *) MallocAndCheck(globals.cx1.max_lv2_items_ * sizeof(uint32_t), __FILE__, __LINE__);
    globals.lv2_read_info = (int64_t *) MallocAndCheck(globals.cx1.max_lv2_items_ * sizeof(int64_t), __FILE__, __LINE__);
    globals.lv2_read_info_db = (int64_t *) MallocAndCheck(globals.cx1.max_lv2_items_ * sizeof(int64_t), __FILE__, __LINE__);
#ifndef USE_GPU
    globals.cpu_sort_space = (uint64_t*) MallocAndCheck(sizeof(uint64_t) * globals.cx1.max_lv2_items_, __FILE__, __LINE__);
#else
    alloc_gpu_buffers(globals.gpu_key_buffer1, globals.gpu_key_buffer2, globals.gpu_value_buffer1, globals.gpu_value_buffer2, (size_t)globals.cx1.max_lv2_items_);
#endif

    // --- initialize output mercy files ---
    globals.num_mercy_files = 1;
    while (globals.num_mercy_files * 10485760LL < globals.num_reads && globals.num_mercy_files < 64) {
        globals.num_mercy_files <<= 1;
    }
    if (cx1_t::kCX1Verbose >= 3) {
        xlog("Number of files for mercy candidate reads: %d\n", globals.num_mercy_files);
    }
    for (int i = 0; i < globals.num_mercy_files; ++i) {
        globals.mercy_files.push_back(OpenFileAndCheck(FormatString("%s.mercy_cand.%d", globals.output_prefix.c_str(), i), "wb"));
    }
    pthread_mutex_init(&globals.lv1_items_scanning_lock, NULL); // init lock

    // --- initialize stat ---
    globals.edge_counting = (int64_t *) MallocAndCheck((kMaxMulti_t + 1) * sizeof(int64_t), __FILE__, __LINE__);
    globals.thread_edge_counting = (int64_t *) MallocAndCheck((kMaxMulti_t + 1) * globals.num_output_threads * sizeof(int64_t), __FILE__, __LINE__);
    memset(globals.edge_counting, 0, (kMaxMulti_t + 1) * sizeof(int64_t));
}

void* s1_lv1_fill_offset(void* _data) {
    readpartition_data_t &rp = *((readpartition_data_t*) _data);
    read2sdbg_global_t &globals = *(rp.globals);
    int64_t *prev_full_offsets = (int64_t *)MallocAndCheck(kNumBuckets * sizeof(int64_t), __FILE__, __LINE__); // temporary array for computing differentials
    assert(prev_full_offsets != NULL);
    for (int b = globals.cx1.lv1_start_bucket_; b < globals.cx1.lv1_end_bucket_; ++b)
        prev_full_offsets[b] = rp.rp_lv1_differential_base;

    // this loop is VERY similar to that in PreprocessScanToFillBucketSizesThread
    GenericKmer k_minus1_mer, rev_k_minus1_mer; // (k+1)-mer and its rc

    int key;
    for (int64_t read_id = rp.rp_start_id; read_id < rp.rp_end_id; ++read_id) {
        int read_length = globals.package.length(read_id);
        if (read_length < globals.kmer_k + 1) {
            continue;
        }

        int64_t which_word = globals.package.get_start_index(read_id) / 16;
        int64_t offset = globals.package.get_start_index(read_id) % 16;
        uint32_t *read_p = &globals.package.packed_seq[which_word];

        k_minus1_mer.init(read_p, offset, globals.kmer_k - 1);
        rev_k_minus1_mer = k_minus1_mer;
        rev_k_minus1_mer.ReverseComplement(globals.kmer_k - 1);

        // ===== this is a macro to save some copy&paste ================
#define CHECK_AND_SAVE_OFFSET(offset, strand)                                                                      \
        do {                                                                                                       \
            assert(offset + globals.kmer_k - 1 <= read_length);                                                    \
            if (((key - globals.cx1.lv1_start_bucket_) ^ (key - globals.cx1.lv1_end_bucket_)) & kSignBitMask) {    \
                int64_t full_offset = EncodeOffset(read_id, offset, strand, globals.offset_num_bits);              \
                int64_t differential = full_offset - prev_full_offsets[key];                                       \
                if (differential > cx1_t::kDifferentialLimit) {                                                    \
                    pthread_mutex_lock(&globals.lv1_items_scanning_lock);                                          \
                    globals.lv1_items[rp.rp_bucket_offsets[key]++] = -globals.cx1.lv1_items_special_.size() - 1;   \
                    globals.cx1.lv1_items_special_.push_back(full_offset);                                         \
                    pthread_mutex_unlock(&globals.lv1_items_scanning_lock);                                        \
                } else {                                                                                           \
                    assert((int) differential >= 0);                                                               \
                    globals.lv1_items[rp.rp_bucket_offsets[key]++] = (int) differential;                           \
                }                                                                                                  \
                prev_full_offsets[key] = full_offset;                                                              \
            }                                                                                                      \
        } while (0)
        // ^^^^^ why is the macro surrounded by a do-while? please ask Google
        // =========== end macro ==========================

        // the first one special handling
        key = k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
        CHECK_AND_SAVE_OFFSET(0, 0);
        key = rev_k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
        CHECK_AND_SAVE_OFFSET(0, 1);

        int last_char_offset = globals.kmer_k - 1;
        int c = globals.package.get_base(read_id, last_char_offset);
        k_minus1_mer.ShiftAppend(c, globals.kmer_k - 1);
        rev_k_minus1_mer.ShiftPreappend(3 - c, globals.kmer_k - 1);

        // shift the key char by char
        while (last_char_offset < read_length - 1) {
            int cmp = k_minus1_mer.cmp(rev_k_minus1_mer, globals.kmer_k - 1);
            if (cmp > 0) {
                key = rev_k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
                CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k + 2, 1);
            } else if (cmp < 0) {
                key = k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
                CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k + 2, 0);
            } else {
                // a not-that-math-correct solution if the edge is palindrome, but works well enough
                int prev = globals.package.get_base(read_id, last_char_offset - (globals.kmer_k - 1));
                int next = globals.package.get_base(read_id, last_char_offset + 1);
                if (prev <= 3 - next) {
                    key = rev_k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
                    CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k + 2, 0);
                } else {
                    key = rev_k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
                    CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k + 2, 1);
                }
            }

            int c = globals.package.get_base(read_id, ++last_char_offset);
            k_minus1_mer.ShiftAppend(c, globals.kmer_k - 1);
            rev_k_minus1_mer.ShiftPreappend(3 - c, globals.kmer_k - 1);
        }

        // the last one special handling
        key = k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
        CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k + 2, 0);
        key = rev_k_minus1_mer.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
        CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k + 2, 1);
    }

#undef CHECK_AND_SAVE_OFFSET

    free(prev_full_offsets);
    return NULL;
}

void* s1_lv2_extract_substr(void* _data) {
    bucketpartition_data_t &bp = *((bucketpartition_data_t*) _data);
    read2sdbg_global_t &globals = *(bp.globals);
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
                int read_length = globals.package.length(read_id);
                int num_chars_to_copy = globals.kmer_k - 1;
                unsigned char prev, next, head, tail; // (k+1)=abScd, prev=a, head=b, tail=c, next=d

                assert(offset < read_length);

                if (offset > 1) {
                    head = globals.package.get_base(read_id, offset - 1);
                    prev = globals.package.get_base(read_id, offset - 2);
                } else {
                    prev = kSentinelValue;
                    if (offset > 0) {
                        head = globals.package.get_base(read_id, offset - 1);
                    } else {
                        head = kSentinelValue;
                    }
                }

                if (offset + globals.kmer_k < read_length) {
                    tail = globals.package.get_base(read_id, offset + globals.kmer_k - 1);
                    next = globals.package.get_base(read_id, offset + globals.kmer_k);
                } else {
                    next = kSentinelValue;
                    if (offset + globals.kmer_k - 1 < read_length) {
                        tail = globals.package.get_base(read_id, offset + globals.kmer_k - 1);
                    } else {
                        tail = kSentinelValue;
                    }
                }

                int64_t which_word = globals.package.get_start_index(read_id) / 16;
                int start_offset = globals.package.get_start_index(read_id) % 16;
                uint32_t *read_p = &globals.package.packed_seq[which_word];
                int words_this_read = DivCeiling(start_offset + read_length, 16);

                if (strand == 0) {
                    CopySubstring(substrings_p, read_p, offset + start_offset, num_chars_to_copy,
                                  globals.cx1.lv2_num_items_, words_this_read, globals.words_per_substring);
                    uint32_t *last_word = substrings_p + int64_t(globals.words_per_substring - 1) * globals.cx1.lv2_num_items_;
                    *last_word |= (head << kBWTCharNumBits) | tail;
                    *read_info_p = (full_offset << 6) | (prev << 3) | next;
                } else {
                    CopySubstringRC(substrings_p, read_p, offset + start_offset, num_chars_to_copy,
                                    globals.cx1.lv2_num_items_, words_this_read, globals.words_per_substring);
                    uint32_t *last_word = substrings_p + int64_t(globals.words_per_substring - 1) * globals.cx1.lv2_num_items_;
                    *last_word |= ((tail == kSentinelValue ? kSentinelValue : 3 - tail) << kBWTCharNumBits) | (head == kSentinelValue ? kSentinelValue : 3 - head);
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

void s1_lv2_sort(read2sdbg_global_t &globals) {
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

void s1_lv2_pre_output_partition(read2sdbg_global_t &globals) {
    // swap double buffers
    globals.lv2_num_items_db = globals.cx1.lv2_num_items_;
    std::swap(globals.lv2_substrings_db, globals.lv2_substrings);
    std::swap(globals.permutation_db, globals.permutation);
    std::swap(globals.lv2_read_info_db, globals.lv2_read_info);

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
                uint32_t *prev_item = globals.lv2_substrings_db + globals.permutation_db[this_end_index - 1];
                uint32_t *item = globals.lv2_substrings_db + globals.permutation_db[this_end_index];
                if (IsDiffKMinusOneMer(prev_item, item, globals.lv2_num_items_db, globals.kmer_k)) {
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

void* s1_lv2_output(void* _op) {
    xtimer_t local_timer;
    if (cx1_t::kCX1Verbose >= 4) {
        local_timer.start();
        local_timer.reset();
    }
    outputpartition_data_t *op = (outputpartition_data_t*) _op;
    read2sdbg_global_t &globals = *(op->globals);
    int64_t op_start_index = op->op_start_index;
    int64_t op_end_index = op->op_end_index;
    int thread_id = op->op_id;
    int end_idx;
    int count_prev_head[5][5];
    int count_tail_next[5][5];
    int count_head_tail[(1 << 2 * kBWTCharNumBits) - 1];
    int64_t offset_mask = (1 << globals.offset_num_bits) - 1;
    int64_t *thread_edge_counting = globals.thread_edge_counting + thread_id * (kMaxMulti_t + 1);

    for (int i = op_start_index; i < op_end_index; i = end_idx) {
        end_idx = i + 1;
        uint32_t *first_item = globals.lv2_substrings_db + (globals.permutation_db[i]);
        memset(count_prev_head, 0, sizeof(count_prev_head));
        memset(count_tail_next, 0, sizeof(count_tail_next));
        memset(count_head_tail, 0, sizeof(count_head_tail));

        {
            uint8_t prev_and_next = ExtractPrevNext(globals.permutation_db[i], globals);
            uint8_t head_and_tail = ExtractHeadTail(globals.lv2_substrings_db + globals.permutation_db[i], globals.lv2_num_items_db, globals.words_per_substring);
            count_prev_head[prev_and_next >> 3][head_and_tail >> 3]++;
            count_tail_next[head_and_tail & 7][prev_and_next & 7]++;
            count_head_tail[head_and_tail]++;
        }

        while (end_idx < op_end_index) {
            if (IsDiffKMinusOneMer(first_item,
                                   globals.lv2_substrings_db + globals.permutation_db[end_idx],
                                   globals.lv2_num_items_db,
                                   globals.kmer_k)) {
                break;
            }

            uint8_t prev_and_next = ExtractPrevNext(globals.permutation_db[end_idx], globals);
            uint8_t head_and_tail = ExtractHeadTail(globals.lv2_substrings_db + globals.permutation_db[end_idx], globals.lv2_num_items_db, globals.words_per_substring);
            count_prev_head[prev_and_next >> 3][head_and_tail >> 3]++;
            count_tail_next[head_and_tail & 7][prev_and_next & 7]++;
            count_head_tail[head_and_tail]++;

            ++end_idx;
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

        int l_has_out = 0, r_has_in = 0;
        for (int j = 0; j < 4; ++j) {
            for (int x = 0; x < 4; ++x) {
                if (count_head_tail[(j << kBWTCharNumBits) | x] >= globals.kmer_freq_threshold) {
                    l_has_out |= 1 << j;
                    r_has_in |= 1 << x;
                }
            }
        }

        while (i < end_idx) {
            uint8_t head_and_tail = ExtractHeadTail(globals.lv2_substrings_db + globals.permutation_db[i], globals.lv2_num_items_db, globals.words_per_substring);
            uint8_t head = head_and_tail >> 3;
            uint8_t tail = head_and_tail & 7;
            if (head != kSentinelValue && tail != kSentinelValue) {
                ++thread_edge_counting[std::min(kMaxMulti_t, count_head_tail[head_and_tail])];
            }

            if (head != kSentinelValue && tail != kSentinelValue && count_head_tail[head_and_tail] >= globals.kmer_freq_threshold) {
                for (int j = 0; j < count_head_tail[head_and_tail]; ++j, ++i) {
                    int64_t read_info = globals.lv2_read_info_db[globals.permutation_db[i]] >> 6;
                    int strand = read_info & 1;
                    int offset = ((read_info >> 1) & offset_mask) - 1;
                    int64_t read_id = read_info >> (1 + globals.offset_num_bits);
                    int l_offset = strand == 0 ? offset : offset + 1;
                    int r_offset = strand == 0 ? offset + 1 : offset;

                    // mark this is a solid edge
                    globals.is_solid.set(globals.num_k1_per_read * read_id + offset);

                    if (!(has_in & (1 << head))) {
                        // no in
                        int64_t packed_mercy_cand = (read_id << (globals.offset_num_bits + 2)) | (l_offset << 2) | (1 + strand);
                        fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1, globals.mercy_files[read_id & (globals.num_mercy_files - 1)]);
                    }

                    if (!(has_out & (1 << tail))) {
                        // no out
                        int64_t packed_mercy_cand = (read_id << (globals.offset_num_bits + 2)) | (r_offset << 2) | (2 - strand);
                        fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1, globals.mercy_files[read_id & (globals.num_mercy_files - 1)]);
                    }
                }
            } else {
                // not solid, but we still need to tell whether its left/right kmer is solid
                for (int j = 0; j < count_head_tail[head_and_tail]; ++j, ++i) {
                    int64_t read_info = globals.lv2_read_info_db[globals.permutation_db[i]] >> 6;
                    int strand = read_info & 1;
                    int offset = ((read_info >> 1) & offset_mask) - 1;
                    int64_t read_id = read_info >> (1 + globals.offset_num_bits);
                    int l_offset = strand == 0 ? offset : offset + 1;
                    int r_offset = strand == 0 ? offset + 1 : offset;

                    if (l_has_out & (1 << head)) {
                        if (has_in & (1 << head)) {
                            // has both in & out
                            int64_t packed_mercy_cand = (read_id << (globals.offset_num_bits + 2)) | (l_offset << 2) | 0;
                            fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1, globals.mercy_files[read_id & (globals.num_mercy_files - 1)]);
                        } else {
                            // has out but no in
                            int64_t packed_mercy_cand = (read_id << (globals.offset_num_bits + 2)) | (l_offset << 2) | (1 + strand);
                            fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1, globals.mercy_files[read_id & (globals.num_mercy_files - 1)]);
                        }
                    } else {
                        if (has_in & (1 << head)) {
                            // has in but no out
                            int64_t packed_mercy_cand = (read_id << (globals.offset_num_bits + 2)) | (l_offset << 2) | (2 - strand);
                            fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1, globals.mercy_files[read_id & (globals.num_mercy_files - 1)]);
                        }
                    }

                    if (r_has_in & (1 << tail)) {
                        if (has_out & (1 << tail)) {
                            // has both in & out
                            int64_t packed_mercy_cand = (read_id << (globals.offset_num_bits + 2)) | (r_offset << 2) | 0;
                            fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1, globals.mercy_files[read_id & (globals.num_mercy_files - 1)]);
                        } else {
                            // has in but no out
                            int64_t packed_mercy_cand = (read_id << (globals.offset_num_bits + 2)) | (r_offset << 2) | (2 - strand);
                            fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1, globals.mercy_files[read_id & (globals.num_mercy_files - 1)]);
                        }
                    } else {
                        if (has_out & (1 << tail)) {
                            // has out but no in
                            int64_t packed_mercy_cand = (read_id << (globals.offset_num_bits + 2)) | (r_offset << 2) | (1 + strand);
                            fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1, globals.mercy_files[read_id & (globals.num_mercy_files - 1)]);
                        }
                    }
                }
            }
        }

        if (i != end_idx) {
            xerr_and_exit("%d %d %d %d\n", i, end_idx, op_start_index, op_end_index);
        }
    }

    if (cx1_t::kCX1Verbose >= 4) {
        local_timer.stop();
        xlog("Counting time elapsed: %.4lfs\n", local_timer.elapsed());
    }
    return NULL;
}

void s1_lv2_post_output(read2sdbg_global_t &globals) {
    for (int t = 0; t < globals.num_output_threads; ++t) {
        for (int i = 1; i <= kMaxMulti_t; ++i) {
            globals.edge_counting[i] += globals.thread_edge_counting[t * (kMaxMulti_t + 1) + i];
        }
    }
}

void s1_post_proc(read2sdbg_global_t &globals) {
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

    // --- cleaning ---
    pthread_mutex_destroy(&globals.lv1_items_scanning_lock);
    free(globals.lv1_items);
    free(globals.lv2_substrings);
    free(globals.permutation);
    free(globals.permutation_db);
    free(globals.lv2_substrings_db);
    free(globals.lv2_read_info);
    free(globals.lv2_read_info_db);
    free(globals.edge_counting);
    free(globals.thread_edge_counting);
    for (int i = 0; i < globals.num_mercy_files; ++i) {
        fclose(globals.mercy_files[i]);
    }

#ifndef USE_GPU
    free(globals.cpu_sort_space);
#else
    free_gpu_buffers(globals.gpu_key_buffer1, globals.gpu_key_buffer2, globals.gpu_value_buffer1, globals.gpu_value_buffer2);
#endif
}

} // s1

} // cx1_read2sdbg