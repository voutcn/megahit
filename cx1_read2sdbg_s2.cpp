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

#include "cx1_read2sdbg.h"

#include <omp.h>
#include <string>
#include <vector>
#include <parallel/algorithm>

#include "utils.h"
#include "kmer.h"
#include "mem_file_checker-inl.h"
#include "packed_reads.h"
#include "read_lib_functions-inl.h"

#include "lv2_cpu_sort.h"
#include "lv2_gpu_functions.h"

extern void kt_dfor(int n_threads, void (*func)(void *, long, int), void *data, long n);

namespace cx1_read2sdbg {

namespace s2 {

typedef CX1<read2sdbg_global_t, kNumBuckets> cx1_t;
typedef CX1<read2sdbg_global_t, kNumBuckets>::readpartition_data_t readpartition_data_t;
typedef CX1<read2sdbg_global_t, kNumBuckets>::bucketpartition_data_t bucketpartition_data_t;
typedef CX1<read2sdbg_global_t, kNumBuckets>::outputpartition_data_t outputpartition_data_t;

// helper functions
inline int64_t EncodeOffset(int64_t read_id, int offset, int strand, SequencePackage &p, int edge_type) {
    // edge_type: 0 left $; 1 solid; 2 right $
    return ((p.get_start_index(read_id) + offset) << 3) | (edge_type << 1) | strand;
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

// helper
inline int ExtractFirstChar(uint32_t *item) {
    return *item >> kTopCharShift;
}

// bS'a
inline int Extract_a(uint32_t *item, int num_words, int64_t spacing, int kmer_k) {
    int non_dollar = (item[(num_words - 1) * spacing] >> kBWTCharNumBits) & 1;

    if (non_dollar) {
        int which_word = (kmer_k - 1) / kCharsPerEdgeWord;
        int word_index = (kmer_k - 1) % kCharsPerEdgeWord;
        return (item[which_word * spacing] >> (kCharsPerEdgeWord - 1 - word_index) * kBitsPerEdgeChar) & kEdgeCharMask;
    }
    else {
        return kSentinelValue;
    }
}

inline int Extract_b(uint32_t *item, int num_words, int64_t spacing) {
    return item[(num_words - 1) * spacing] & ((1 << kBWTCharNumBits) - 1);
}


// cx1 core functions
int64_t s2_encode_lv1_diff_base(int64_t read_id, read2sdbg_global_t &globals) {
    return EncodeOffset(read_id, 0, 0, globals.package, 0);
}

void s2_read_mercy_prepare(read2sdbg_global_t &globals) {
    if (!globals.need_mercy || globals.kmer_freq_threshold == 1) return;

    xtimer_t timer;

    if (cx1_t::kCX1Verbose >= 3) {
        timer.reset();
        timer.start();
        xlog("Adding mercy edges...\n");
    }

    std::vector<uint64_t> mercy_cand;
    uint64_t num_mercy = 0;
    AtomicBitVector read_marker;
    read_marker.reset(globals.num_short_reads);

    for (int fid = 0; fid < globals.num_mercy_files; ++fid) {
        FILE *fp = OpenFileAndCheck(FormatString("%s.mercy_cand.%d", globals.output_prefix.c_str(), fid), "rb");
        mercy_cand.clear();

        int num_read = 0;
        uint64_t buf[4096];

        while ((num_read = fread(buf, sizeof(uint64_t), 4096, fp)) > 0) {
            mercy_cand.insert(mercy_cand.end(), buf, buf + num_read);
        }

        if (cx1_t::kCX1Verbose >= 4) {
            xlog("Mercy file: %s, %lu\n", FormatString("%s.mercy_cand.%d", globals.output_prefix.c_str(), fid), mercy_cand.size());
        }

        omp_set_num_threads(globals.num_cpu_threads);
        __gnu_parallel::sort(mercy_cand.begin(), mercy_cand.end());

        // multi threading
        uint64_t avg = DivCeiling(mercy_cand.size(), globals.num_cpu_threads);
        std::vector<uint64_t> start_idx(globals.num_cpu_threads);
        std::vector<uint64_t> end_idx(globals.num_cpu_threads);

        // manually distribute threads
        for (int tid = 0; tid < globals.num_cpu_threads; ++tid) {
            if (tid == 0) {
                start_idx[tid] = 0;
            }
            else {
                start_idx[tid] = end_idx[tid - 1];
            }

            uint64_t this_end = std::min(start_idx[tid] + avg, (uint64_t)mercy_cand.size());
            uint64_t read_id = globals.package.get_id(mercy_cand[this_end] >> 2);

            while (this_end < mercy_cand.size() && globals.package.get_id(mercy_cand[this_end] >> 2) == read_id) {
                ++this_end;
            }

            end_idx[tid] = this_end;
        }

        #pragma omp parallel for reduction(+:num_mercy)

        for (int tid = 0; tid < globals.num_cpu_threads; ++tid) {
            std::vector<bool> no_in(globals.max_read_length);
            std::vector<bool> no_out(globals.max_read_length);
            std::vector<bool> has_solid_kmer(globals.max_read_length);

            uint64_t i = start_idx[tid];

            // go read by read
            while (i != end_idx[tid]) {
                uint64_t read_id = globals.package.get_id(mercy_cand[i] >> 2);
                assert(!read_marker.get(read_id));
                read_marker.set(read_id);
                int first_0_out = globals.max_read_length + 1;
                int last_0_in = -1;

                std::fill(no_in.begin(), no_in.end(), false);
                std::fill(no_out.begin(), no_out.end(), false);
                std::fill(has_solid_kmer.begin(), has_solid_kmer.end(), false);

                while (i != end_idx[tid] && globals.package.get_id(mercy_cand[i] >> 2) == read_id) {
                    int offset = (mercy_cand[i] >> 2) - globals.package.get_start_index(read_id);
                    if ((mercy_cand[i] & 3) == 2) {
                        no_out[offset] = true;
                        first_0_out = std::min(first_0_out, offset);
                    }
                    else if ((mercy_cand[i] & 3) == 1) {
                        no_in[offset] = true;
                        last_0_in = std::max(last_0_in, offset);
                    }

                    has_solid_kmer[offset] = true;
                    ++i;
                }

                if (last_0_in < first_0_out) {
                    continue;
                }

                int read_length = globals.package.length(read_id);
                int last_no_out = -1;

                for (int i = 0; i + globals.kmer_k < read_length; ++i) {
                    if (globals.is_solid.get(globals.package.get_start_index(read_id) + i)) {
                        has_solid_kmer[i] = has_solid_kmer[i + 1] = true;
                    }
                }

                for (int i = 0; i + globals.kmer_k <= read_length; ++i) {
                    if (no_in[i] && last_no_out != -1) {
                        for (int j = last_no_out; j < i; ++j) {
                            globals.is_solid.set(globals.package.get_start_index(read_id) + j);
                        }

                        num_mercy += i - last_no_out;
                    }

                    if (has_solid_kmer[i]) {
                        last_no_out = -1;
                    }

                    if (no_out[i]) {
                        last_no_out = i;
                    }
                }

            }
        }

        fclose(fp);
    }

    if (cx1_t::kCX1Verbose >= 3) {
        timer.stop();
        xlog("Adding mercy Done. Time elapsed: %.4lf\n", timer.elapsed());
        xlog("Number mercy: %llu\n", (unsigned long long)num_mercy);
    }

    // set cx1 param
    globals.cx1.num_cpu_threads_ = globals.num_cpu_threads;
    globals.cx1.num_output_threads_ = globals.num_output_threads;
    globals.cx1.num_items_ = globals.num_reads;
}

void *s2_lv0_calc_bucket_size(void *_data) {
    readpartition_data_t &rp = *((readpartition_data_t *) _data);
    read2sdbg_global_t &globals = *(rp.globals);
    int64_t *bucket_sizes = rp.rp_bucket_sizes;
    memset(bucket_sizes, 0, kNumBuckets * sizeof(int64_t));
    GenericKmer edge, rev_edge; // (k+1)-mer and its rc

    for (int64_t read_id = rp.rp_start_id; read_id < rp.rp_end_id; ++read_id) {
        int read_length = globals.package.length(read_id);

        if (read_length < globals.kmer_k + 1) {
            continue;
        }

        int64_t which_word = globals.package.get_start_index(read_id) / 16;
        int64_t offset = globals.package.get_start_index(read_id) % 16;
        uint32_t *read_p = &globals.package.packed_seq[which_word];

        edge.init(read_p, offset, globals.kmer_k + 1);
        rev_edge = edge;
        rev_edge.ReverseComplement(globals.kmer_k + 1);

        int last_char_offset = globals.kmer_k;
        int64_t full_offset = globals.package.get_start_index(read_id);
        bool is_solid = globals.kmer_freq_threshold == 1 || read_id >= globals.num_short_reads;

        while (true) {
            if (is_solid || globals.is_solid.get(full_offset)) {
                bool is_palindrome = rev_edge.cmp(edge, globals.kmer_k + 1) == 0;
                bucket_sizes[(edge.data_[0] << 2) >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;

                if (!is_palindrome)
                    bucket_sizes[(rev_edge.data_[0] << 2) >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;

                if (last_char_offset == globals.kmer_k || !(is_solid || globals.is_solid.get(full_offset - 1))) {
                    bucket_sizes[edge.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;

                    if (!is_palindrome)
                        bucket_sizes[(rev_edge.data_[0] << 4) >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
                }

                if (last_char_offset == read_length - 1 || !(is_solid || globals.is_solid.get(full_offset + 1))) {
                    bucket_sizes[(edge.data_[0] << 4) >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;

                    if (!is_palindrome)
                        bucket_sizes[rev_edge.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
                }
            }

            ++full_offset;

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

void s2_init_global_and_set_cx1(read2sdbg_global_t &globals) {

    globals.max_bucket_size = *std::max_element(globals.cx1.bucket_sizes_, globals.cx1.bucket_sizes_ + kNumBuckets);
    globals.tot_bucket_size = 0;
    int num_non_empty = 0;

    for (int i = 0; i < kNumBuckets; ++i) {
        globals.tot_bucket_size += globals.cx1.bucket_sizes_[i];

        if (globals.cx1.bucket_sizes_[i] > 0) {
            num_non_empty++;
        }
    }

    globals.words_per_substring = DivCeiling(globals.kmer_k * kBitsPerEdgeChar + kBWTCharNumBits + 1, kBitsPerEdgeWord);
    globals.words_per_dummy_node = DivCeiling(globals.kmer_k * kBitsPerEdgeChar, kBitsPerEdgeWord);

    if (cx1_t::kCX1Verbose >= 2) {
        xlog("%d words per substring, words per dummy node ($v): %d\n", globals.words_per_substring, globals.words_per_dummy_node);
    }

    // --- calculate lv2 memory ---
#ifdef USE_GPU
    int64_t lv2_mem = globals.gpu_mem - 1073741824; // should reserver ~1G for GPU sorting
    globals.cx1.max_lv2_items_ = std::min(lv2_mem / cx1_t::kGPUBytePerItem, std::max(globals.max_bucket_size, kMinLv2BatchSizeGPU));

    if (globals.max_bucket_size > globals.cx1.max_lv2_items_) {
        xerr_and_exit("Bucket too large for GPU: contains %lld items. Please try CPU version.\n", globals.max_bucket_size);
        // TODO: auto switch to CPU version
    }

    // lv2 bytes: substring (double buffer), permutation, aux
    int64_t lv2_bytes_per_item = (globals.words_per_substring * sizeof(uint32_t) + sizeof(uint32_t)) * 2 + sizeof(int64_t);

    // --- memory stuff ---
    int64_t mem_remained = globals.host_mem
                           - globals.mem_packed_reads
                           - kNumBuckets * sizeof(int64_t) * (globals.num_cpu_threads * 3 + 1);

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
    globals.permutation = (uint32_t *) MallocAndCheck(globals.cx1.max_lv2_items_ * sizeof(uint32_t), __FILE__, __LINE__);
    globals.lv2_substrings_db = (uint32_t *) MallocAndCheck(globals.cx1.max_lv2_items_ * globals.words_per_substring * sizeof(uint32_t), __FILE__, __LINE__);
    globals.permutation_db = (uint32_t *) MallocAndCheck(globals.cx1.max_lv2_items_ * sizeof(uint32_t), __FILE__, __LINE__);
    alloc_gpu_buffers(globals.gpu_key_buffer1, globals.gpu_key_buffer2, globals.gpu_value_buffer1, globals.gpu_value_buffer2, (size_t)globals.cx1.max_lv2_items_);

#else

    num_non_empty = std::max(1, num_non_empty);

    for (int i = 0; i < kNumBuckets; ++i) {
        if (globals.cx1.bucket_sizes_[i] > 2 * globals.tot_bucket_size / num_non_empty) {
            // xlog("Bucket %d size = %lld > %lld = 2 * avg\n", i, (long long)globals.cx1.bucket_sizes_[i], (long long)2 * globals.tot_bucket_size / num_non_empty);
        }
    }

    int64_t lv2_bytes_per_item = globals.words_per_substring * sizeof(uint32_t) + sizeof(uint32_t) + sizeof(uint32_t);

    globals.max_sorting_items = std::max(globals.tot_bucket_size * globals.num_cpu_threads / num_non_empty, globals.max_bucket_size * 2);
    globals.cx1.lv1_just_go_ = true;
    globals.num_output_threads = globals.num_cpu_threads;

    int64_t mem_remained = globals.host_mem
                           - globals.mem_packed_reads
                           - globals.num_cpu_threads * 65536 * sizeof(uint64_t) // radix sort buckets
                           - kNumBuckets * sizeof(int64_t) * (globals.num_cpu_threads * 3 + 1);
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

    globals.cx1.max_mem_remain_ = globals.cx1.max_lv1_items_ * sizeof(int) + globals.max_sorting_items * lv2_bytes_per_item;
    globals.cx1.bytes_per_sorting_item_ = lv2_bytes_per_item;

    globals.lv1_items = (int *) MallocAndCheck(globals.cx1.max_mem_remain_ + globals.num_cpu_threads * sizeof(uint64_t) * 65536, __FILE__, __LINE__);

#endif

    if (cx1_t::kCX1Verbose >= 2) {
        xlog("Memory for sequence: %lld\n", globals.mem_packed_reads);
        xlog("max # lv.1 items = %lld\n", globals.cx1.max_lv1_items_);
#ifdef USE_GPU
        xlog("max # lv.2 items = %lld\n", globals.cx1.max_lv2_items_);
#endif
    }

    pthread_mutex_init(&globals.lv1_items_scanning_lock, NULL); // init lock
    // --- init output ---
    globals.sdbg_writer.set_num_threads(globals.num_output_threads);
    globals.sdbg_writer.set_kmer_size(globals.kmer_k);
    globals.sdbg_writer.set_num_buckets(kNumBuckets);
    globals.sdbg_writer.set_file_prefix(globals.output_prefix);
    globals.sdbg_writer.init_files();
}

void *s2_lv1_fill_offset(void *_data) {
    readpartition_data_t &rp = *((readpartition_data_t *) _data);
    read2sdbg_global_t &globals = *(rp.globals);
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

        int64_t which_word = globals.package.get_start_index(read_id) / 16;
        int64_t offset = globals.package.get_start_index(read_id) % 16;
        uint32_t *read_p = &globals.package.packed_seq[which_word];

        edge.init(read_p, offset, globals.kmer_k + 1);
        rev_edge = edge;
        rev_edge.ReverseComplement(globals.kmer_k + 1);

        // ===== this is a macro to save some copy&paste ================
#define CHECK_AND_SAVE_OFFSET(offset, strand, edge_type)                                                                \
    do {                                                                                                                \
        if (globals.cx1.cur_lv1_buckets_[key]) {                                                                        \
            int key_ = globals.cx1.bucket_rank_[key];                                                                   \
            int64_t full_offset = EncodeOffset(read_id, offset, strand, globals.package, edge_type);                    \
            int64_t differential = full_offset - prev_full_offsets[key_];                                               \
            if (differential > cx1_t::kDifferentialLimit) {                                                             \
                pthread_mutex_lock(&globals.lv1_items_scanning_lock);                                                   \
                globals.lv1_items[rp.rp_bucket_offsets[key_]++] = -globals.cx1.lv1_items_special_.size() - 1;           \
                globals.cx1.lv1_items_special_.push_back(full_offset);                                                  \
                pthread_mutex_unlock(&globals.lv1_items_scanning_lock);                                                 \
            } else {                                                                                                    \
                assert ((int) differential >= 0);                                                                       \
                globals.lv1_items[rp.rp_bucket_offsets[key_]++] = (int) differential;                                   \
            }                                                                                                           \
            prev_full_offsets[key_] = full_offset;                                                                      \
        }                                                                                                               \
    } while (0)
        // ^^^^^ why is the macro surrounded by a do-while? please ask Google
        // =========== end macro ==========================

        // shift the key char by char
        int last_char_offset = globals.kmer_k;
        int64_t full_offset = globals.package.get_start_index(read_id);
        bool is_solid = globals.kmer_freq_threshold == 1 || read_id >= globals.num_short_reads;

        while (true) {
            if (is_solid || globals.is_solid.get(full_offset)) {
                bool is_palindrome = rev_edge.cmp(edge, globals.kmer_k + 1) == 0;

                // left $
                if (last_char_offset == globals.kmer_k || !(is_solid || globals.is_solid.get(full_offset - 1))) {
                    key = edge.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
                    CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 0, 0);

                    if (!is_palindrome) {
                        key = (rev_edge.data_[0] << 4) >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
                        CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 1, 0);
                    }
                }

                // solid
                key = (edge.data_[0] << 2) >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
                CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 0, 1);

                if (!is_palindrome) {
                    key = (rev_edge.data_[0] << 2) >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
                    CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 1, 1);
                }

                // right $
                if (last_char_offset == read_length - 1 || !(is_solid || globals.is_solid.get(full_offset + 1))) {
                    key = (edge.data_[0] << 4) >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
                    CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 0, 2);

                    if (!is_palindrome) {
                        key = rev_edge.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
                        CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 1, 2);
                    }
                }
            }

            ++full_offset;

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

void s2_lv2_extract_substr_(int bp_from, int bp_to, read2sdbg_global_t &globals, uint32_t *substr, int num_items) {
    int *lv1_p = globals.lv1_items + globals.cx1.rp_[0].rp_bucket_offsets[bp_from];

    for (int b = bp_from; b < bp_to; ++b) {
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

                int64_t read_id = globals.package.get_id(full_offset >> 3);
                int offset = (full_offset >> 3) - globals.package.get_start_index(read_id);
                int strand = full_offset & 1;
                int edge_type = (full_offset >> 1) & 3;
                int read_length = globals.package.length(read_id);

                int64_t which_word = globals.package.get_start_index(read_id) / 16;
                int start_offset = globals.package.get_start_index(read_id) % 16;
                uint32_t *read_p = &globals.package.packed_seq[which_word];
                int words_this_read = DivCeiling(start_offset + read_length, 16);

                if (strand == 0) {
                    int num_chars_to_copy = globals.kmer_k;
                    uint8_t prev = kSentinelValue;

                    switch (edge_type) {
                    case 0:
                        break;

                    case 1:
                        prev = globals.package.get_base(read_id, offset);
                        offset++;
                        break;

                    case 2:
                        prev = globals.package.get_base(read_id, offset + 1);
                        offset += 2;
                        num_chars_to_copy--;
                        break;

                    default:
                        assert(false);
                    }

                    CopySubstring(substr, read_p, offset + start_offset, num_chars_to_copy,
                                  num_items, words_this_read, globals.words_per_substring);

                    uint32_t *last_word = substr + int64_t(globals.words_per_substring - 1) * num_items;
                    *last_word |= int(num_chars_to_copy == globals.kmer_k) << kBWTCharNumBits;
                    *last_word |= prev;
                }
                else {
                    int num_chars_to_copy = globals.kmer_k;
                    uint8_t prev = kSentinelValue;

                    switch (edge_type) {
                    case 0:
                        num_chars_to_copy--;
                        prev = 3 - globals.package.get_base(read_id, offset + globals.kmer_k - 1);
                        break;

                    case 1:
                        prev = 3 - globals.package.get_base(read_id, offset + globals.kmer_k);
                        break;

                    case 2:
                        offset++;
                        break;

                    default:
                        assert(false);
                    }

                    CopySubstringRC(substr, read_p, offset + start_offset, num_chars_to_copy,
                                    num_items, words_this_read, globals.words_per_substring);

                    uint32_t *last_word = substr + int64_t(globals.words_per_substring - 1) * num_items;
                    *last_word |= int(num_chars_to_copy == globals.kmer_k) << kBWTCharNumBits;
                    *last_word |= prev;
                }

                substr++;
            }
        }
    }
}

void *s2_lv2_extract_substr(void *_data) {
    bucketpartition_data_t &bp = *((bucketpartition_data_t *) _data);
    read2sdbg_global_t &globals = *(bp.globals);
    uint32_t *substr = globals.lv2_substrings +
                       (globals.cx1.rp_[0].rp_bucket_offsets[bp.bp_start_bucket] - globals.cx1.rp_[0].rp_bucket_offsets[globals.cx1.lv2_start_bucket_]);
    s2_lv2_extract_substr_(bp.bp_start_bucket, bp.bp_end_bucket, globals, substr, globals.cx1.lv2_num_items_);
    return NULL;
}

void s2_lv2_sort(read2sdbg_global_t &globals) {
    xtimer_t local_timer;
#ifdef USE_GPU

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

void s2_lv2_pre_output_partition(read2sdbg_global_t &globals) {
    // swap double buffers
    globals.lv2_num_items_db = globals.cx1.lv2_num_items_;
    std::swap(globals.lv2_substrings_db, globals.lv2_substrings);
    std::swap(globals.permutation_db, globals.permutation);

    // distribute threads
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

void output_(int64_t from, int64_t to, read2sdbg_global_t &globals, uint32_t *substr, uint32_t *permutation, int tid, int num_items) {
    int start_idx, end_idx;
    int has_solid_a = 0; // has solid (k+1)-mer aSb
    int has_solid_b = 0; // has solid aSb
    int last_a[4], outputed_b;
    uint32_t tip_label[32];

    for (start_idx = from; start_idx < to; start_idx = end_idx) {
        end_idx = start_idx + 1;
        uint32_t *item = substr + permutation[start_idx];

        while (end_idx < to &&
                !IsDiffKMinusOneMer(
                    item,
                    substr + permutation[end_idx],
                    num_items,
                    globals.kmer_k)) {
            ++end_idx;
        }

        // clean marking
        has_solid_a = has_solid_b = 0;
        outputed_b = 0;

        for (int i = start_idx; i < end_idx; ++i) {
            uint32_t *cur_item = substr + permutation[i];
            int a = Extract_a(cur_item, globals.words_per_substring, num_items, globals.kmer_k);
            int b = Extract_b(cur_item, globals.words_per_substring, num_items);

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
            uint32_t *cur_item = substr + permutation[i];
            int a = Extract_a(cur_item, globals.words_per_substring, num_items, globals.kmer_k);
            int b = Extract_b(cur_item, globals.words_per_substring, num_items);

            j = i + 1;

            while (j < end_idx) {
                uint32_t *next_item = substr + permutation[j];

                if (Extract_a(next_item, globals.words_per_substring, num_items, globals.kmer_k) != a ||
                        Extract_b(next_item, globals.words_per_substring, num_items) != b) {
                    break;
                }
                else {
                    ++j;
                }
            }

            int w, last, is_dollar = 0;
            int count = std::min(j - i, kMaxMulti_t);;

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
            last = (a == kSentinelValue) ? 0 : ((last_a[a] == j - 1) ? 1 : 0);
            outputed_b |= 1 << b;

            if (is_dollar) {
                for (int64_t i = 0; i < globals.words_per_dummy_node; ++i) {
                    tip_label[i] = cur_item[i * num_items];
                }
            }

            globals.sdbg_writer.write(tid, cur_item[0] >> (32 - kBucketPrefixLength * 2), w, last, is_dollar, count, tip_label);
        }
    }
}

void *s2_lv2_output(void *_op) {
    outputpartition_data_t *op = (outputpartition_data_t *) _op;
    read2sdbg_global_t &globals = *(op->globals);
    int64_t op_start_index = op->op_start_index;
    int64_t op_end_index = op->op_end_index;

    output_(op_start_index, op_end_index, globals, globals.lv2_substrings_db, globals.permutation_db, op->op_id, globals.lv2_num_items_db);

    return NULL;
}

void s2_lv2_post_output(read2sdbg_global_t &globals) {
}

struct kt_sort_t {
    read2sdbg_global_t *globals;
    int64_t acc_size;
    int activated_threads;
    std::vector<int64_t> thread_offset;
    std::vector<int> tid_map;
    volatile int lock_;
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

    s2_lv2_extract_substr_(b, b + 1, *(kg->globals), substr_ptr, kg->globals->cx1.bucket_sizes_[b]);
    lv2_cpu_radix_sort_st(substr_ptr, permutation_ptr, cpu_sort_space_ptr, bucket, kg->globals->words_per_substring, kg->globals->cx1.bucket_sizes_[b]);
    output_(0, kg->globals->cx1.bucket_sizes_[b], *(kg->globals), substr_ptr, permutation_ptr, tid, kg->globals->cx1.bucket_sizes_[b]);
}

void s2_lv1_direct_sort_and_proc(read2sdbg_global_t &globals) {
    kt_sort_t kg;
    kg.globals = &globals;
    int64_t acc_size = 0;

    for (int i = 0, b = globals.cx1.lv1_start_bucket_; i < globals.num_cpu_threads && b < globals.cx1.lv1_end_bucket_; ++i, ++b) {
        kg.thread_offset.push_back(acc_size);
        acc_size += globals.cx1.bucket_sizes_[b];
    }

    kt_dfor(globals.num_cpu_threads, kt_sort, &kg, globals.cx1.lv1_end_bucket_ - globals.cx1.lv1_start_bucket_);
}

void s2_post_proc(read2sdbg_global_t &globals) {
    if (cx1_t::kCX1Verbose >= 2) {
        xlog("Number of $ A C G T A- C- G- T-:\n");
    }

    xlog("");

    for (int i = 0; i < 9; ++i) {
        xlog_ext("%lld ", (long long)globals.sdbg_writer.num_w(i));
    }

    xlog_ext("\n");

    if (cx1_t::kCX1Verbose >= 2) {
        xlog("Total number of edges: %lld\n", (long long)globals.sdbg_writer.num_edges());
        xlog("Total number of ONEs: %lld\n", (long long)globals.sdbg_writer.num_last1());
        xlog("Total number of $v edges: %lld\n", (long long)globals.sdbg_writer.num_tips());
    }

    // --- clean ---
    pthread_mutex_destroy(&globals.lv1_items_scanning_lock);
    free(globals.lv1_items);
#ifdef USE_GPU
    free(globals.lv2_substrings);
    free(globals.permutation);
    free(globals.lv2_substrings_db);
    free(globals.permutation_db);
    free_gpu_buffers(globals.gpu_key_buffer1, globals.gpu_key_buffer2, globals.gpu_value_buffer1, globals.gpu_value_buffer2);
#endif
}

} // s2

} // cx1_read2sdbg