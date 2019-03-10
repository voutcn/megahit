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
#include <mutex>

#include "safe_alloc_open-inl.h"
#include "sequence/kseq.h"
#include "utils.h"
#include "sequence/kmer.h"
#include "packed_reads.h"
#include "read_lib_functions-inl.h"

#include "sorting.h"

namespace cx1_kmer_count {

// helpers
typedef CX1<count_global_t, kNumBuckets> cx1_t;
typedef CX1<count_global_t, kNumBuckets>::ReadPartition readpartition_data_t;

/**
 * @brief encode read_id and its offset in one int64_t
 */
inline int64_t EncodeOffset(int64_t read_id, int offset, int strand, SeqPackage &p) {
    return ((p.StartPos(read_id) + offset) << 1) | strand;
}

inline bool IsDifferentEdges(uint32_t *item1, uint32_t *item2, int num_words, int64_t spacing) {
    for (int i = num_words - 1; i >= 0; --i) {
        if (*(item1 + i * spacing) != *(item2 + i * spacing)) {
            return true;
        }
    }

    return false;
}

/**
 * @brief pack an edge and its multiplicity to word-aligned spaces
 */
inline void PackEdge(uint32_t *dest, uint32_t *item, int64_t counting, struct count_global_t &globals) {
    for (int i = 0; i < globals.words_per_edge && i < globals.words_per_substring; ++i) {
        dest[i] = *(item + i);
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

    dest[globals.words_per_edge - 1] |= std::min(int64_t(kMaxMul), counting);
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
        FILE *assist_seq_info = xfopen((globals.assist_seq_file + ".info").c_str(), "r");
        long long num_ass_bases, num_ass_seq;
        if (fscanf(assist_seq_info, "%lld%lld", &num_ass_seq, &num_ass_bases) != 2) {
            xfatal("Invalid format\n");
        }
        fclose(assist_seq_info);

        num_bases += num_ass_bases;
        num_reads += num_ass_seq;
    }

  globals.package.ReserveSequences(num_reads);
  globals.package.ReserveBases(num_bases);

    ReadBinaryLibs(globals.read_lib_file, globals.package, globals.lib_info, is_reverse);

    if (globals.assist_seq_file != "") {
        SequenceManager seq_manager;
        seq_manager.set_file_type(SequenceManager::kFastxReads);
        seq_manager.set_file(globals.assist_seq_file);

        bool reverse_read = true;
        bool append_to_package = true;
        bool trimN = false;

        seq_manager.ReadShortReads(1ULL << 60, 1ULL << 60, append_to_package, reverse_read, trimN);
        seq_manager.clear();
    }

    globals.package.BuildIndex();
    globals.max_read_length = globals.package.MaxSequenceLength();
    globals.num_reads = globals.package.size();

    xinfo("%ld reads, %d max read length\n", globals.num_reads, globals.max_read_length);

    globals.mem_packed_reads = globals.package.SizeInByte();

    int64_t mem_low_bound = globals.mem_packed_reads
                            + globals.num_reads * sizeof(unsigned char) * 2 // first_in0 & last_out0
                            + kNumBuckets * sizeof(int64_t) * (globals.num_cpu_threads * 3 + 1)
                            + (kMaxMul + 1) * (globals.num_output_threads + 1) * sizeof(int64_t);
    mem_low_bound *= 1.05;

    if (mem_low_bound > globals.host_mem) {
        xfatal("%lld bytes is not enough for CX1 sorting, please set -m parameter to at least %lld\n",
               globals.host_mem, mem_low_bound);
    }

    // set cx1 param
    globals.cx1.num_cpu_threads_ = globals.num_cpu_threads;
    globals.cx1.num_items_ = globals.num_reads;
}

void *lv0_calc_bucket_size(void *_data) {
    readpartition_data_t &rp = *((readpartition_data_t *) _data);
    count_global_t &globals = *(rp.globals);
    std::array<int64_t, kNumBuckets> &bucket_sizes = rp.rp_bucket_sizes;
    std::fill(bucket_sizes.begin(), bucket_sizes.end(), 0);
    GenericKmer edge, rev_edge; // (k+1)-mer and its rc

    for (int64_t read_id = rp.rp_start_id; read_id < rp.rp_end_id; ++read_id) {
        int read_length = globals.package.SequenceLength(read_id);

        if (read_length < globals.kmer_k + 1) {
            continue;
        }

        auto start_ptr_and_offset = globals.package.WordPtrAndOffset(read_id);
        edge.InitFromPtr(start_ptr_and_offset.first, start_ptr_and_offset.second, globals.kmer_k + 1);
        rev_edge = edge;
        rev_edge.ReverseComplement(globals.kmer_k + 1);

        int last_char_offset = globals.kmer_k;

        while (true) {
            if (rev_edge.cmp(edge, globals.kmer_k + 1) < 0) {
                bucket_sizes[rev_edge.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
            }
            else {
                bucket_sizes[edge.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
            }

            if (++last_char_offset >= read_length) {
                break;
            }
            else {
                int c = globals.package.GetBase(read_id, last_char_offset);
                edge.ShiftAppend(c, globals.kmer_k + 1);
                rev_edge.ShiftPreappend(3 - c, globals.kmer_k + 1);
            }
        }
    }
    return nullptr;
}

void init_global_and_set_cx1(count_global_t &globals) {
    // --- calculate lv2 memory ---
    globals.max_bucket_size = *std::max_element(globals.cx1.bucket_sizes_.begin(), globals.cx1.bucket_sizes_.end());
    globals.tot_bucket_size = 0;
    int num_non_empty = 0;

    for (int i = 0; i < kNumBuckets; ++i) {
        globals.tot_bucket_size += globals.cx1.bucket_sizes_[i];

        if (globals.cx1.bucket_sizes_[i] > 0) {
            num_non_empty++;
        }
    }

    globals.words_per_substring = DivCeiling((globals.kmer_k + 1) * kBitsPerEdgeChar, kBitsPerEdgeWord);
    globals.words_per_edge = DivCeiling((globals.kmer_k + 1) * kBitsPerEdgeChar + kBitsPerMul, kBitsPerEdgeWord);

    if (cx1_t::kCX1Verbose >= 2) {
        xinfo("%d words per substring, %d words per edge\n", globals.words_per_substring, globals.words_per_edge);
    }

    globals.cx1.lv1_just_go_ = true;
    globals.num_output_threads = globals.num_cpu_threads;

    num_non_empty = std::max(1, num_non_empty);

    for (int i = 0; i < kNumBuckets; ++i) {
        if (globals.cx1.bucket_sizes_[i] > 2 * globals.tot_bucket_size / num_non_empty) {
            // xinfo("Bucket %d size = %lld > %lld = 2 * avg\n", i, (long long)globals.cx1.bucket_sizes_[i], (long long)2 * globals.tot_bucket_size / num_non_empty);
        }
    }

    globals.max_sorting_items = std::max(3 * globals.tot_bucket_size / num_non_empty * globals.num_cpu_threads, globals.max_bucket_size);

    // lv2 bytes: substring, readinfo
    int64_t lv2_bytes_per_item = globals.words_per_substring * sizeof(uint32_t) + sizeof(int64_t);

    int64_t mem_remained = globals.host_mem
                           - globals.mem_packed_reads
                           - globals.num_cpu_threads * 65536 * sizeof(uint64_t) // radix sort buckets
                           - globals.num_reads * sizeof(unsigned char) * 2 // first_in0 & last_out0
                           - kNumBuckets * sizeof(int64_t) * (globals.num_cpu_threads * 3 + 1)
                           - (kMaxMul + 1) * (globals.num_output_threads + 1) * sizeof(int64_t);
    int64_t min_lv1_items = globals.tot_bucket_size / (kMaxLv1ScanTime - 0.5);


    if (globals.mem_flag == 1) {
        // auto set memory
        globals.cx1.max_lv1_items_ = int64_t(globals.tot_bucket_size / (kDefaultLv1ScanTime - 0.5));
        globals.cx1.max_lv1_items_ = std::max(globals.cx1.max_lv1_items_, globals.max_bucket_size);
        int64_t mem_needed = globals.cx1.max_lv1_items_ * cx1_t::kLv1BytePerItem + globals.max_sorting_items * lv2_bytes_per_item;
        xinfo("Set: %lld, %lld\n", mem_needed, mem_remained);
        if (mem_needed > mem_remained) {
            globals.cx1.adjust_mem_just_go(mem_remained, lv2_bytes_per_item, min_lv1_items, globals.max_bucket_size,
                                           globals.max_sorting_items, globals.cx1.max_lv1_items_, globals.max_sorting_items);
        }

    }
    else if (globals.mem_flag == 0) {
        // min memory
        globals.cx1.max_lv1_items_ = int64_t(globals.tot_bucket_size / (kMaxLv1ScanTime - 0.5));
        globals.cx1.max_lv1_items_ = std::max(globals.cx1.max_lv1_items_, globals.max_bucket_size);
        int64_t mem_needed = globals.cx1.max_lv1_items_ * cx1_t::kLv1BytePerItem + globals.max_sorting_items * lv2_bytes_per_item;

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
        xfatal("No enough memory to process.");
    }

    globals.cx1.max_mem_remain_ = globals.cx1.max_lv1_items_ * cx1_t::kLv1BytePerItem + globals.max_sorting_items * lv2_bytes_per_item;
    xinfo("%lld, %lld %lld %lld\n", globals.cx1.max_lv1_items_, globals.max_sorting_items, globals.cx1.max_mem_remain_, mem_remained);
    globals.cx1.bytes_per_sorting_item_ = lv2_bytes_per_item;
    globals.lv1_items = (int32_t *) xmalloc(
        globals.cx1.max_mem_remain_ + globals.num_cpu_threads * sizeof(uint64_t) * 65536);

    if (cx1_t::kCX1Verbose >= 2) {
        xinfo("Memory for reads: %lld\n", globals.mem_packed_reads);
        xinfo("max # lv.1 items = %lld\n", globals.cx1.max_lv1_items_);
    }

    // --- malloc read first_in / last_out ---
#ifdef LONG_READS
    globals.first_0_out = (uint32_t *) xmalloc(globals.num_reads * sizeof(uint32_t));
    globals.last_0_in = (uint32_t *) xmalloc(globals.num_reads * sizeof(uint32_t));
#else
    globals.first_0_out = (uint8_t *) MallocAndCheck(globals.num_reads * sizeof(uint8_t));
    globals.last_0_in = (uint8_t *) MallocAndCheck(globals.num_reads * sizeof(uint8_t));
#endif
    memset(globals.first_0_out, 0xFF, globals.num_reads * sizeof(globals.first_0_out[0]));
    memset(globals.last_0_in, 0xFF, globals.num_reads * sizeof(globals.last_0_in[0]));

    // --- initialize stat ---
    globals.edge_counting = (int64_t *) xmalloc((kMaxMul + 1) * sizeof(int64_t));
    globals.thread_edge_counting = (int64_t *) xmalloc((kMaxMul + 1) * globals.num_output_threads * sizeof(int64_t));
    memset(globals.edge_counting, 0, (kMaxMul + 1) * sizeof(int64_t));

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
    std::array<int64_t, kNumBuckets> prev_full_offsets;

    for (int b = globals.cx1.lv1_start_bucket_; b < globals.cx1.lv1_end_bucket_; ++b)
        prev_full_offsets[b] = rp.rp_lv1_differential_base;

    // this loop is VERY similar to that in PreprocessScanToFillBucketSizesThread
    GenericKmer edge, rev_edge; // (k+1)-mer and its rc
    int key;

    for (int64_t read_id = rp.rp_start_id; read_id < rp.rp_end_id; ++read_id) {
        int read_length = globals.package.SequenceLength(read_id);

        if (read_length < globals.kmer_k + 1) {
            continue;
        }

        auto ptr_and_offset = globals.package.WordPtrAndOffset(read_id);
        edge.InitFromPtr(ptr_and_offset.first, ptr_and_offset.second, globals.kmer_k + 1);
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
                std::lock_guard<std::mutex> lk(globals.lv1_items_scanning_lock);                                \
                globals.lv1_items[rp.rp_bucket_offsets[key_]++] = -globals.cx1.lv1_items_special_.size() - 1;   \
                globals.cx1.lv1_items_special_.push_back(full_offset);                                          \
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
                key = rev_edge.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
                CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 1);
            }
            else {
                key = edge.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
                CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 0);
            }

            if (++last_char_offset >= read_length) {
                break;
            }
            else {
                int c = globals.package.GetBase(read_id, last_char_offset);
                edge.ShiftAppend(c, globals.kmer_k + 1);
                rev_edge.ShiftPreappend(3 - c, globals.kmer_k + 1);
            }
        }
    }

#undef CHECK_AND_SAVE_OFFSET
    return nullptr;
}

inline void lv2_extract_substr_(uint32_t *substrings_p, count_global_t &globals, int start_bucket, int end_bucket) {
    int *lv1_p = globals.lv1_items + globals.cx1.rp_[0].rp_bucket_offsets[start_bucket];

    for (int b = start_bucket; b < end_bucket; ++b) {
        for (int t = 0; t < globals.num_cpu_threads; ++t) {
            int64_t full_offset = globals.cx1.rp_[t].rp_lv1_differential_base;
            int64_t num = globals.cx1.rp_[t].rp_bucket_sizes[b];

            for (int64_t i = 0; i < num; ++i) {
                if (*lv1_p >= 0) {
                    full_offset += *(lv1_p++);
                }
                else {
                    full_offset = globals.cx1.lv1_items_special_[-1 - * (lv1_p++)];
                }

                int64_t read_id = globals.package.GetSeqID(full_offset >> 1);
                int strand = full_offset & 1;
                int offset = (full_offset >> 1) - globals.package.StartPos(read_id);
                int num_chars_to_copy = globals.kmer_k + 1;

                int read_length = globals.package.SequenceLength(read_id);
                auto ptr_and_offset = globals.package.WordPtrAndOffset(read_id);
                int start_offset = ptr_and_offset.second;
                int words_this_seq = DivCeiling(start_offset + read_length, 16);
                const uint32_t *read_p = ptr_and_offset.first;

                unsigned char prev, next;

                if (offset > 0) {
                    prev = globals.package.GetBase(read_id, offset - 1);
                }
                else {
                    prev = kSentinelValue;
                }

                if (offset + globals.kmer_k + 1 < read_length) {
                    next = globals.package.GetBase(read_id, offset + globals.kmer_k + 1);
                }
                else {
                    next = kSentinelValue;
                }
                uint64_t *read_info_p = reinterpret_cast<uint64_t*>(substrings_p + globals.words_per_substring);
                if (strand == 0) {
                    CopySubstring(substrings_p, read_p, offset + start_offset, num_chars_to_copy,
                                  1, words_this_seq, globals.words_per_substring);
                    *read_info_p = (full_offset << 6) | (prev << 3) | next;
                }
                else {
                    CopySubstringRC(substrings_p, read_p, offset + start_offset, num_chars_to_copy,
                                    1, words_this_seq, globals.words_per_substring);
                    *read_info_p = (full_offset << 6) | ((next == kSentinelValue ? kSentinelValue : (3 - next)) << 3)
                                   | (prev == kSentinelValue ? kSentinelValue : (3 - prev));
                }

                substrings_p += globals.words_per_substring + 2;
            }
        }
    }
}

void lv2_output_(int64_t start_index, int64_t end_index, int thread_id, count_global_t &globals,
                 uint32_t *substrings) {
    uint32_t packed_edge[32];
    int64_t count_prev[5], count_next[5];
    int64_t *thread_edge_counting = globals.thread_edge_counting + thread_id * (kMaxMul + 1);

    int64_t from_;
    int64_t to_;

    EdgeWriter::Snapshot snapshot;

    for (int64_t i = start_index; i < end_index; i = to_) {
        from_ = i;
        to_ = i + 1;
        uint32_t *first_item = substrings + i * (globals.words_per_substring + 2);

        while (to_ < end_index) {
            if (IsDifferentEdges(first_item,
                                 substrings + to_ * (globals.words_per_substring + 2),
                                 globals.words_per_substring, 1)) {
                break;
            }

            ++to_;
        }

        int64_t count = to_ - from_;

        // update read's first and last

        memset(count_prev, 0, sizeof(count_prev[0]) * 4);
        memset(count_next, 0, sizeof(count_next[0]) * 4);
        bool has_in = false;
        bool has_out = false;

        for (int64_t j = from_; j < to_; ++j) {
          auto *read_info = reinterpret_cast<uint64_t*>(
                substrings + j * (globals.words_per_substring + 2) + globals.words_per_substring);
            int prev_and_next = *read_info & ((1 << 6) - 1);
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
            for (int64_t j = from_; j < to_; ++j) {
                auto *read_info_ptr = reinterpret_cast<uint64_t*>(
                  substrings + j * (globals.words_per_substring + 2) + globals.words_per_substring);
                int64_t read_info = *read_info_ptr >> 6;
                int64_t read_id = globals.package.GetSeqID(read_info >> 1);
                int strand = read_info & 1;
                uint32_t offset = (read_info >> 1) - globals.package.StartPos(read_id);

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
            for (int64_t j = from_; j < to_; ++j) {
                auto *read_info_ptr = reinterpret_cast<uint64_t*>(
                    substrings + j * (globals.words_per_substring + 2) + globals.words_per_substring);
                int64_t read_info = *read_info_ptr >> 6;
                int64_t read_id = globals.package.GetSeqID(read_info >> 1);
                int strand = read_info & 1;
                uint32_t offset = (read_info >> 1) - globals.package.StartPos(read_id);

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

        ++thread_edge_counting[std::min(count, int64_t(kMaxMul))];

        if (count >= globals.kmer_freq_threshold) {
            PackEdge(packed_edge, first_item, count, globals);
            globals.edge_writer.write(packed_edge, packed_edge[0] >> (32 - 2 * kBucketPrefixLength), thread_id,
                &snapshot);
        }
    }

    globals.edge_writer.SaveSnapshot(snapshot);
}

struct kt_sort_t {
    count_global_t *globals;
    std::vector<int64_t> thread_offset;
    std::vector<int64_t> rank;
    std::mutex mutex;
    int64_t acc = 0;
    int seen = 0;
};

void kt_sort(void *_g, long i, int tid) {
    kt_sort_t *kg = (kt_sort_t *)_g;
    int b = kg->globals->cx1.lv1_start_bucket_ + i;

    if (kg->thread_offset[tid] == -1) {
        std::lock_guard<std::mutex> lk(kg->mutex);
        kg->thread_offset[tid] = kg->acc;
        kg->acc += kg->globals->cx1.bucket_sizes_[b];
        kg->rank[tid] = kg->seen;
        kg->seen++;
    }

    if (kg->globals->cx1.bucket_sizes_[b] == 0) {
        return;
    }

    size_t offset = kg->globals->cx1.lv1_num_items_ * sizeof(int32_t) +
                    kg->thread_offset[tid] * kg->globals->cx1.bytes_per_sorting_item_ +
                    kg->rank[tid] * sizeof(uint64_t) * 65536;

    auto *substr_ptr = reinterpret_cast<uint32_t*>((char *)kg->globals->lv1_items + offset);
    lv2_extract_substr_(substr_ptr, *(kg->globals), b, b + 1);
    SortSubStr(substr_ptr, kg->globals->words_per_substring, kg->globals->cx1.bucket_sizes_[b], 2);
    lv2_output_(0, kg->globals->cx1.bucket_sizes_[b], tid, *(kg->globals), substr_ptr);
}

void lv1_direct_sort_and_count(count_global_t &globals) {
    kt_sort_t kg;
    kg.globals = &globals;

    kg.thread_offset.resize(globals.num_cpu_threads, -1);
    kg.rank.resize(globals.num_cpu_threads, 0);
    omp_set_num_threads(globals.num_cpu_threads);
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < globals.cx1.lv1_end_bucket_ - globals.cx1.lv1_start_bucket_; ++i) {
        kt_sort(&kg, i, omp_get_thread_num());
    }
}

void post_proc(count_global_t &globals) {
    for (int t = 0; t < globals.num_output_threads; ++t) {
        for (int i = 1; i <= kMaxMul; ++i) {
            globals.edge_counting[i] += globals.thread_edge_counting[t * (kMaxMul + 1) + i];
        }
    }

    // --- output reads for mercy ---
    int64_t num_candidate_reads = 0;
    int64_t num_has_tips = 0;
    FILE *candidate_file = xfopen((globals.output_prefix + ".cand").c_str(), "wb");
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
        xinfo("Total number of candidate reads: %lld(%lld)\n", num_candidate_reads, num_has_tips);
    }

    // --- stat ---
    int64_t num_solid_edges = 0;

    for (int i = globals.kmer_freq_threshold; i <= kMaxMul; ++i) {
        num_solid_edges += globals.edge_counting[i];
    }

    if (cx1_t::kCX1Verbose >= 2) {
        xinfo("Total number of solid edges: %llu\n", num_solid_edges);
    }

    FILE *counting_file = xfopen((globals.output_prefix + ".counting").c_str(), "w");

    for (int64_t i = 1, acc = 0; i <= kMaxMul; ++i) {
        acc += globals.edge_counting[i];
        fprintf(counting_file, "%lld %lld\n", (long long)i, (long long)acc);
    }

    fclose(counting_file);

    // --- cleaning ---
    free(globals.lv1_items);
    free(globals.first_0_out);
    free(globals.last_0_in);
    free(globals.edge_counting);
    free(globals.thread_edge_counting);
    globals.edge_writer.destroy();
}

} // namespace::cx1_kmer_count