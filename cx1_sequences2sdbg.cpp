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


#include "cx1_sequences2sdbg.h"

#include <omp.h>
#include <string>
#include <vector>

#include "utils.h"
#include "packed_reads.h"
#include "mac_pthread_barrier.h"
#include "kmer.h"

#ifndef USE_GPU
#include "lv2_cpu_sort.h"
#else
#include "lv2_gpu_functions.h"
#endif

namespace cx1_sequences2sdbg {

// helpers
typedef CX1<sequences2sdbg_global_t, kNumBuckets> cx1_t;
typedef CX1<sequences2sdbg_global_t, kNumBuckets>::readpartition_data_t readpartition_data_t;
typedef CX1<sequences2sdbg_global_t, kNumBuckets>::bucketpartition_data_t bucketpartition_data_t;
typedef CX1<sequences2sdbg_global_t, kNumBuckets>::outputpartition_data_t outputpartition_data_t;

/**
 * @brief encode seq_id and its offset in one int64_t
 */
inline int64_t EncodeEdgeOffset(int64_t seq_id, int offset, int strand, int bits_for_offset) {
    return (seq_id << (bits_for_offset + 1)) | (strand << bits_for_offset) | offset;
}

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

inline int ExtractFirstChar(uint32_t *item) {
    return *item >> kTopCharShift;
}

inline int Extract_a(uint32_t *item, int num_words, int64_t spacing, int kmer_k) {
    int non_dollar = (item[(num_words - 1) * spacing] >> (kBWTCharNumBits + kBitsPerMulti_t)) & 1;
    if (non_dollar) {
        int which_word = (kmer_k - 1) / kCharsPerEdgeWord;
        int word_index = (kmer_k - 1) % kCharsPerEdgeWord;
        return (item[which_word * spacing] >> (kCharsPerEdgeWord - 1 - word_index) * kBitsPerEdgeChar) & kEdgeCharMask;
    } else {
        return kSentinelValue;
    }
}

inline int Extract_b(uint32_t *item, int num_words, int64_t spacing) {
    return (item[(num_words - 1) * spacing] >> kBitsPerMulti_t) & ((1 << kBWTCharNumBits) - 1);
}

inline int ExtractCounting(uint32_t *item, int num_words, int64_t spacing) {
    return item[(num_words - 1) * spacing] & kMaxMulti_t;
}

// cx1 core functions
int64_t encode_lv1_diff_base(int64_t read_id, sequences2sdbg_global_t &g) {
    return EncodeEdgeOffset(read_id, 0, 0, g.bits_for_offset);
}

void read_seq_and_prepare(sequences2sdbg_global_t &globals) {
    // --- init reader ---
    SequenceManager seq_manager(&globals.package);
    seq_manager.set_multiplicity_vector(&globals.multiplicity);

    if (globals.contig_file_name != "") {
        seq_manager.set_file_type(SequenceManager::kMegahitContigs);
        seq_manager.set_file(globals.contig_file_name);
        seq_manager.ReadMegahitContigs(1LL << 60, 1LL << 60, true, true, globals.kmer_from, globals.kmer_k, false, true);
        seq_manager.clear();
    }

    if (globals.add_contig_file_name != "") {
        seq_manager.set_file_type(SequenceManager::kMegahitContigs);
        seq_manager.set_file(globals.add_contig_file_name);
        seq_manager.ReadMegahitContigs(1LL << 60, 1LL << 60, true, true, globals.kmer_from, globals.kmer_k, false, false);
        seq_manager.clear();
    }

    if (globals.input_prefix != "") {
        seq_manager.set_file_type(SequenceManager::kMegahitEdges);
        seq_manager.set_edge_files(globals.input_prefix + ".edges", 1);
        seq_manager.ReadEdges(1LL << 60, true);
        seq_manager.clear();
    }

    // --- compute bits_for_offset ---
    {
        globals.bits_for_offset = 0;
        int len = 1;
        while (len < (int)globals.package.max_read_len() + 1) {
            globals.bits_for_offset++;
            len *= 2;
        }
    }

    globals.num_seq = globals.package.size();

    // --- set cx1 param ---
    globals.cx1.num_cpu_threads_ = globals.num_cpu_threads;
    globals.cx1.num_output_threads_ = globals.num_output_threads;
    globals.cx1.num_items_ = globals.num_seq;
}

void* lv0_calc_bucket_size(void* _data) {
    readpartition_data_t &rp = *((readpartition_data_t*) _data);
    sequences2sdbg_global_t &globals = *(rp.globals);
    int64_t *bucket_sizes = rp.rp_bucket_sizes;
    memset(bucket_sizes, 0, kNumBuckets * sizeof(int64_t));

    for (int64_t seq_id = rp.rp_start_id; seq_id < rp.rp_end_id; ++seq_id) {
        int seq_len = globals.package.length(seq_id);
        if (seq_len < globals.kmer_k + 1) { continue; }

        uint32_t key = 0; // $$$$$$$$
        // build initial partial key
        for (int i = 0; i < kBucketPrefixLength - 1; ++i) {
            key = key * kBucketBase + globals.package.get_base(seq_id, i) + 1;
        }
        // sequence = xxxxxxxxx
        // edges = $xxxx, xxxxx, ..., xxxx$
        for (int i = kBucketPrefixLength - 1; i - (kBucketPrefixLength - 1) + globals.kmer_k - 1 <= seq_len; ++i) {
            key = (key * kBucketBase + globals.package.get_base(seq_id, i) + 1) % kNumBuckets;
            bucket_sizes[key]++;
        }

        // reverse complement
        key = 0;
        for (int i = 0; i < kBucketPrefixLength - 1; ++i) {
            key = key * kBucketBase + (3 - globals.package.get_base(seq_id, seq_len - 1 - i)) + 1; // complement
        }
        for (int i = kBucketPrefixLength - 1; i - (kBucketPrefixLength - 1) + globals.kmer_k - 1 <= seq_len; ++i) {
            key = key * kBucketBase + (3 - globals.package.get_base(seq_id, seq_len - 1 - i)) + 1;
            key %= kNumBuckets;
            bucket_sizes[key]++;
        }
    }
    return NULL;
}

void init_global_and_set_cx1(sequences2sdbg_global_t &globals) {
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
        err("[ERROR B::%s] Bucket too large for GPU: contains %lld items. Please try CPU version.\n", __func__, globals.max_bucket_size);
        // TODO: auto switch to CPU version
        exit(1);
    }
#endif
    globals.words_per_substring = DivCeiling(globals.kmer_k * kBitsPerEdgeChar + kBWTCharNumBits + 1 + kBitsPerMulti_t, kBitsPerEdgeWord);
    globals.words_per_dummy_node = DivCeiling(globals.kmer_k * kBitsPerEdgeChar, kBitsPerEdgeWord);
    // lv2 bytes: substring (double buffer), permutation, aux
    int64_t lv2_bytes_per_item = (globals.words_per_substring * sizeof(uint32_t) + sizeof(uint32_t)) * 2 + sizeof(unsigned char);
#ifndef USE_GPU
    lv2_bytes_per_item += sizeof(uint64_t) * 2; // simulate GPU
#endif

    if (cx1_t::kCX1Verbose >= 2) {
        log("[B::%s] %d words per substring, bits_for_offset: %d, words per dummy node ($v): %d\n", __func__, globals.words_per_substring, globals.bits_for_offset, globals.words_per_dummy_node);
    }

    // --- memory stuff ---
    globals.mem_packed_seq = globals.package.size_in_byte();
    int64_t mem_remained = globals.host_mem
                           - globals.mem_packed_seq
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
        log("[B::%s] Memory for sequence: %lld\n", __func__, globals.mem_packed_seq);
        log("[B::%s] max # lv.1 items = %lld\n", __func__, globals.cx1.max_lv1_items_);
        log("[B::%s] max # lv.2 items = %lld\n", __func__, globals.cx1.max_lv2_items_);
    }

    // --- alloc memory ---
    globals.lv1_items = (int*) MallocAndCheck(globals.cx1.max_lv1_items_ * sizeof(int), __FILE__, __LINE__);
    globals.lv2_substrings = (uint32_t*) MallocAndCheck(globals.cx1.max_lv2_items_ * globals.words_per_substring * sizeof(uint32_t), __FILE__, __LINE__);
    globals.permutation = (uint32_t *) MallocAndCheck(globals.cx1.max_lv2_items_ * sizeof(uint32_t), __FILE__, __LINE__);
    globals.lv2_substrings_db = (uint32_t*) MallocAndCheck(globals.cx1.max_lv2_items_ * globals.words_per_substring * sizeof(uint32_t), __FILE__, __LINE__);
    globals.permutation_db = (uint32_t *) MallocAndCheck(globals.cx1.max_lv2_items_ * sizeof(uint32_t), __FILE__, __LINE__);
    globals.lv2_aux = (unsigned char*) MallocAndCheck(globals.cx1.max_lv2_items_ * sizeof(unsigned char), __FILE__, __LINE__);
#ifndef USE_GPU
    globals.cpu_sort_space = (uint64_t*) MallocAndCheck(sizeof(uint64_t) * globals.cx1.max_lv2_items_, __FILE__, __LINE__);
#else
    alloc_gpu_buffers(globals.gpu_key_buffer1, globals.gpu_key_buffer2, globals.gpu_value_buffer1, globals.gpu_value_buffer2, (size_t)globals.cx1.max_lv2_items_);
#endif

    // --- init lock ---
    pthread_mutex_init(&globals.lv1_items_scanning_lock, NULL);

    // --- init stat ---
    globals.cur_prefix = -1;
    globals.cur_suffix_first_char = -1;
    globals.num_ones_in_last = 0;
    globals.total_number_edges = 0;
    globals.num_dollar_nodes = 0;
    memset(globals.num_chars_in_w, 0, sizeof(globals.num_chars_in_w));

    // --- init output ---
    globals.sdbg_writer.init((std::string(globals.output_prefix)+".w").c_str(),
                             (std::string(globals.output_prefix)+".last").c_str(),
                             (std::string(globals.output_prefix)+".isd").c_str());
    globals.dummy_nodes_writer.init((std::string(globals.output_prefix)+".dn").c_str());
    globals.output_f_file = OpenFileAndCheck((std::string(globals.output_prefix)+".f").c_str(), "w");
    globals.output_multiplicity_file = OpenFileAndCheck((std::string(globals.output_prefix)+".mul").c_str(), "wb");
    globals.output_multiplicity_file2 = OpenFileAndCheck((std::string(globals.output_prefix)+".mul2").c_str(), "wb");
    // --- write header ---
    fprintf(globals.output_f_file, "-1\n");
    globals.dummy_nodes_writer.output(globals.words_per_dummy_node);
}

void* lv1_fill_offset(void* _data) {
    readpartition_data_t &rp = *((readpartition_data_t*) _data);
    sequences2sdbg_global_t &globals = *(rp.globals);
    int64_t *prev_full_offsets = (int64_t *)MallocAndCheck(kNumBuckets * sizeof(int64_t), __FILE__, __LINE__); // temporary array for computing differentials
    assert(prev_full_offsets != NULL);
    for (int b = globals.cx1.lv1_start_bucket_; b < globals.cx1.lv1_end_bucket_; ++b)
        prev_full_offsets[b] = rp.rp_lv1_differential_base;
    // this loop is VERY similar to that in PreprocessScanToFillBucketSizesThread

    // ===== this is a macro to save some copy&paste ================
#define CHECK_AND_SAVE_OFFSET(offset, strand)                                                                   \
    do {                                                                                                        \
        if (((key - globals.cx1.lv1_start_bucket_) ^ (key - globals.cx1.lv1_end_bucket_)) & kSignBitMask) {     \
            int64_t full_offset = EncodeEdgeOffset(seq_id, offset, strand, globals.bits_for_offset);            \
            int64_t differential = full_offset - prev_full_offsets[key];                                        \
            if (differential > cx1_t::kDifferentialLimit) {                                                     \
                pthread_mutex_lock(&globals.lv1_items_scanning_lock);                                           \
                globals.lv1_items[rp.rp_bucket_offsets[key]++] = -globals.cx1.lv1_items_special_.size() - 1;    \
                globals.cx1.lv1_items_special_.push_back(full_offset);                                          \
                pthread_mutex_unlock(&globals.lv1_items_scanning_lock);                                         \
            } else {                                                                                            \
                assert(differential >= 0);                                                                      \
                globals.lv1_items[rp.rp_bucket_offsets[key]++] = (int) differential;                            \
            }                                                                                                   \
            prev_full_offsets[key] = full_offset;                                                               \
        }                                                                                                       \
    } while (0)
    // ^^^^^ why is the macro surrounded by a do-while? please ask Google
    // =========== end macro ==========================

    for (int64_t seq_id = rp.rp_start_id; seq_id < rp.rp_end_id; ++seq_id) {
        int seq_len = globals.package.length(seq_id);
        if (seq_len < globals.kmer_k + 1) { continue; }

        uint32_t key = 0; // $$$$$$$$
        // build initial partial key
        for (int i = 0; i < kBucketPrefixLength - 1; ++i) {
            key = key * kBucketBase + globals.package.get_base(seq_id, i) + 1;
        }
        // sequence = xxxxxxxxx
        // edges = $xxxx, xxxxx, ..., xxxx$
        for (int i = kBucketPrefixLength - 1; i - (kBucketPrefixLength - 1) + globals.kmer_k - 1 <= seq_len; ++i) {
            key = (key * kBucketBase + globals.package.get_base(seq_id, i) + 1) % kNumBuckets;
            CHECK_AND_SAVE_OFFSET(i - kBucketPrefixLength + 1, 0);
        }

        // reverse complement
        key = 0;
        for (int i = 0; i < kBucketPrefixLength - 1; ++i) {
            key = key * kBucketBase + (3 - globals.package.get_base(seq_id, seq_len - 1 - i)) + 1; // complement
        }
        for (int i = kBucketPrefixLength - 1; i - (kBucketPrefixLength - 1) + globals.kmer_k - 1 <= seq_len; ++i) {
            key = key * kBucketBase + (3 - globals.package.get_base(seq_id, seq_len - 1 - i)) + 1;
            key %= kNumBuckets;
            CHECK_AND_SAVE_OFFSET(i - kBucketPrefixLength + 1, 1);
        }
    }

#undef CHECK_AND_SAVE_OFFSET

    free(prev_full_offsets);
    return NULL;
}

// inline int BucketToPrefix(int x) {
//     int y = 0;
//     for (int i=0; i < kBucketPrefixLength; ++i) {
//         int z = x % kBucketBase;
//         if (z > 0) { --z; }
//         y |= (z << (i * kBitsPerEdgeChar));
//         x /= kBucketBase;
//     }
//     return y;
// }


void* lv2_extract_substr(void* _data) {
    bucketpartition_data_t &bp = *((bucketpartition_data_t*) _data);
    sequences2sdbg_global_t &globals = *(bp.globals);
    int *lv1_p = globals.lv1_items + globals.cx1.rp_[0].rp_bucket_offsets[bp.bp_start_bucket];
    int64_t offset_mask = (1 << globals.bits_for_offset) - 1; // 0000....00011..11
    uint32_t *substrings_p = globals.lv2_substrings +
                                (globals.cx1.rp_[0].rp_bucket_offsets[bp.bp_start_bucket] - globals.cx1.rp_[0].rp_bucket_offsets[globals.cx1.lv2_start_bucket_]);

    for (int bucket = bp.bp_start_bucket; bucket < bp.bp_end_bucket; ++bucket) {
        for (int t = 0; t < globals.num_cpu_threads; ++t) {
            int64_t full_offset = globals.cx1.rp_[t].rp_lv1_differential_base;
            int num = globals.cx1.rp_[t].rp_bucket_sizes[bucket];
            for (int i = 0; i < num; ++i) {
                if (*lv1_p >= 0) {
                    full_offset += *(lv1_p++);
                } else {
                    full_offset = globals.cx1.lv1_items_special_[-1 - *(lv1_p++)];
                }
                int64_t seq_id = full_offset >> (1 + globals.bits_for_offset);
                int offset = full_offset & offset_mask;
                int strand = (full_offset >> globals.bits_for_offset) & 1;
                int seq_len = globals.package.length(seq_id);
                int num_chars_to_copy = globals.kmer_k - (offset + globals.kmer_k > seq_len);
                int counting = 0;
                if (offset > 0 && offset + globals.kmer_k <= seq_len) {
                    counting = globals.multiplicity[seq_id];
                }

                int64_t which_word = globals.package.start_idx[seq_id] / 16;
                int start_offset = globals.package.start_idx[seq_id] % 16;
                int words_this_seq = DivCeiling(start_offset + seq_len, 16);
                uint32_t *edge_p = &globals.package.packed_seq[which_word];

                if (strand == 0) {
                    // copy counting and W char
                    int prev_char;
                    if (offset == 0) {
                        assert(num_chars_to_copy == globals.kmer_k);
                        prev_char = kSentinelValue;
                    } else {
                        prev_char = globals.package.get_base(seq_id, offset - 1);
                    }

                    CopySubstring(substrings_p, edge_p, offset + start_offset, num_chars_to_copy,
                                  globals.cx1.lv2_num_items_, words_this_seq, globals.words_per_substring);

                    uint32_t *last_word = substrings_p + int64_t(globals.words_per_substring - 1) * globals.cx1.lv2_num_items_;
                    *last_word |= int(num_chars_to_copy == globals.kmer_k) << (kBWTCharNumBits + kBitsPerMulti_t);
                    *last_word |= prev_char << kBitsPerMulti_t;
                    *last_word |= std::max(0, kMaxMulti_t - counting); // then larger counting come first after sorting
                } else {
                    int prev_char;
                    if (offset == 0) {
                        assert(num_chars_to_copy == globals.kmer_k);
                        prev_char = kSentinelValue;
                    } else {
                        prev_char = 3 - globals.package.get_base(seq_id, seq_len - 1 - offset + 1);
                    }

                    offset = seq_len - 1 - offset - (globals.kmer_k - 1); // switch to normal strand
                    if (offset < 0) { assert(num_chars_to_copy == globals.kmer_k - 1); offset = 0; }
                    CopySubstringRC(substrings_p, edge_p, offset + start_offset, num_chars_to_copy,
                                    globals.cx1.lv2_num_items_, words_this_seq, globals.words_per_substring);

                    uint32_t *last_word = substrings_p + int64_t(globals.words_per_substring - 1) * globals.cx1.lv2_num_items_;
                    *last_word |= int(num_chars_to_copy == globals.kmer_k) << (kBWTCharNumBits + kBitsPerMulti_t);
                    *last_word |= prev_char << kBitsPerMulti_t;
                    *last_word |= std::max(0, kMaxMulti_t - counting);
                }

                // if ((*substrings_p >> (32 - kBucketPrefixLength * 2)) != BucketToPrefix(bucket)) {
                //     err("WRONG substring wrong:%d right:%d read_id:%lld offset:%d strand: %d num_chars_to_copy:%d\n", *substrings_p >> (32 - kBucketPrefixLength * 2),  BucketToPrefix(bucket), seq_id, offset, strand, num_chars_to_copy);
                //     Kmer<6, uint32_t> kmer;
                //     kmer.init(edge_p, offset + start_offset, kBucketPrefixLength);
                //     err("%d\n", kmer.data_[0] >> (32 - kBucketPrefixLength * 2));
                // }
                substrings_p++;
            }
        }
    }
    return NULL;
}

void lv2_sort(sequences2sdbg_global_t &globals) {
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
        log("[B::%s] Sorting substrings with CPU...done. Time elapsed: %.4lf\n", __func__, local_timer.elapsed());
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
        log("[B::%s] Sorting substrings with GPU...done. Time elapsed: %.4lf\n", __func__, local_timer.elapsed());
    }
#endif
}

void lv2_pre_output_partition(sequences2sdbg_global_t &globals) {
    // swap double buffers
    globals.lv2_num_items_db = globals.cx1.lv2_num_items_;
    std::swap(globals.lv2_substrings_db, globals.lv2_substrings);
    std::swap(globals.permutation_db, globals.permutation);

    // distribute threads
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
                if (IsDiffKMinusOneMer(item, prev_item, globals.lv2_num_items_db, globals.kmer_k)) {
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

    memset(globals.lv2_aux, 0, sizeof(globals.lv2_aux[0]) * globals.lv2_num_items_db);
    pthread_barrier_init(&globals.output_barrier, NULL, globals.num_output_threads);
}

void* lv2_output(void* _op) {
    outputpartition_data_t *op = (outputpartition_data_t*) _op;
    sequences2sdbg_global_t &globals = *(op->globals);
    int64_t op_start_index = op->op_start_index;
    int64_t op_end_index = op->op_end_index;
    int start_idx, end_idx;
    int has_solid_a = 0; // has solid (k+1)-mer aSb
    int has_solid_b = 0; // has solid aSb
    int last_a[4], outputed_b;

    for (start_idx = op_start_index; start_idx < op_end_index; start_idx = end_idx) {
        end_idx = start_idx + 1;
        uint32_t *item = globals.lv2_substrings_db + globals.permutation_db[start_idx];
        while (end_idx < op_end_index &&
                !IsDiffKMinusOneMer(
                    item,
                    globals.lv2_substrings_db + globals.permutation_db[end_idx],
                    globals.lv2_num_items_db,
                    globals.kmer_k)) {
            ++end_idx;
        }

        // clean marking
        has_solid_a = has_solid_b = 0;
        outputed_b = 0;
        for (int i = start_idx; i < end_idx; ++i) {
            uint32_t *cur_item = globals.lv2_substrings_db + globals.permutation_db[i];
            int a = Extract_a(cur_item, globals.words_per_substring, globals.lv2_num_items_db, globals.kmer_k);
            int b = Extract_b(cur_item, globals.words_per_substring, globals.lv2_num_items_db);

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
            uint32_t *cur_item = globals.lv2_substrings_db + globals.permutation_db[i];
            int a = Extract_a(cur_item, globals.words_per_substring, globals.lv2_num_items_db, globals.kmer_k);
            int b = Extract_b(cur_item, globals.words_per_substring, globals.lv2_num_items_db);

            j = i + 1;
            while (j < end_idx) {
                uint32_t *next_item = globals.lv2_substrings_db + globals.permutation_db[j];
                if (Extract_a(next_item, globals.words_per_substring, globals.lv2_num_items_db, globals.kmer_k) != a ||
                        Extract_b(next_item, globals.words_per_substring, globals.lv2_num_items_db) != b) {
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
        for (int i = 0; i < globals.lv2_num_items_db; ++i) {
            if (globals.lv2_aux[i] & (1 << 7)) {
                uint32_t *item = globals.lv2_substrings_db + globals.permutation_db[i];
                while (ExtractFirstChar(item) > globals.cur_suffix_first_char) {
                    ++globals.cur_suffix_first_char;
                    fprintf(globals.output_f_file, "%lld\n", (long long)globals.total_number_edges);
                }

                multi_t counting_db = kMaxMulti_t - std::min(kMaxMulti_t,
                                               ExtractCounting(item, globals.words_per_substring, globals.lv2_num_items_db));
                // output
                globals.sdbg_writer.outputW(globals.lv2_aux[i] & 0xF);
                globals.sdbg_writer.outputLast((globals.lv2_aux[i] >> 4) & 1);
                globals.sdbg_writer.outputIsDollar((globals.lv2_aux[i] >> 5) & 1);
                if (counting_db <= kMaxMulti2_t) {
                    multi2_t c = counting_db;
                    fwrite(&c, sizeof(multi2_t), 1, globals.output_multiplicity_file);
                } else {
                    int64_t c = counting_db | (globals.total_number_edges << 16);
                    fwrite(&c, sizeof(int64_t), 1, globals.output_multiplicity_file2);
                    fwrite(&kMulti2Sp, sizeof(multi2_t), 1, globals.output_multiplicity_file);
                }

                globals.total_number_edges++;
                globals.num_chars_in_w[globals.lv2_aux[i] & 0xF]++;
                globals.num_ones_in_last += (globals.lv2_aux[i] >> 4) & 1;

                if ((globals.lv2_aux[i] >> 5) & 1) {
                    globals.num_dollar_nodes++;
                    if (globals.num_dollar_nodes >= kMaxDummyEdges) {
                        err("[ERROR B::%s] Too many dummy nodes (>= %lld)! The graph contains too many tips!\n", __func__, (long long)kMaxDummyEdges);
                        exit(1);
                    }
                    for (int64_t i = 0; i < globals.words_per_dummy_node; ++i) {
                        globals.dummy_nodes_writer.output(item[i * globals.lv2_num_items_db]);
                    }
                }
                if ((globals.lv2_aux[i] & 0xF) == 0) {
                    globals.num_dummy_edges++;
                }
            }
        }
        local_timer.stop();

        if (cx1_t::kCX1Verbose >= 4) {
            log("[B::%s] SdBG calc linear part: %lf\n", __func__, local_timer.elapsed());
        }
    }
    return NULL;
}

void lv2_post_output(sequences2sdbg_global_t &globals) {
    pthread_barrier_destroy(&globals.output_barrier);
}

void post_proc(sequences2sdbg_global_t &globals) {
    if (cx1_t::kCX1Verbose >= 2) {
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

    if (cx1_t::kCX1Verbose >= 2) {
        log("[B::%s] Total number of edges: %llu\n", __func__, globals.total_number_edges);
        log("[B::%s] Total number of ONEs: %llu\n", __func__, globals.num_ones_in_last);
        log("[B::%s] Total number of v$ edges: %llu\n", __func__, globals.num_dummy_edges);
        log("[B::%s] Total number of $v edges: %llu\n", __func__, globals.num_dollar_nodes);
    }

    // --- cleaning ---
    pthread_mutex_destroy(&globals.lv1_items_scanning_lock);
    free(globals.lv1_items);
    free(globals.lv2_substrings);
    free(globals.permutation);
    free(globals.lv2_substrings_db);
    free(globals.permutation_db);
    free(globals.lv2_aux);
    fclose(globals.output_f_file);
    fclose(globals.output_multiplicity_file);
    fclose(globals.output_multiplicity_file2);
    globals.dummy_nodes_writer.destroy();
#ifndef USE_GPU
    free(globals.cpu_sort_space);
#else
    free_gpu_buffers(globals.gpu_key_buffer1, globals.gpu_key_buffer2, globals.gpu_value_buffer1, globals.gpu_value_buffer2);
#endif
}

} // namespace