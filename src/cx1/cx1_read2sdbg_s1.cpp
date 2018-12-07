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
#include "sequence/sequence_package.h"
#include "read_lib_functions-inl.h"

#include "sorting.h"
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
inline int64_t EncodeOffset(int64_t read_id, int offset, int strand, SeqPackage &p) {
    return ((p.StartPosition(read_id) + offset) << 1) | strand;
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

inline uint8_t ExtractPrevNext(int i, int64_t *readinfo) {
    return readinfo[i] & ((1 << 2 * kBWTCharNumBits) - 1);
}

// cx1 core functions

int64_t s1_encode_lv1_diff_base(int64_t read_id, read2sdbg_global_t &g) {
    return EncodeOffset(read_id, 0, 0, g.package);
}

void s1_read_input_prepare(read2sdbg_global_t &globals) {
    bool is_reverse = true;
    int64_t num_bases, num_reads;
    GetBinaryLibSize(globals.read_lib_file, num_bases, num_reads);

    globals.num_short_reads = num_reads;
    globals.num_short_read_bases = num_bases;

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
    // set up these figures before reading assist seq
    globals.max_read_length = globals.package.MaxSequenceLength();

    if (globals.assist_seq_file != "") {
        SequenceManager seq_manager;
        seq_manager.set_readlib_type(SequenceManager::kSingle);
        seq_manager.set_file_type(SequenceManager::kFastxReads);
        seq_manager.set_file(globals.assist_seq_file);
        seq_manager.set_package(&globals.package);

        bool reverse_read = true;
        bool append_to_package = true;
        bool trimN = false;

        seq_manager.ReadShortReads(1LL << 60, 1LL << 60, append_to_package, reverse_read, trimN);
        seq_manager.clear();
    }

  globals.package.BuildIndex();
    globals.num_reads = globals.package.size();

    xinfo("%ld reads, %d max read length, %lld total bases\n", globals.num_reads, globals.max_read_length,
          globals.package.BaseCount());

    int bits_read_length = 1; // bit needed to store read_length

    while ((1 << bits_read_length) - 1 < globals.max_read_length) {
        ++bits_read_length;
    }

    globals.read_length_mask = (1 << bits_read_length) - 1;
    globals.offset_num_bits = bits_read_length;

    // --- allocate memory for is_solid bit_vector
    globals.num_k1_per_read = globals.max_read_length - globals.kmer_k;

    if (globals.kmer_freq_threshold == 1) {
        // do not need to count solid kmers
        globals.mem_packed_reads = globals.package.SizeInByte();
    }
    else {
        globals.is_solid.reset(globals.package.BaseCount());
        globals.mem_packed_reads = DivCeiling(globals.package.BaseCount(), 8) + globals.package.SizeInByte();
    }

    int64_t mem_low_bound = globals.mem_packed_reads
                            + kNumBuckets * sizeof(int64_t) * (globals.num_cpu_threads * 3 + 1)
                            + (kMaxMul + 1) * (globals.num_output_threads + 1) * sizeof(int64_t);
    mem_low_bound *= 1.05;

    if (mem_low_bound > globals.host_mem) {
        xfatal("%lld bytes is not enough for CX1 sorting, please set -m parameter to at least %lld\n", globals.host_mem, mem_low_bound);
    }

    // set cx1 param
    globals.cx1.num_cpu_threads_ = globals.num_cpu_threads;
    globals.cx1.num_output_threads_ = globals.num_output_threads;
    globals.cx1.num_items_ = globals.num_reads;
}

void *s1_lv0_calc_bucket_size(void *_data) {
    readpartition_data_t &rp = *((readpartition_data_t *) _data);
    read2sdbg_global_t &globals = *(rp.globals);
    int64_t *bucket_sizes = rp.rp_bucket_sizes;
    memset(bucket_sizes, 0, kNumBuckets * sizeof(int64_t));
    GenericKmer k_minus1_mer, rev_k_minus1_mer; // (k-1)-mer and its rc

    for (int64_t read_id = rp.rp_start_id; read_id < rp.rp_end_id; ++read_id) {
        int read_length = globals.package.SequenceLength(read_id);

        if (read_length < globals.kmer_k + 1) {
            continue;
        }

        auto ptr_and_offset = globals.package.WordPtrAndOffset(read_id);
        int64_t offset = ptr_and_offset.second;
        const uint32_t *read_p = ptr_and_offset.first;

        k_minus1_mer.InitFromPtr(read_p, offset, globals.kmer_k - 1);
        rev_k_minus1_mer = k_minus1_mer;
        rev_k_minus1_mer.ReverseComplement(globals.kmer_k - 1);

        // the first one special handling
        bucket_sizes[k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
        bucket_sizes[rev_k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;

        int last_char_offset = globals.kmer_k - 1;
        int c = globals.package.GetBase(read_id, last_char_offset);
        k_minus1_mer.ShiftAppend(c, globals.kmer_k - 1);
        rev_k_minus1_mer.ShiftPreappend(3 - c, globals.kmer_k - 1);

        while (last_char_offset < read_length - 1) {
            int cmp = k_minus1_mer.cmp(rev_k_minus1_mer, globals.kmer_k - 1);

            if (cmp > 0) {
                bucket_sizes[rev_k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
            }
            else {
                bucket_sizes[k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
            }

            int c = globals.package.GetBase(read_id, ++last_char_offset);
            k_minus1_mer.ShiftAppend(c, globals.kmer_k - 1);
            rev_k_minus1_mer.ShiftPreappend(3 - c, globals.kmer_k - 1);
        }

        // last one special handling
        bucket_sizes[k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
        bucket_sizes[rev_k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
    }

    return NULL;
}

void s1_init_global_and_set_cx1(read2sdbg_global_t &globals) {
    globals.max_bucket_size = *std::max_element(globals.cx1.bucket_sizes_, globals.cx1.bucket_sizes_ + kNumBuckets);
    globals.tot_bucket_size = 0;
    int num_non_empty = 0;

    for (int i = 0; i < kNumBuckets; ++i) {
        globals.tot_bucket_size += globals.cx1.bucket_sizes_[i];
        num_non_empty += (globals.cx1.bucket_sizes_[i] > 0);
    }

    // to count (k+1)-mers, sort by the internal (k-1)-mer
    // (k+1)-mer = abS[0..k-2]cd
    // is solid: number of bSc >= threshold
    // bS has in coming: for some a, num of abS >= threshold
    // Sc has outgoing: for some a, num of Scd >= threshold
    globals.words_per_substring = DivCeiling((globals.kmer_k - 1) * kBitsPerEdgeChar + 2 * kBWTCharNumBits, kBitsPerEdgeWord);

    if (cx1_t::kCX1Verbose >= 2) {
        xinfo("%d words per substring\n", globals.words_per_substring);
    }

    // lv2 bytes: substring, permutation, readinfo
    int64_t lv2_bytes_per_item = (globals.words_per_substring) * sizeof(uint32_t) + sizeof(uint32_t) + sizeof(int64_t);

    globals.cx1.lv1_just_go_ = true;
    globals.num_output_threads = globals.num_cpu_threads;
    num_non_empty = std::max(1, num_non_empty);

    for (int i = 0; i < kNumBuckets; ++i) {
        if (globals.cx1.bucket_sizes_[i] > 2 * globals.tot_bucket_size / num_non_empty) {
            // xinfo("Bucket %d size = %lld > %lld = 2 * avg\n", i, (long long)globals.cx1.bucket_sizes_[i], (long long)2 * globals.tot_bucket_size / num_non_empty);
        }
    }

    globals.max_sorting_items = std::max(3 * globals.tot_bucket_size / num_non_empty * globals.num_cpu_threads, globals.max_bucket_size);

    lv2_bytes_per_item += sizeof(uint32_t); // CPU memory is used to simulate GPU

    int64_t mem_remained = globals.host_mem
                           - globals.mem_packed_reads
                           - globals.num_cpu_threads * 65536 * sizeof(uint64_t) // radix sort buckets
                           - kNumBuckets * sizeof(int64_t) * (globals.num_cpu_threads * 3 + 1)
                           - (kMaxMul + 1) * (globals.num_output_threads + 1) * sizeof(int64_t);
    int64_t min_lv1_items = globals.tot_bucket_size / (kMaxLv1ScanTime - 0.5);

    if (globals.mem_flag == 1) {
        // auto set memory
        globals.cx1.max_lv1_items_ = int64_t(globals.tot_bucket_size / (kDefaultLv1ScanTime - 0.5));
        globals.cx1.max_lv1_items_ = std::max(globals.cx1.max_lv1_items_, globals.max_bucket_size);
        int64_t mem_needed = globals.cx1.max_lv1_items_ * cx1_t::kLv1BytePerItem + globals.max_sorting_items * lv2_bytes_per_item;

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

    globals.cx1.max_mem_remain_ = globals.cx1.max_lv1_items_ * sizeof(int) + globals.max_sorting_items * lv2_bytes_per_item;
    globals.cx1.bytes_per_sorting_item_ = lv2_bytes_per_item;

    globals.lv1_items = (int *) xmalloc(
        globals.cx1.max_mem_remain_ + globals.num_cpu_threads * sizeof(uint64_t) * 65536, __FILE__, __LINE__);

    if (cx1_t::kCX1Verbose >= 2) {
        xinfo("Memory for reads: %lld\n", globals.mem_packed_reads);
        xinfo("max # lv.1 items = %lld\n", globals.cx1.max_lv1_items_);
    }

    // --- initialize output mercy files ---
    globals.num_mercy_files = 1;

    while (globals.num_mercy_files * 10485760LL < globals.num_short_reads && globals.num_mercy_files < 64) {
        globals.num_mercy_files <<= 1;
    }

    if (cx1_t::kCX1Verbose >= 3) {
        xinfo("Number of files for mercy candidate reads: %d\n", globals.num_mercy_files);
    }

    for (int i = 0; i < globals.num_mercy_files; ++i) {
        globals.mercy_files.push_back(xfopen(FormatString("%s.mercy_cand.%d", globals.output_prefix.c_str(), i), "wb"));
    }

    pthread_mutex_init(&globals.lv1_items_scanning_lock, NULL); // init lock

    // --- initialize stat ---
    globals.edge_counting = (int64_t *) xmalloc((kMaxMul + 1) * sizeof(int64_t), __FILE__, __LINE__);
    globals.thread_edge_counting = (int64_t *) xmalloc((kMaxMul + 1) * globals.num_output_threads * sizeof(int64_t),
                                                       __FILE__,
                                                       __LINE__);
    memset(globals.edge_counting, 0, (kMaxMul + 1) * sizeof(int64_t));
    memset(globals.thread_edge_counting, 0, sizeof(int64_t) * (kMaxMul + 1) * globals.num_output_threads);
}

void *s1_lv1_fill_offset(void *_data) {
    readpartition_data_t &rp = *((readpartition_data_t *) _data);
    read2sdbg_global_t &globals = *(rp.globals);
    int64_t *prev_full_offsets = (int64_t *) xmalloc(kNumBuckets * sizeof(int64_t), __FILE__, __LINE__); // temporary array for computing differentials
    assert(prev_full_offsets != NULL);

    for (int b = globals.cx1.lv1_start_bucket_; b < globals.cx1.lv1_end_bucket_; ++b)
        prev_full_offsets[b] = rp.rp_lv1_differential_base;

    // this loop is VERY similar to that in PreprocessScanToFillBucketSizesThread
    GenericKmer k_minus1_mer, rev_k_minus1_mer; // (k+1)-mer and its rc

    int key;

    for (int64_t read_id = rp.rp_start_id; read_id < rp.rp_end_id; ++read_id) {
        int read_length = globals.package.SequenceLength(read_id);

        if (read_length < globals.kmer_k + 1) {
            continue;
        }

        auto ptr_and_offset = globals.package.WordPtrAndOffset(read_id);
        int64_t offset = ptr_and_offset.second;
        const uint32_t *read_p = ptr_and_offset.first;

        k_minus1_mer.InitFromPtr(read_p, offset, globals.kmer_k - 1);
        rev_k_minus1_mer = k_minus1_mer;
        rev_k_minus1_mer.ReverseComplement(globals.kmer_k - 1);

        // ===== this is a macro to save some copy&paste ================
#define CHECK_AND_SAVE_OFFSET(offset, strand)                                                                  \
    do {                                                                                                       \
        if (globals.cx1.cur_lv1_buckets_[key]) {                                                               \
            int key_ = globals.cx1.bucket_rank_[key];                                                          \
            int64_t full_offset = EncodeOffset(read_id, offset, strand, globals.package);                      \
            int64_t differential = full_offset - prev_full_offsets[key_];                                      \
            if (differential > cx1_t::kDifferentialLimit) {                                                    \
                pthread_mutex_lock(&globals.lv1_items_scanning_lock);                                          \
                globals.lv1_items[rp.rp_bucket_offsets[key_]++] = -globals.cx1.lv1_items_special_.size() - 1;  \
                globals.cx1.lv1_items_special_.push_back(full_offset);                                         \
                pthread_mutex_unlock(&globals.lv1_items_scanning_lock);                                        \
            } else {                                                                                           \
                assert((int) differential >= 0);                                                               \
                globals.lv1_items[rp.rp_bucket_offsets[key_]++] = (int) differential;                          \
            }                                                                                                  \
            prev_full_offsets[key_] = full_offset;                                                             \
        }                                                                                                      \
    } while (0)
        // ^^^^^ why is the macro surrounded by a do-while? please ask Google
        // =========== end macro ==========================

        // the first one special handling
        key = k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
        CHECK_AND_SAVE_OFFSET(0, 0);
        key = rev_k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
        CHECK_AND_SAVE_OFFSET(0, 1);

        int last_char_offset = globals.kmer_k - 1;
        int c = globals.package.GetBase(read_id, last_char_offset);
        k_minus1_mer.ShiftAppend(c, globals.kmer_k - 1);
        rev_k_minus1_mer.ShiftPreappend(3 - c, globals.kmer_k - 1);

        // shift the key char by char
        while (last_char_offset < read_length - 1) {
            int cmp = k_minus1_mer.cmp(rev_k_minus1_mer, globals.kmer_k - 1);

            if (cmp > 0) {
                key = rev_k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
                CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k + 2, 1);
            }
            else if (cmp < 0) {
                key = k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
                CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k + 2, 0);
            }
            else {
                // a not-that-math-correct solution if the edge is palindrome, but works well enough
                int prev = globals.package.GetBase(read_id, last_char_offset - (globals.kmer_k - 1));
                int next = globals.package.GetBase(read_id, last_char_offset + 1);

                if (prev <= 3 - next) {
                    key = rev_k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
                    CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k + 2, 0);
                }
                else {
                    key = rev_k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
                    CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k + 2, 1);
                }
            }

            int c = globals.package.GetBase(read_id, ++last_char_offset);
            k_minus1_mer.ShiftAppend(c, globals.kmer_k - 1);
            rev_k_minus1_mer.ShiftPreappend(3 - c, globals.kmer_k - 1);
        }

        // the last one special handling
        key = k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
        CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k + 2, 0);
        key = rev_k_minus1_mer.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
        CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k + 2, 1);
    }

#undef CHECK_AND_SAVE_OFFSET

    free(prev_full_offsets);
    return NULL;
}

void s1_extract_subtstr_(int bp_from, int bp_to, read2sdbg_global_t &globals, uint32_t *substr, int64_t *readinfo_ptr) {
    int *lv1_p = globals.lv1_items + globals.cx1.rp_[0].rp_bucket_offsets[bp_from];

    for (int b = bp_from; b < bp_to; ++b) {
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
                int offset = (full_offset >> 1) - globals.package.StartPosition(read_id);
                int read_length = globals.package.SequenceLength(read_id);
                int num_chars_to_copy = globals.kmer_k - 1;
                unsigned char prev, next, head, tail; // (k+1)=abScd, prev=a, head=b, tail=c, next=d

                assert(offset < read_length);

                if (offset > 1) {
                    head = globals.package.GetBase(read_id, offset - 1);
                    prev = globals.package.GetBase(read_id, offset - 2);
                }
                else {
                    prev = kSentinelValue;

                    if (offset > 0) {
                        head = globals.package.GetBase(read_id, offset - 1);
                    }
                    else {
                        head = kSentinelValue;
                    }
                }

                if (offset + globals.kmer_k < read_length) {
                    tail = globals.package.GetBase(read_id, offset + globals.kmer_k - 1);
                    next = globals.package.GetBase(read_id, offset + globals.kmer_k);
                }
                else {
                    next = kSentinelValue;

                    if (offset + globals.kmer_k - 1 < read_length) {
                        tail = globals.package.GetBase(read_id, offset + globals.kmer_k - 1);
                    }
                    else {
                        tail = kSentinelValue;
                    }
                }

                auto ptr_and_offset = globals.package.WordPtrAndOffset(read_id);
                int start_offset = ptr_and_offset.second;
                const uint32_t *read_p = ptr_and_offset.first;
                int words_this_read = DivCeiling(start_offset + read_length, 16);

                if (strand == 0) {
                    CopySubstring(substr, read_p, offset + start_offset, num_chars_to_copy,
                                  1, words_this_read, globals.words_per_substring);
                    uint32_t *last_word = substr + int64_t(globals.words_per_substring - 1) * 1;
                    *last_word |= (head << kBWTCharNumBits) | tail;
                    *readinfo_ptr = (full_offset << 6) | (prev << 3) | next;
                }
                else {
                    CopySubstringRC(substr, read_p, offset + start_offset, num_chars_to_copy,
                                    1, words_this_read, globals.words_per_substring);
                    uint32_t *last_word = substr + int64_t(globals.words_per_substring - 1) * 1;
                    *last_word |= ((tail == kSentinelValue ? kSentinelValue : 3 - tail) << kBWTCharNumBits) | (head == kSentinelValue ? kSentinelValue : 3 - head);
                    *readinfo_ptr = (full_offset << 6) | ((next == kSentinelValue ? kSentinelValue : (3 - next)) << 3)
                                    | (prev == kSentinelValue ? kSentinelValue : (3 - prev));
                }

                substr += globals.words_per_substring;
                readinfo_ptr++;
            }
        }
    }
}

void s1_lv2_output_(int64_t from, int64_t to, int tid, read2sdbg_global_t &globals, uint32_t *substr, int64_t *readinfo_ptr) {
    int64_t end_idx;
    int64_t count_prev_head[5][5];
    int64_t count_tail_next[5][5];
    int64_t count_head_tail[(1 << 2 * kBWTCharNumBits) - 1];
    int64_t *thread_edge_counting = globals.thread_edge_counting + tid * (kMaxMul + 1);

    for (int64_t i = from; i < to; i = end_idx) {
        end_idx = i + 1;
        uint32_t *first_item = substr + i * globals.words_per_substring;
        memset(count_prev_head, 0, sizeof(count_prev_head));
        memset(count_tail_next, 0, sizeof(count_tail_next));
        memset(count_head_tail, 0, sizeof(count_head_tail));

        {
            uint8_t prev_and_next = ExtractPrevNext(i, readinfo_ptr);
            uint8_t head_and_tail = ExtractHeadTail(first_item, 1, globals.words_per_substring);
            count_prev_head[prev_and_next >> 3][head_and_tail >> 3]++;
            count_tail_next[head_and_tail & 7][prev_and_next & 7]++;
            count_head_tail[head_and_tail]++;
        }

        while (end_idx < to) {
            if (IsDiffKMinusOneMer(first_item,
                                   substr + end_idx * globals.words_per_substring,
                                   1,
                                   globals.kmer_k)) {
                break;
            }

            uint8_t prev_and_next = ExtractPrevNext(end_idx, readinfo_ptr);
            uint8_t head_and_tail = ExtractHeadTail(substr + end_idx * globals.words_per_substring, 1, globals.words_per_substring);
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
            uint8_t head_and_tail = ExtractHeadTail(substr + i * globals.words_per_substring, 1, globals.words_per_substring);
            uint8_t head = head_and_tail >> 3;
            uint8_t tail = head_and_tail & 7;

            if (head != kSentinelValue && tail != kSentinelValue) {
                ++thread_edge_counting[std::min(int64_t(kMaxMul), count_head_tail[head_and_tail])];
            }

            if (head != kSentinelValue && tail != kSentinelValue && count_head_tail[head_and_tail] >= globals.kmer_freq_threshold) {
                for (int64_t j = 0; j < count_head_tail[head_and_tail]; ++j, ++i) {
                    int64_t read_info = readinfo_ptr[i] >> 6;
                    int strand = read_info & 1;
                    int64_t read_id = globals.package.GetSeqID(read_info >> 1);
                    int offset = (read_info >> 1) - globals.package.StartPosition(read_id) - 1;
                    int l_offset = strand == 0 ? offset : offset + 1;
                    int r_offset = strand == 0 ? offset + 1 : offset;

                    // mark this is a solid edge
                    globals.is_solid.set((read_info >> 1) - 1);

                    if (!(has_in & (1 << head))) {
                        // no in
                        int64_t packed_mercy_cand = ((globals.package.StartPosition(read_id) + l_offset) << 2) | (1 + strand);
                        fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1, globals.mercy_files[read_id & (globals.num_mercy_files - 1)]);
                    }

                    if (!(has_out & (1 << tail))) {
                        // no out
                        int64_t packed_mercy_cand = ((globals.package.StartPosition(read_id) + r_offset) << 2) | (2 - strand);
                        fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1, globals.mercy_files[read_id & (globals.num_mercy_files - 1)]);
                    }
                }
            }
            else {
                // not solid, but we still need to tell whether its left/right kmer is solid
                for (int64_t j = 0; j < count_head_tail[head_and_tail]; ++j, ++i) {
                    int64_t read_info = readinfo_ptr[i] >> 6;
                    int strand = read_info & 1;
                    int64_t read_id = globals.package.GetSeqID(read_info >> 1);
                    int offset = (read_info >> 1) - globals.package.StartPosition(read_id) - 1;
                    int l_offset = strand == 0 ? offset : offset + 1;
                    int r_offset = strand == 0 ? offset + 1 : offset;

                    if (l_has_out & (1 << head)) {
                        if (has_in & (1 << head)) {
                            // has both in & out
                            int64_t packed_mercy_cand = ((globals.package.StartPosition(read_id) + l_offset) << 2) | 0;
                            fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1, globals.mercy_files[read_id & (globals.num_mercy_files - 1)]);
                        }
                        else {
                            // has out but no in
                            int64_t packed_mercy_cand = ((globals.package.StartPosition(read_id) + l_offset) << 2) | (1 + strand);
                            fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1, globals.mercy_files[read_id & (globals.num_mercy_files - 1)]);
                        }
                    }
                    else {
                        if (has_in & (1 << head)) {
                            // has in but no out
                            int64_t packed_mercy_cand = ((globals.package.StartPosition(read_id) + l_offset) << 2) | (2 - strand);
                            fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1, globals.mercy_files[read_id & (globals.num_mercy_files - 1)]);
                        }
                    }

                    if (r_has_in & (1 << tail)) {
                        if (has_out & (1 << tail)) {
                            // has both in & out
                            int64_t packed_mercy_cand = ((globals.package.StartPosition(read_id) + r_offset) << 2) | 0;
                            fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1, globals.mercy_files[read_id & (globals.num_mercy_files - 1)]);
                        }
                        else {
                            // has in but no out
                            int64_t packed_mercy_cand = ((globals.package.StartPosition(read_id) + r_offset) << 2) | (2 - strand);
                            fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1, globals.mercy_files[read_id & (globals.num_mercy_files - 1)]);
                        }
                    }
                    else {
                        if (has_out & (1 << tail)) {
                            // has out but no in
                            int64_t packed_mercy_cand = ((globals.package.StartPosition(read_id) + r_offset) << 2) | (1 + strand);
                            fwrite(&packed_mercy_cand, sizeof(packed_mercy_cand), 1, globals.mercy_files[read_id & (globals.num_mercy_files - 1)]);
                        }
                    }
                }
            }
        }
    }
}

struct kt_sort_t {
    read2sdbg_global_t *globals;
    std::vector<int64_t> thread_offset;
    std::vector<int> rank;
    int64_t acc = 0;
    int seen = 0;
    std::mutex mutex;
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

    uint32_t *substr_ptr = (uint32_t *) ((char *)kg->globals->lv1_items + offset);
    uint64_t *bucket = (uint64_t *)(substr_ptr + kg->globals->cx1.bucket_sizes_[b] * kg->globals->words_per_substring);
    uint32_t *permutation_ptr = (uint32_t *)(bucket + 65536);
    uint32_t *cpu_sort_space_ptr = permutation_ptr + kg->globals->cx1.bucket_sizes_[b];
    int64_t *readinfo_ptr = (int64_t *) (cpu_sort_space_ptr + kg->globals->cx1.bucket_sizes_[b]);

    s1_extract_subtstr_(b, b + 1, *(kg->globals), substr_ptr, readinfo_ptr);
    lv2_cpu_radix_sort_st2(substr_ptr, kg->globals->words_per_substring, kg->globals->cx1.bucket_sizes_[b]);
    s1_lv2_output_(0, kg->globals->cx1.bucket_sizes_[b], tid, *(kg->globals), substr_ptr, readinfo_ptr);
}

void s1_lv1_direct_sort_and_count(read2sdbg_global_t &globals) {
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

void s1_post_proc(read2sdbg_global_t &globals) {
    for (int t = 0; t < globals.num_output_threads; ++t) {
        for (int i = 1; i <= kMaxMul; ++i) {
            globals.edge_counting[i] += globals.thread_edge_counting[t * (kMaxMul + 1) + i];
        }
    }

    // --- stat ---
    int64_t num_solid_edges = 0;

    for (int i = globals.kmer_freq_threshold; i <= kMaxMul; ++i) {
        num_solid_edges += globals.edge_counting[i];
    }

    if (cx1_t::kCX1Verbose >= 2) {
        xinfo("Total number of solid edges: %llu\n", num_solid_edges);
    }

    FILE *counting_file = xfopen((std::string(globals.output_prefix) + ".counting").c_str(), "w");

    for (int64_t i = 1, acc = 0; i <= kMaxMul; ++i) {
        acc += globals.edge_counting[i];
        fprintf(counting_file, "%lld %lld\n", (long long)i, (long long)acc);
    }

    fclose(counting_file);

    // --- cleaning ---
    pthread_mutex_destroy(&globals.lv1_items_scanning_lock);
    free(globals.lv1_items);

    for (int i = 0; i < globals.num_mercy_files; ++i) {
        fclose(globals.mercy_files[i]);
    }
}

} // s1

} // cx1_read2sdbg