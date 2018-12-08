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


#include "cx1_seq2sdbg.h"

#include <omp.h>
#include <string>
#include <vector>
#include <mutex>

#include "utils.h"
#include "packed_reads.h"
#include "sequence/kmer.h"
#include "iterate/async_sequence_reader.h"

#include "sorting.h"

namespace cx1_seq2sdbg {

// helpers
typedef CX1<seq2sdbg_global_t, kNumBuckets> cx1_t;
typedef CX1<seq2sdbg_global_t, kNumBuckets>::readpartition_data_t readpartition_data_t;

/**
 * @brief encode seq_id and its offset in one int64_t
 */
inline int64_t EncodeEdgeOffset(int64_t seq_id, int offset, int strand, SeqPackage &p) {
    return ((p.StartPosition(seq_id) + offset) << 1) | strand;
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

inline int Extract_a(uint32_t *item, int num_words, int64_t spacing, int kmer_k) {
    int non_dollar = (item[(num_words - 1) * spacing] >> (kBWTCharNumBits + kBitsPerMul)) & 1;

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
    return (item[(num_words - 1) * spacing] >> kBitsPerMul) & ((1 << kBWTCharNumBits) - 1);
}

inline int ExtractCounting(uint32_t *item, int num_words, int64_t spacing) {
    return item[(num_words - 1) * spacing] & kMaxMul;
}

// cx1 core functions
int64_t encode_lv1_diff_base(int64_t read_id, seq2sdbg_global_t &g) {
    assert(read_id < (int64_t)g.package.size());
    return EncodeEdgeOffset(read_id, 0, 0, g.package);
}

/**
 * @brief build lkt for faster binary search for mercy
 */
void InitLookupTable(int64_t *lookup_table, SeqPackage &p) {
    memset(lookup_table, 0xFF, sizeof(int64_t) * kLookUpSize * 2);

    if (p.size() == 0) {
        return;
    }

    Kmer<1, uint32_t> kmer;
    kmer.InitFromPtr(p.WordPtrAndOffset(0).first, 0, 16);

    uint32_t cur_prefix = kmer.data()[0] >> kLookUpShift;
    lookup_table[cur_prefix * 2] = 0;

    for (int64_t i = 1, num_edges = p.size(); i < num_edges; ++i) {
        auto ptr_and_offset = p.WordPtrAndOffset(i);
        kmer.InitFromPtr(ptr_and_offset.first, ptr_and_offset.second, 16);

        if ((kmer.data()[0] >> kLookUpShift) > cur_prefix) {
            lookup_table[cur_prefix * 2 + 1] = i - 1;
            cur_prefix = kmer.data()[0] >> kLookUpShift;
            lookup_table[cur_prefix * 2] = i;
        }
        else {
            assert(cur_prefix == (kmer.data()[0] >> kLookUpShift));
        }
    }

    lookup_table[cur_prefix * 2 + 1] = p.size() - 1;
}

/**
 * @brief search mercy kmer
 */
int64_t BinarySearchKmer(GenericKmer &kmer, int64_t *lookup_table, SeqPackage &p, int kmer_size) {
    // --- first look up ---
    int64_t l = lookup_table[(kmer.data()[0] >> kLookUpShift) * 2];

    if (l == -1) {
        return -1;
    }

    int64_t r = lookup_table[(kmer.data()[0] >> kLookUpShift) * 2 + 1];
    GenericKmer mid_kmer;

    while (l <= r) {
        int64_t mid = (l + r) / 2;
        auto ptr_and_offset = p.WordPtrAndOffset(mid);
        mid_kmer.InitFromPtr(ptr_and_offset.first, ptr_and_offset.second, kmer_size);
        int cmp = kmer.cmp(mid_kmer, kmer_size);

        if (cmp > 0) {
            l = mid + 1;
        }
        else if (cmp < 0) {
            r = mid - 1;
        }
        else {
            return mid;
        }
    }

    return -1;
}

static void *MercyInputThread(void *seq_manager) {
    SequenceManager *sm = (SequenceManager *)seq_manager;
    int64_t kMaxReads = 1 << 22;
    int64_t kMaxBases = 1 << 28;
    bool append = false;
    bool reverse = false;
    sm->ReadShortReads(kMaxReads, kMaxBases, append, reverse);
    return nullptr;
}

void GenMercyEdges(seq2sdbg_global_t &globals) {
    int64_t *edge_lookup = (int64_t *) xmalloc(kLookUpSize * 2 * sizeof(int64_t), __FILE__, __LINE__);
    InitLookupTable(edge_lookup, globals.package);

    std::vector<GenericKmer> mercy_edges;
    AsyncReadReader reader(globals.input_prefix + ".cand");

    int num_threads = globals.num_cpu_threads - 1;
    omp_set_num_threads(num_threads);
    omp_lock_t mercy_lock;
    omp_init_lock(&mercy_lock);

    int64_t num_mercy_edges = 0;
    int64_t num_mercy_reads = 0;

    while (true) {
        SeqPackage &rp = reader.Next();

        if (rp.size() == 0) {
            break;
        }

        num_mercy_reads += rp.size();
        mercy_edges.clear();
        #pragma omp parallel for reduction(+:num_mercy_edges)
        for (unsigned read_id = 0; read_id < rp.size(); ++read_id) {

            int read_len = rp.SequenceLength(read_id);

            if (read_len < globals.kmer_k + 2) {
                continue;
            }

            std::vector<bool> has_in, has_out;
            GenericKmer kmer, rev_kmer;

            has_in.resize(read_len);
            has_out.resize(read_len);
            std::fill(has_in.begin(), has_in.end(), false);
            std::fill(has_out.begin(), has_out.end(), false);

            auto ptr_and_offset = rp.WordPtrAndOffset(read_id);
            kmer.InitFromPtr(ptr_and_offset.first, ptr_and_offset.second, globals.kmer_k);
            rev_kmer = kmer;
            rev_kmer.ReverseComplement(globals.kmer_k);

            // mark those positions with in/out
            for (int i = 0; i + globals.kmer_k <= read_len; ++i) {
                if (!has_in[i]) {
                    // search rc
                    if (BinarySearchKmer(rev_kmer, edge_lookup, globals.package, globals.kmer_k) != -1) {
                        has_in[i] = true;
                    }
                    else {
                        // left append ACGT to kmer, if the (k+1)-mer exist, the kmer has in
                        rev_kmer.SetBase(globals.kmer_k, 3); // rev kmer is used to compare to kmer, if it's smaller, kmer would not exist in the table
                        kmer.ShiftPreappend(0, globals.kmer_k + 1);

                        for (int c = 0; c < 4; ++c) {
                            kmer.SetBase(0, c);

                            if (kmer.cmp(rev_kmer, globals.kmer_k + 1) > 0) {
                                break;
                            }

                            if (BinarySearchKmer(kmer, edge_lookup, globals.package, globals.kmer_k + 1) != -1) {
                                has_in[i] = true;
                                break;
                            }
                        }

                        rev_kmer.SetBase(globals.kmer_k, 0);
                        kmer.ShiftAppend(0, globals.kmer_k + 1); // clean the k+1-th char
                    }
                }

                // check whether has out
                int64_t edge_id = BinarySearchKmer(kmer, edge_lookup, globals.package, globals.kmer_k);

                if (edge_id != -1) {
                    has_out[i] = true;

                    // BWT see whether the next has in too
                    if (i + globals.kmer_k < read_len &&
                        globals.package.GetBase(edge_id, globals.kmer_k) == rp.GetBase(read_id, i + globals.kmer_k)) {
                        has_in[i + 1] = true;
                    }
                }
                else {
                    // search the rc
                    kmer.SetBase(globals.kmer_k, 3);
                    int next_char = i + globals.kmer_k < read_len ? 3 - rp.GetBase(read_id, i + globals.kmer_k) : 0;
                    rev_kmer.ShiftPreappend(next_char, globals.kmer_k + 1);

                    if (rev_kmer.cmp(kmer, globals.kmer_k + 1) <= 0 && BinarySearchKmer(rev_kmer, edge_lookup, globals.package, globals.kmer_k + 1) != -1) {
                        has_out[i] = true;
                        has_in[i + 1] = true;
                    }
                    else {
                        for (int c = 0; c < 4; ++c) {
                            if (c == next_char) {
                                continue;
                            }

                            rev_kmer.SetBase(0, c);

                            if (rev_kmer.cmp(kmer, globals.kmer_k + 1) > 0) {
                                break;
                            }

                            if (BinarySearchKmer(rev_kmer, edge_lookup, globals.package, globals.kmer_k + 1) != -1) {
                                has_out[i] = true;
                                break;
                            }
                        }
                    }

                    kmer.SetBase(globals.kmer_k, 0);
                    rev_kmer.ShiftAppend(0, globals.kmer_k + 1);
                }

                // shift kmer and rev_kmer
                if (i + globals.kmer_k < read_len) {
                    int next_char = rp.GetBase(read_id, i + globals.kmer_k);
                    kmer.ShiftAppend(next_char, globals.kmer_k);
                    rev_kmer.ShiftPreappend(3 - next_char, globals.kmer_k);
                }
            }

            // adding mercy edges
            int last_no_out = -1;

            for (int i = 0; i + globals.kmer_k <= read_len; ++i) {
                switch (has_in[i] | (int(has_out[i]) << 1)) {
                case 1: { // has incoming only
                    last_no_out = i;
                    break;
                }

                case 2: { // has outgoing only
                    if (last_no_out >= 0) {
                        for (int j = last_no_out; j < i; ++j) {
                            omp_set_lock(&mercy_lock);
                            auto ptr_and_offset = rp.WordPtrAndOffset(read_id);
                            mercy_edges.push_back(GenericKmer(ptr_and_offset.first, ptr_and_offset.second + j, globals.kmer_k + 1));
                            omp_unset_lock(&mercy_lock);
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
        }

        for (unsigned i = 0; i < mercy_edges.size(); ++i) {
            globals.package.AppendCompactSequence(mercy_edges[i].data(), globals.kmer_k + 1);
        }
    }

    omp_destroy_lock(&mercy_lock);
    free(edge_lookup);

    globals.multiplicity.insert(globals.multiplicity.end(), num_mercy_edges, 1);

    if (cx1_t::kCX1Verbose >= 2) {
        xinfo("Number of reads: %ld, Number of mercy edges: %ld\n", num_mercy_reads, num_mercy_edges);
    }
}

void read_seq_and_prepare(seq2sdbg_global_t &globals) {
    // --- init reader ---
    SequenceManager seq_manager(&globals.package);
    seq_manager.set_multiplicity_vector(&globals.multiplicity);

    // reserve space
    {
        long long bases_to_reserve = 0;
        long long num_contigs_to_reserve = 0;
        long long num_multiplicities_to_reserve = 0;

        if (globals.input_prefix != "") {
            EdgeReader edge_reader;
            edge_reader.set_file_prefix(globals.input_prefix);
            edge_reader.read_info();
            int64_t num_edges = edge_reader.num_edges();
            xinfo("Number edges: %lld\n", (long long)num_edges);

            if (globals.need_mercy) {
                num_edges *= 1.25;    // it is rare that # mercy > 25%
            }

            bases_to_reserve += num_edges * (edge_reader.kmer_size() + 1);
            num_multiplicities_to_reserve += num_edges;
        }

        if (globals.contig != "") {
            long long num_contigs, num_bases;
            FILE *contig_info = xfopen((globals.contig + ".info").c_str(), "r");
            if (fscanf(contig_info, "%lld%lld", &num_contigs, &num_bases) != 2) {
              xfatal("Invalid format\n");
            }
            bases_to_reserve += num_bases;
            num_contigs_to_reserve += num_contigs;
            num_multiplicities_to_reserve += num_contigs;
            fclose(contig_info);
        }

        if (globals.addi_contig != "") {
            long long num_contigs, num_bases;
            FILE *contig_info = xfopen((globals.addi_contig + ".info").c_str(), "r");
            if (fscanf(contig_info, "%lld%lld", &num_contigs, &num_bases) != 2) {
                xfatal("Invalid format\n");
            }
            bases_to_reserve += num_bases;
            num_contigs_to_reserve += num_contigs;
            num_multiplicities_to_reserve += num_contigs;
            fclose(contig_info);
        }

        if (globals.local_contig != "") {
            long long num_contigs, num_bases;
            FILE *contig_info = xfopen((globals.local_contig + ".info").c_str(), "r");
            if (fscanf(contig_info, "%lld%lld", &num_contigs, &num_bases) != 2) {
                xfatal("Invalid format\n");
            }
            bases_to_reserve += num_bases;
            num_contigs_to_reserve += num_contigs;
            num_multiplicities_to_reserve += num_contigs;
            fclose(contig_info);
        }

        xinfo("Bases to reserve: %lld, number contigs: %lld, number multiplicity: %lld\n", bases_to_reserve, num_contigs_to_reserve, num_multiplicities_to_reserve);
        globals.package.ReserveSequences(num_contigs_to_reserve);
        globals.package.ReserveBases(bases_to_reserve);
        globals.multiplicity.reserve(num_multiplicities_to_reserve);
    }

    xinfo("Before reading, sizeof seq_package: %lld, multiplicity vector: %lld\n",
          globals.package.SizeInByte(), globals.multiplicity.capacity());

    if (globals.input_prefix != "") {
        seq_manager.set_file_type(globals.need_mercy ? SequenceManager::kSortedEdges : SequenceManager::kMegahitEdges);
        seq_manager.set_edge_files(globals.input_prefix);
        seq_manager.ReadEdgesWithFixedLen(1LL << 60, true);
        seq_manager.clear();
    }

    if (globals.need_mercy) {
        SimpleTimer timer;

        if (cx1_t::kCX1Verbose >= 3) {
            timer.reset();
            timer.start();
            xinfo("Adding mercy edges...\n");
        }

        GenMercyEdges(globals);

        if (cx1_t::kCX1Verbose >= 3) {
            timer.stop();
            xinfo("Done. Time elapsed: %.4lf\n", timer.elapsed());
        }
    }

    if (globals.contig != "") {
        seq_manager.set_file_type(SequenceManager::kMegahitContigs);
        seq_manager.set_file(globals.contig);
        seq_manager.set_kmer_size(globals.kmer_from, globals.kmer_k);
        seq_manager.set_min_len(globals.kmer_k + 1);

        bool contig_reverse = true;
        bool append_to_package = true;
        int discard_flag = 0;
        bool extend_loop = true;
        bool calc_depth = false;

        seq_manager.ReadMegahitContigs(1LL << 60, 1LL << 60, append_to_package, contig_reverse, discard_flag, extend_loop, calc_depth);
        seq_manager.clear();

        // read bubble
        seq_manager.set_file_type(SequenceManager::kMegahitContigs);
        seq_manager.set_file(globals.bubble_seq);
        seq_manager.set_kmer_size(globals.kmer_from, globals.kmer_k);
        seq_manager.set_min_len(globals.kmer_k + 1);

        contig_reverse = true;
        append_to_package = true;
        discard_flag = 0;
        extend_loop = true;
        calc_depth = false;

        seq_manager.ReadMegahitContigs(1LL << 60, 1LL << 60, append_to_package, contig_reverse, discard_flag, extend_loop, calc_depth);
        seq_manager.clear();
    }

    if (globals.addi_contig != "") {
        seq_manager.set_file_type(SequenceManager::kMegahitContigs);
        seq_manager.set_file(globals.addi_contig);
        seq_manager.set_kmer_size(globals.kmer_from, globals.kmer_k);
        seq_manager.set_min_len(globals.kmer_k + 1);

        bool contig_reverse = true;
        bool append_to_package = true;
        int discard_flag = 0;
        bool extend_loop = true;
        bool calc_depth = false;

        seq_manager.ReadMegahitContigs(1LL << 60, 1LL << 60, append_to_package, contig_reverse, discard_flag, extend_loop, calc_depth);
        seq_manager.clear();
    }

    if (globals.local_contig != "") {
        seq_manager.set_file_type(SequenceManager::kMegahitContigs);
        seq_manager.set_file(globals.local_contig);
        seq_manager.set_kmer_size(globals.kmer_from, globals.kmer_k);
        seq_manager.set_min_len(globals.kmer_k + 1);

        bool contig_reverse = true;
        bool append_to_package = true;
        int discard_flag = 0;
        bool extend_loop = true;
        bool calc_depth = false;

        seq_manager.ReadMegahitContigs(1LL << 60, 1LL << 60, append_to_package, contig_reverse, discard_flag, extend_loop, calc_depth);
        seq_manager.clear();
    }

    xinfo("After reading, sizeof seq_package: %lld, multiplicity vector: %lld\n",
          globals.package.SizeInByte(), globals.multiplicity.capacity());

    globals.package.BuildIndex();
    globals.num_seq = globals.package.size();

    globals.mem_packed_seq = globals.package.SizeInByte() + globals.multiplicity.size() * sizeof(mul_t);
    int64_t mem_low_bound = globals.mem_packed_seq
                            + kNumBuckets * sizeof(int64_t) * (globals.num_cpu_threads * 3 + 1);
    mem_low_bound *= 1.05;

    if (mem_low_bound > globals.host_mem) {
        xfatal("%lld bytes is not enough for CX1 sorting, please set -m parameter to at least %lld\n", globals.host_mem, mem_low_bound);
    }

    // --- set cx1 param ---
    globals.cx1.num_cpu_threads_ = globals.num_cpu_threads;
    globals.cx1.num_items_ = globals.num_seq;
}

void *lv0_calc_bucket_size(void *_data) {
    readpartition_data_t &rp = *((readpartition_data_t *) _data);
    seq2sdbg_global_t &globals = *(rp.globals);
    int64_t *bucket_sizes = rp.rp_bucket_sizes;
    memset(bucket_sizes, 0, kNumBuckets * sizeof(int64_t));

    for (int64_t seq_id = rp.rp_start_id; seq_id < rp.rp_end_id; ++seq_id) {
        int seq_len = globals.package.SequenceLength(seq_id);

        if (seq_len < globals.kmer_k + 1) {
            continue;
        }

        uint32_t key = 0; // $$$$$$$$

        // build initial partial key
        for (int i = 0; i < kBucketPrefixLength - 1; ++i) {
            key = key * kBucketBase + globals.package.GetBase(seq_id, i);
        }

        // sequence = xxxxxxxxx
        // edges = $xxxx, xxxxx, ..., xxxx$
        for (int i = kBucketPrefixLength - 1; i - (kBucketPrefixLength - 1) + globals.kmer_k - 1 <= seq_len; ++i) {
            key = (key * kBucketBase + globals.package.GetBase(seq_id, i)) % kNumBuckets;
            bucket_sizes[key]++;
        }

        // reverse complement
        key = 0;

        for (int i = 0; i < kBucketPrefixLength - 1; ++i) {
            key = key * kBucketBase + (3 - globals.package.GetBase(seq_id, seq_len - 1 - i)); // complement
        }

        for (int i = kBucketPrefixLength - 1; i - (kBucketPrefixLength - 1) + globals.kmer_k - 1 <= seq_len; ++i) {
            key = key * kBucketBase + (3 - globals.package.GetBase(seq_id, seq_len - 1 - i));
            key %= kNumBuckets;
            bucket_sizes[key]++;
        }
    }
    return nullptr;
}

void init_global_and_set_cx1(seq2sdbg_global_t &globals) {
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

    globals.words_per_substring = DivCeiling(globals.kmer_k * kBitsPerEdgeChar + kBWTCharNumBits + 1 + kBitsPerMul, kBitsPerEdgeWord);
    globals.words_per_dummy_node = DivCeiling(globals.kmer_k * kBitsPerEdgeChar, kBitsPerEdgeWord);

    num_non_empty = std::max(1, num_non_empty);

    for (int i = 0; i < kNumBuckets; ++i) {
        if (globals.cx1.bucket_sizes_[i] > 2 * globals.tot_bucket_size / num_non_empty) {
            // xinfo("Bucket %d size = %lld > %lld = 2 * avg\n", i, (long long)globals.cx1.bucket_sizes_[i], (long long)2 * globals.tot_bucket_size / num_non_empty);
        }
    }

    int64_t lv2_bytes_per_item = globals.words_per_substring * sizeof(uint32_t);

    globals.max_sorting_items = std::max(3 * globals.tot_bucket_size * globals.num_cpu_threads / num_non_empty, globals.max_bucket_size);
    globals.cx1.lv1_just_go_ = true;
    globals.num_output_threads = globals.num_cpu_threads;

    int64_t mem_remained = globals.host_mem
                           - globals.mem_packed_seq
                           - globals.num_cpu_threads * 65536 * sizeof(uint64_t) // radix sort buckets
                           - kNumBuckets * sizeof(int64_t) * (globals.num_cpu_threads * 3 + 1);
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

    globals.lv1_items = (int32_t *) xmalloc(
        globals.cx1.max_mem_remain_ + globals.num_cpu_threads * sizeof(uint64_t) * 65536, __FILE__, __LINE__);

    if (cx1_t::kCX1Verbose >= 2) {
        xinfo("Memory for sequence: %lld\n", globals.mem_packed_seq);
        xinfo("max # lv.1 items = %lld\n", globals.cx1.max_lv1_items_);
    }

    // --- init output ---
    globals.sdbg_writer.set_num_threads(globals.num_output_threads);
    globals.sdbg_writer.set_kmer_size(globals.kmer_k);
    globals.sdbg_writer.set_num_buckets(kNumBuckets);
    globals.sdbg_writer.set_file_prefix(globals.output_prefix);
    globals.sdbg_writer.InitFiles();
}

void *lv1_fill_offset(void *_data) {
    readpartition_data_t &rp = *((readpartition_data_t *) _data);
    seq2sdbg_global_t &globals = *(rp.globals);
    std::array<int64_t, kNumBuckets> prev_full_offsets;

    for (int b = globals.cx1.lv1_start_bucket_; b < globals.cx1.lv1_end_bucket_; ++b)
        prev_full_offsets[b] = rp.rp_lv1_differential_base;

    // this loop is VERY similar to that in PreprocessScanToFillBucketSizesThread

    // ===== this is a macro to save some copy&paste ================
#define CHECK_AND_SAVE_OFFSET(key__, offset, strand)                                                            \
    do {                                                                                                        \
        if (globals.cx1.cur_lv1_buckets_[key__]) {                                                              \
            int key_ = globals.cx1.bucket_rank_[key__];                                                         \
            int64_t full_offset = EncodeEdgeOffset(seq_id, offset, strand, globals.package);                    \
            int64_t differential = full_offset - prev_full_offsets[key_];                                       \
            if (differential > cx1_t::kDifferentialLimit) {                                                     \
                std::lock_guard<std::mutex> lk(globals.lv1_items_scanning_lock);                                \
                globals.lv1_items[rp.rp_bucket_offsets[key_]++] = -globals.cx1.lv1_items_special_.size() - 1;   \
                globals.cx1.lv1_items_special_.push_back(full_offset);                                          \
            } else {                                                                                            \
                assert(differential >= 0);                                                                      \
                globals.lv1_items[rp.rp_bucket_offsets[key_]++] = (int) differential;                           \
            }                                                                                                   \
            assert(rp.rp_bucket_offsets[key_] <= globals.cx1.lv1_num_items_);                                   \
            prev_full_offsets[key_] = full_offset;                                                              \
        }                                                                                                       \
    } while (0)
    // ^^^^^ why is the macro surrounded by a do-while? please ask Google
    // =========== end macro ==========================

    for (int64_t seq_id = rp.rp_start_id; seq_id < rp.rp_end_id; ++seq_id) {
        int seq_len = globals.package.SequenceLength(seq_id);

        if (seq_len < globals.kmer_k + 1) {
            continue;
        }

        // build initial partial key
        Kmer<1, uint32_t> kmer, rev_kmer;
        auto ptr_and_offset = globals.package.WordPtrAndOffset(seq_id);
        kmer.InitFromPtr(ptr_and_offset.first, ptr_and_offset.second, kBucketPrefixLength);
        auto rev_ptr_and_offset = globals.package.WordPtrAndOffset(seq_id, seq_len - kBucketPrefixLength);
        rev_kmer.InitFromPtr(rev_ptr_and_offset.first, rev_ptr_and_offset.second, kBucketPrefixLength);
        rev_kmer.ReverseComplement(kBucketPrefixLength);


        int key = kmer.data()[0] >> (32 - kBucketPrefixLength * 2);
        int rev_key = rev_kmer.data()[0] >> (32 - kBucketPrefixLength * 2);
        CHECK_AND_SAVE_OFFSET(key, 0, 0);
        CHECK_AND_SAVE_OFFSET(rev_key, 0, 1);

        // sequence = xxxxxxxxx
        // edges = $xxxx, xxxxx, ..., xxxx$
        for (int i = kBucketPrefixLength; i - (kBucketPrefixLength - 1) + globals.kmer_k - 1 <= seq_len; ++i) {
            key = (key * kBucketBase + globals.package.GetBase(seq_id, i)) % kNumBuckets;
            rev_key = rev_key * kBucketBase + (3 - globals.package.GetBase(seq_id, seq_len - 1 - i));
            rev_key %= kNumBuckets;
            CHECK_AND_SAVE_OFFSET(key, i - kBucketPrefixLength + 1, 0);
            CHECK_AND_SAVE_OFFSET(rev_key, i - kBucketPrefixLength + 1, 1);
        }
    }

#undef CHECK_AND_SAVE_OFFSET
    return nullptr;
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

void lv2_extract_substr_(int from_bucket, int to_bucket, seq2sdbg_global_t &globals, uint32_t *substr) {
    int *lv1_p = globals.lv1_items + globals.cx1.rp_[0].rp_bucket_offsets[from_bucket];

    for (int bucket = from_bucket; bucket < to_bucket; ++bucket) {
        for (int t = 0; t < globals.num_cpu_threads; ++t) {
            int64_t full_offset = globals.cx1.rp_[t].rp_lv1_differential_base;
            int64_t num = globals.cx1.rp_[t].rp_bucket_sizes[bucket];

            for (int64_t i = 0; i < num; ++i) {
                if (*lv1_p >= 0) {
                    full_offset += *(lv1_p++);
                }
                else {
                    full_offset = globals.cx1.lv1_items_special_[-1 - * (lv1_p++)];
                }

                int64_t seq_id = globals.package.GetSeqID(full_offset >> 1);
                int offset = (full_offset >> 1) - globals.package.StartPosition(seq_id);
                int strand = full_offset & 1;

                int seq_len = globals.package.SequenceLength(seq_id);
                int num_chars_to_copy = globals.kmer_k - (offset + globals.kmer_k > seq_len);
                int counting = 0;

                if (offset > 0 && offset + globals.kmer_k <= seq_len) {
                    counting = globals.multiplicity[seq_id];
                }

                auto ptr_and_offset = globals.package.WordPtrAndOffset(seq_id);
                int start_offset = ptr_and_offset.second;
                int words_this_seq = DivCeiling(start_offset + seq_len, 16);
                const uint32_t *edge_p = ptr_and_offset.first;

                if (strand == 0) {
                    // copy counting and W char
                    int prev_char;

                    if (offset == 0) {
                        assert(num_chars_to_copy == globals.kmer_k);
                        prev_char = kSentinelValue;
                    }
                    else {
                        prev_char = globals.package.GetBase(seq_id, offset - 1);
                    }

                    CopySubstring(substr, edge_p, offset + start_offset, num_chars_to_copy,
                                  1, words_this_seq, globals.words_per_substring);

                    uint32_t *last_word = substr + globals.words_per_substring - 1;
                    *last_word |= int(num_chars_to_copy == globals.kmer_k) << (kBWTCharNumBits + kBitsPerMul);
                    *last_word |= prev_char << kBitsPerMul;
                    *last_word |= std::max(0, kMaxMul - counting); // then larger counting come first after sorting
                }
                else {
                    int prev_char;

                    if (offset == 0) {
                        assert(num_chars_to_copy == globals.kmer_k);
                        prev_char = kSentinelValue;
                    }
                    else {
                        prev_char = 3 - globals.package.GetBase(seq_id, seq_len - 1 - offset + 1);
                    }

                    offset = seq_len - 1 - offset - (globals.kmer_k - 1); // switch to normal strand

                    if (offset < 0) {
                        assert(num_chars_to_copy == globals.kmer_k - 1);
                        offset = 0;
                    }

                    CopySubstringRC(substr, edge_p, offset + start_offset, num_chars_to_copy,
                                    1, words_this_seq, globals.words_per_substring);

                    uint32_t *last_word = substr + globals.words_per_substring - 1;
                    *last_word |= int(num_chars_to_copy == globals.kmer_k) << (kBWTCharNumBits + kBitsPerMul);
                    *last_word |= prev_char << kBitsPerMul;
                    *last_word |= std::max(0, kMaxMul - counting);
                }

                substr += globals.words_per_substring;
            }
        }
    }
}

void output_(int64_t from, int64_t to, seq2sdbg_global_t &globals, uint32_t *substr, int tid) {
    int64_t start_idx, end_idx;
    int has_solid_a = 0; // has solid (k+1)-mer aSb
    int has_solid_b = 0; // has solid aSb
    int64_t last_a[4], outputed_b;
    uint32_t tip_label[32];

    for (start_idx = from; start_idx < to; start_idx = end_idx) {
        end_idx = start_idx + 1;
        uint32_t *item = substr + start_idx * globals.words_per_substring;

        while (end_idx < to &&
                !IsDiffKMinusOneMer(
                    item,
                    substr + end_idx * globals.words_per_substring,
                    1,
                    globals.kmer_k)) {
            ++end_idx;
        }

        // clean marking
        has_solid_a = has_solid_b = 0;
        outputed_b = 0;

        for (int64_t i = start_idx; i < end_idx; ++i) {
            uint32_t *cur_item = substr + i * globals.words_per_substring;
            int a = Extract_a(cur_item, globals.words_per_substring, 1, globals.kmer_k);
            int b = Extract_b(cur_item, globals.words_per_substring, 1);

            if (a != kSentinelValue && b != kSentinelValue) {
                has_solid_a |= 1 << a;
                has_solid_b |= 1 << b;
            }

            if (a != kSentinelValue &&
                    (b != kSentinelValue || !(has_solid_a & (1 << a)))) {
                last_a[a] = i;
            }
        }

        for (int64_t i = start_idx, j; i < end_idx; i = j) {
            uint32_t *cur_item = substr + i * globals.words_per_substring;
            int a = Extract_a(cur_item, globals.words_per_substring, 1, globals.kmer_k);
            int b = Extract_b(cur_item, globals.words_per_substring, 1);

            j = i + 1;

            while (j < end_idx) {
                uint32_t *next_item = substr + j * globals.words_per_substring;

                if (Extract_a(next_item, globals.words_per_substring, 1, globals.kmer_k) != a ||
                        Extract_b(next_item, globals.words_per_substring, 1) != b) {
                    break;
                }
                else {
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
            last = (a == kSentinelValue) ? 0 : ((last_a[a] == j - 1) ? 1 : 0);
            outputed_b |= 1 << b;

            if (is_dollar) {
                for (int64_t i = 0; i < globals.words_per_dummy_node; ++i) {
                    tip_label[i] = cur_item[i];
                }
            }

            globals.sdbg_writer.Write(tid,
                                      cur_item[0] >> (32 - kBucketPrefixLength * 2),
                                      w,
                                      last,
                                      is_dollar,
                                      kMaxMul - ExtractCounting(cur_item, globals.words_per_substring, 1),
                                      tip_label);
        }
    }
}

struct kt_sort_t {
    seq2sdbg_global_t *globals;
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
    lv2_extract_substr_(b, b + 1, *(kg->globals), substr_ptr);
    SortSubStr(substr_ptr, kg->globals->words_per_substring, kg->globals->cx1.bucket_sizes_[b]);
    output_(0, kg->globals->cx1.bucket_sizes_[b], *(kg->globals), substr_ptr, tid);
}

void lv1_direct_sort_and_proc(seq2sdbg_global_t &globals) {
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

void post_proc(seq2sdbg_global_t &globals) {
    globals.sdbg_writer.Finalize();
    if (cx1_t::kCX1Verbose >= 2) {
        xinfo("Number of $ A C G T A- C- G- T-:\n");
    }
    xinfo("");
    for (int i = 0; i < 9; ++i) {
        xinfoc("%lld ", (long long)globals.sdbg_writer.final_meta().w_count(i));
    }

    xinfoc("\n");

    if (cx1_t::kCX1Verbose >= 2) {
        xinfo("Total number of edges: %lld\n", (long long) globals.sdbg_writer.final_meta().item_count());
        xinfo("Total number of ONEs: %lld\n", (long long)globals.sdbg_writer.final_meta().ones_in_last());
        xinfo("Total number of $v edges: %lld\n", (long long)globals.sdbg_writer.final_meta().tip_count());
    }

    // --- cleaning ---
    free(globals.lv1_items);
}

} // namespace