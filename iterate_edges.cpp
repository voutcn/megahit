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

#include <stdio.h>
#include <pthread.h>
#include <omp.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <zlib.h>

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <stdexcept>

#include "definitions.h"
#include "options_description.h"
#include "atomic_bit_vector.h"
#include "utils.h"
#include "kmer_plus.h"
#include "hash_table.h"
#include "sequence_manager.h"
#include "sequence_package.h"
#include "edge_io.h"

using std::string;
using std::vector;

struct IterateGlobalData {
    char dna_map[256];

    std::string contig_file;
    std::string bubble_file;
    std::string read_file;
    std::string read_format;
    std::string output_prefix;

    int kmer_k;
    int step;
    int next_k1; // = next_k + 1
    int num_cpu_threads;

    // stat
    int64_t num_of_reads;
    int64_t num_of_contigs;
    int64_t num_of_iterative_edges;
    int64_t num_of_remaining_reads;
};


struct iter_opt_t {
    string contig_file;
    string bubble_file;
    string read_file;
    string read_format;
    int num_cpu_threads;
    int kmer_k;
    int step;
    string output_prefix;

    iter_opt_t() {
        read_format = "";
        num_cpu_threads = 0;
        kmer_k = 0;
        step = 0;
    }
};

struct MarginalKmer {
    int64_t s_seq;
    float mul;
} __attribute__((packed));

static iter_opt_t opt;

static void ParseIterOptions(int argc, char *argv[]) {
    OptionsDescription desc;

    desc.AddOption("contig_file", "c", opt.contig_file, "(*) contigs file, fasta/fastq format, output by assembler");
    desc.AddOption("bubble_file", "b", opt.bubble_file, "(*) bubble file, fasta/fastq format, output by assembler");
    desc.AddOption("read_file", "r", opt.read_file, "(*) reads to be aligned. \"-\" for stdin. Can be gzip'ed.");
    desc.AddOption("read_format", "f", opt.read_format, "(*) reads' format. fasta, fastq or binary.");
    desc.AddOption("num_cpu_threads", "t", opt.num_cpu_threads, "number of cpu threads, at least 2. 0 for auto detect.");
    desc.AddOption("kmer_k", "k", opt.kmer_k, "(*) current kmer size.");
    desc.AddOption("step", "s", opt.step, "(*) step for iteration (<= 29). i.e. this iteration is from kmer_k to (kmer_k + step)");
    desc.AddOption("output_prefix", "o", opt.output_prefix, "(*) output_prefix.edges.0 and output_prefix.rr.pb will be created.");

    try {
        desc.Parse(argc, argv);

        if (opt.step + opt.kmer_k >= std::max((int)Kmer<4>::max_size(), (int)GenericKmer::max_size())) {
            std::ostringstream os;
            os << "kmer_k + step must less than " << std::max((int)Kmer<4>::max_size(), (int)GenericKmer::max_size());
            throw std::logic_error(os.str());
        }
        else if (opt.contig_file == "") {
            throw std::logic_error("No contig file!");
        }
        else if (opt.bubble_file == "") {
            throw std::logic_error("No bubble file!");
        }
        else if (opt.read_file == "") {
            throw std::logic_error("No reads file!");
        }
        else if (opt.kmer_k <= 0) {
            throw std::logic_error("Invalid kmer size!");
        }
        else if (opt.step <= 0 || opt.step > 28 || opt.step % 2 == 1) {
            throw std::logic_error("Invalid step size!");
        }
        else if (opt.output_prefix == "") {
            throw std::logic_error("No output prefix!");
        }
        else if (opt.read_format != "binary" && opt.read_format != "fasta" && opt.read_format != "fastq") {
            throw std::logic_error("Invalid read format!");
        }

        if (opt.num_cpu_threads == 0) {
            opt.num_cpu_threads = omp_get_max_threads();
        }

        // must set the number of threads before the parallel hash table declared
        if (opt.num_cpu_threads > 1) {
            omp_set_num_threads(opt.num_cpu_threads - 1);
        }
        else {
            omp_set_num_threads(1);
        }
    }
    catch (std::exception &e) {
        std::cerr << e.what() << std::endl;
        std::cerr << "Usage: " << argv[0] << " [opt]" << std::endl;
        std::cerr << "opt with (*) are must" << std::endl;
        std::cerr << "opt:" << std::endl;
        std::cerr << desc << std::endl;
        exit(1);
    }
}

static void InitGlobalData(IterateGlobalData &globals) {
    globals.kmer_k = opt.kmer_k;
    globals.step = opt.step;
    globals.next_k1 = globals.kmer_k + globals.step + 1;
    globals.num_cpu_threads = opt.num_cpu_threads;
    globals.read_format = opt.read_format;
    globals.contig_file = opt.contig_file;
    globals.bubble_file = opt.bubble_file;
    globals.read_file = opt.read_file;
    globals.output_prefix = opt.output_prefix;
}

static void *ReadContigsThread(void *seq_manager) {
    SequenceManager *sm = (SequenceManager *)seq_manager;

    int64_t kMaxNumContigs = 1 << 22;
    int64_t kMaxNumBases = 1 << 28;
    bool append = false;
    bool reverse = false;
    int discard_flag = contig_flag::kLoop | contig_flag::kIsolated;
    bool extend_loop = false;
    bool calc_depth = false;
    sm->ReadMegahitContigs(kMaxNumContigs, kMaxNumBases, append, reverse, discard_flag, extend_loop, calc_depth);

    return NULL;
}

static void *ReadReadsThread(void *seq_manager) {
    SequenceManager *sm = (SequenceManager *)seq_manager;

    int64_t kMaxNumReads = 1 << 22;
    int64_t kMaxNumBases = 1 << 28;
    bool append = false;
    bool reverse = false;
    sm->ReadShortReads(kMaxNumReads, kMaxNumBases, append, reverse);

    // if (sm->package_->size() > 0) {
    //     for (unsigned i = 0; i < sm->package_->length(0); ++i) {
    //         putchar("ACGT"[sm->package_->get_base(0, i)]);
    //     }
    //     puts("");
    // }

    return NULL;
}

template<uint32_t kNumKmerWord_n, typename kmer_word_n_t>
class WriteEdgeFunc {
  public:
    typedef KmerPlus<kNumKmerWord_n, kmer_word_n_t, uint16_t> kp_t;

    void set_edge_writer(EdgeWriter *edge_writer) {
        edge_writer_ = edge_writer;
    }
    void set_nextk1_size(int nextk1) {
        next_k1_ = nextk1;
        next_k_ = nextk1 - 1;
        last_shift_ = nextk1 % 16;
        last_shift_ = (last_shift_ == 0 ? 0 : 16 - last_shift_) * 2;
        words_per_edge_ = DivCeiling(nextk1 * 2 + kBitsPerMulti_t, 32);
        packed_edges_.resize(omp_get_max_threads() * words_per_edge_);
    }

    void operator() (kp_t &kp) {
        int tid = omp_get_thread_num();
        uint32_t *packed_edge = &packed_edges_[0] + tid * words_per_edge_;
        memset(packed_edge, 0, sizeof(uint32_t) * words_per_edge_);

        int w = 0;
        int end_word = 0;

        for (int j = 0; j < next_k1_; ) {
            w = (w << 2) | kp.kmer.get_base(next_k_ - j);
            ++j;

            if (j % 16 == 0) {
                packed_edge[end_word] = w;
                w = 0;
                end_word++;
            }
        }

        packed_edge[end_word] = (w << last_shift_);
        assert((packed_edge[words_per_edge_ - 1] & kMaxMulti_t) == 0);
        packed_edge[words_per_edge_ - 1] |= kp.ann;
        edge_writer_->write_unsorted(packed_edge, tid);
    }

  private:
    EdgeWriter *edge_writer_;
    int last_shift_;
    int words_per_edge_;
    int next_k1_, next_k_;
    std::vector<uint32_t> packed_edges_;
};

template<uint32_t kNumKmerWord_n, typename kmer_word_n_t, uint32_t kNumKmerWord_p, typename kmer_word_p_t>
static bool ReadReadsAndProcessKernel(IterateGlobalData &globals,
                                      HashTable<KmerPlus<kNumKmerWord_p, kmer_word_p_t, MarginalKmer>, Kmer<kNumKmerWord_p, kmer_word_p_t> > &crusial_kmers) {
    if (Kmer<kNumKmerWord_n, kmer_word_n_t>::max_size() < (unsigned)globals.kmer_k + globals.step + 1) {
        return false;
    }

    HashTable<KmerPlus<kNumKmerWord_n, kmer_word_n_t, uint16_t>, Kmer<kNumKmerWord_n, kmer_word_n_t> > iterative_edges;
    SequencePackage packages[2];
    SequenceManager seq_manager;
    pthread_t input_thread;
    int input_thread_index = 0;
    int64_t num_aligned_reads = 0;
    int64_t num_total_reads = 0;

    if (globals.read_format == "binary") {
        seq_manager.set_file_type(SequenceManager::kBinaryReads);
    }
    else {
        seq_manager.set_file_type(SequenceManager::kFastxReads);
    }

    seq_manager.set_file(globals.read_file);
    seq_manager.set_readlib_type(SequenceManager::kSingle); // PE info not used
    seq_manager.set_package(&packages[input_thread_index]);

    pthread_create(&input_thread, NULL, ReadReadsThread, &seq_manager);
    iterative_edges.reserve(crusial_kmers.size() * 4); // tunable

    while (true) {
        pthread_join(input_thread, NULL);
        SequencePackage &rp = packages[input_thread_index];

        if (rp.size() == 0) {
            break;
        }

        input_thread_index ^= 1;
        seq_manager.set_package(&packages[input_thread_index]);
        pthread_create(&input_thread, NULL, ReadReadsThread, &seq_manager);

        #pragma omp parallel for reduction(+: num_aligned_reads)

        for (unsigned i = 0; i < (unsigned)rp.size(); ++i) {
            int length = rp.length(i);

            if (length < globals.kmer_k + globals.step + 1) {
                continue;
            }

            vector<bool> kmer_exist(length, false);
            vector<float> kmer_mul(length, 0);
            int cur_pos = 0;
            int last_marked_pos = -1;
            Kmer<kNumKmerWord_p, kmer_word_p_t> kmer;

            for (int j = 0; j < globals.kmer_k + 1; ++j) {
                kmer.ShiftAppend(rp.get_base(i, j), globals.kmer_k + 1);
            }

            Kmer<kNumKmerWord_p, kmer_word_p_t> rev_kmer(kmer);
            rev_kmer.ReverseComplement(globals.kmer_k + 1);

            while (cur_pos + globals.kmer_k + 1 <= length) {
                int next_pos = cur_pos + 1;

                if (!kmer_exist[cur_pos]) {
                    auto iter = crusial_kmers.find(kmer);

                    if (iter != crusial_kmers.end()) {
                        kmer_exist[cur_pos] = true;
                        uint64_t s_seq = iter->ann.s_seq;
                        int s_seq_length = s_seq & 63;
                        float mul = iter->ann.mul;
                        kmer_mul[cur_pos] = mul;

                        int j;

                        for (j = 0; j < s_seq_length && cur_pos + globals.kmer_k + 1 + j < length; ++j) {
                            if (rp.get_base(i, cur_pos + globals.kmer_k + 1 + j) == int((s_seq >> (31 - j) * 2) & 3)) {
                                kmer_exist[cur_pos + j + 1] = true;
                                kmer_mul[cur_pos + j + 1] = mul;
                            }
                            else {
                                break;
                            }
                        }

                        last_marked_pos = cur_pos + j;
                        next_pos = last_marked_pos + 1;
                    }
                    else if ((iter = crusial_kmers.find(rev_kmer)) != crusial_kmers.end()) {
                        kmer_exist[cur_pos] = true;
                        uint64_t s_seq = iter->ann.s_seq;
                        int s_seq_length = s_seq & 63;
                        float mul = iter->ann.mul;
                        kmer_mul[cur_pos] = mul;

                        int j;

                        for (j = 0; j < s_seq_length && cur_pos - 1 - j >= 0; ++j) {
                            if (3 - rp.get_base(i, cur_pos - 1 - j) == int((s_seq >> (31 - j) * 2) & 3)) {
                                kmer_exist[cur_pos - 1 - j] = true;
                                kmer_mul[cur_pos - 1 - j] = mul;
                            }
                            else {
                                break;
                            }
                        }
                    }
                }

                if (next_pos + globals.kmer_k + 1 <= length) {
                    while (cur_pos < next_pos) {
                        ++cur_pos;
                        uint8_t c = rp.get_base(i, cur_pos + globals.kmer_k);
                        kmer.ShiftAppend(c, globals.kmer_k + 1);
                        rev_kmer.ShiftPreappend(3 - c, globals.kmer_k + 1);
                    }
                }
                else {
                    break;
                }
            }

            bool aligned = false;
            int acc_exist = 0;
            KmerPlus<kNumKmerWord_n, kmer_word_n_t, uint16_t> kmer_p;
            KmerPlus<kNumKmerWord_n, kmer_word_n_t, uint16_t> rev_kmer_p;

            for (int j = 1; j + globals.kmer_k + 1 <= length; ++j) {
                kmer_mul[j] += kmer_mul[j - 1];
            }

            for (int j = 0, last_j = -globals.kmer_k - 1; j + globals.kmer_k + 1 <= length; ++j) {
                acc_exist = kmer_exist[j] ? acc_exist + 1 : 0;

                if (acc_exist >= globals.step + 1) {
                    if (j - last_j < 8) { // tunable
                        for (int x = last_j + 1; x <= j; ++x) {
                            uint8_t c = rp.get_base(i, x + globals.kmer_k);
                            kmer_p.kmer.ShiftAppend(c, globals.next_k1);
                            rev_kmer_p.kmer.ShiftPreappend(3 - c, globals.next_k1);
                        }
                    }
                    else if (j - last_j < globals.kmer_k + globals.step + 1) {
                        for (int x = last_j + 1; x <= j; ++x) {
                            kmer_p.kmer.ShiftAppend(rp.get_base(i, x + globals.kmer_k), globals.next_k1);
                        }

                        rev_kmer_p.kmer = kmer_p.kmer;
                        rev_kmer_p.kmer.ReverseComplement(globals.next_k1);
                    }
                    else {
                        for (int k = j - globals.step; k < j + globals.kmer_k + 1; ++k) {
                            kmer_p.kmer.ShiftAppend(rp.get_base(i, k), globals.next_k1);
                        }

                        rev_kmer_p.kmer = kmer_p.kmer;
                        rev_kmer_p.kmer.ReverseComplement(globals.next_k1);
                    }

                    float mul = (kmer_mul[j] - (j - (globals.step + 1) >= 0 ? kmer_mul[j - (globals.step + 1)] : 0)) / (globals.step + 1);
                    assert(mul <= kMaxMulti_t + 1);

                    if (kmer_p.kmer < rev_kmer_p.kmer) {
                        kmer_p.ann = std::min(kMaxMulti_t, int(mul + 0.5));
                        iterative_edges.find_or_insert(kmer_p);
                    }
                    else {
                        rev_kmer_p.ann = std::min(kMaxMulti_t, int(mul + 0.5));
                        iterative_edges.find_or_insert(rev_kmer_p);
                    }

                    last_j = j;
                    aligned = true;
                }
            }

            if (aligned) {
                ++num_aligned_reads;
            }
        }

        num_total_reads += rp.size();

        if (num_total_reads % (16 << 22) == 0) {
            xlog("Processed: %lld, aligned: %lld. Iterative edges: %llu\n", (long long)num_total_reads, (long long)num_aligned_reads, (unsigned long long)iterative_edges.size());
        }
    }

    xlog("Total: %lld, aligned: %lld. Iterative edges: %llu\n", (long long)num_total_reads, (long long)num_aligned_reads, (unsigned long long)iterative_edges.size());

    // write iterative edges
    if (iterative_edges.size() > 0) {
        xlog("Writing iterative edges...\n");
        EdgeWriter edge_writer;

        uint32_t next_k = globals.kmer_k + globals.step;
        omp_set_num_threads(globals.num_cpu_threads);

        edge_writer.set_num_threads(globals.num_cpu_threads);
        edge_writer.set_file_prefix(globals.output_prefix);
        edge_writer.set_unsorted();
        edge_writer.set_kmer_size(next_k);
        edge_writer.init_files();

        WriteEdgeFunc<kNumKmerWord_n, kmer_word_n_t> write_edge_func;
        write_edge_func.set_nextk1_size(next_k + 1);
        write_edge_func.set_edge_writer(&edge_writer);

        iterative_edges.for_each(write_edge_func);
    }

    return true;
}

template<uint32_t kNumKmerWord_p, typename kmer_word_p_t>
static void ReadReadsAndProcess(IterateGlobalData &globals,
                                HashTable<KmerPlus<kNumKmerWord_p, kmer_word_p_t, MarginalKmer>, Kmer<kNumKmerWord_p, kmer_word_p_t> > &crusial_kmers) {
    if (ReadReadsAndProcessKernel<1, uint64_t>(globals, crusial_kmers)) return;

    if (ReadReadsAndProcessKernel<3, uint32_t>(globals, crusial_kmers)) return;

    if (ReadReadsAndProcessKernel<2, uint64_t>(globals, crusial_kmers)) return;

    if (ReadReadsAndProcessKernel<5, uint32_t>(globals, crusial_kmers)) return;

    if (ReadReadsAndProcessKernel<3, uint64_t>(globals, crusial_kmers)) return;

    if (ReadReadsAndProcessKernel<7, uint32_t>(globals, crusial_kmers)) return;

    if (ReadReadsAndProcessKernel<4, uint64_t>(globals, crusial_kmers)) return;

    if (ReadReadsAndProcessKernel<kUint32PerKmerMaxK, uint32_t>(globals, crusial_kmers)) return;

    assert (false);
}

template<uint32_t kNumKmerWord_p, typename kmer_word_p_t>
static void ReadContigsAndBuildHash(IterateGlobalData &globals, std::string file_name, 
                                    HashTable<KmerPlus<kNumKmerWord_p, kmer_word_p_t, MarginalKmer>, Kmer<kNumKmerWord_p, kmer_word_p_t> > &crusial_kmers) {
    SequencePackage packages[2];
    std::vector<float> f_muls[2];
    SequenceManager seq_manager;
    int input_thread_index = 0;
    pthread_t input_thread;

    seq_manager.set_file_type(SequenceManager::kMegahitContigs);
    seq_manager.set_package(&packages[input_thread_index]);
    seq_manager.set_float_multiplicity_vector(&f_muls[input_thread_index]);
    seq_manager.set_file(file_name);

    pthread_create(&input_thread, NULL, ReadContigsThread, &seq_manager);

    while (true) {
        pthread_join(input_thread, NULL);
        SequencePackage &cp = packages[input_thread_index];
        std::vector<float> &f_mul = f_muls[input_thread_index];

        if (cp.size() == 0) {
            break;
        }

        input_thread_index ^= 1;
        seq_manager.set_float_multiplicity_vector(&f_muls[input_thread_index]);
        seq_manager.set_package(&packages[input_thread_index]);
        pthread_create(&input_thread, NULL, ReadContigsThread, &seq_manager);

        #pragma omp parallel for

        for (unsigned i = 0; i < cp.size(); ++i) {
            if ((int)cp.length(i) < globals.kmer_k + 1) {
                continue;
            }

            KmerPlus<kNumKmerWord_p, kmer_word_p_t, MarginalKmer> kmer_p;
            Kmer<kNumKmerWord_p, kmer_word_p_t> &kmer = kmer_p.kmer;

            for (int j = 0; j < globals.kmer_k + 1; ++j) {
                kmer.ShiftAppend(cp.get_base(i, j), globals.kmer_k + 1);
            }

            uint64_t s_seq = 0;
            int s_length = std::min(globals.step, (int)cp.length(i) - (globals.kmer_k + 1));

            for (int j = 0; j < globals.step && j < s_length; ++j) {
                s_seq |= uint64_t(cp.get_base(i, j + globals.kmer_k + 1)) << (31 - j) * 2;
            }

            s_seq |= s_length;
            kmer_p.ann.s_seq = s_seq;
            kmer_p.ann.mul = f_mul[i];
            crusial_kmers.find_or_insert(kmer_p);

            if ((int)cp.length(i) > globals.kmer_k + 1) {
                for (int j = 0; j < globals.kmer_k + 1; ++j) {
                    kmer.ShiftAppend(3 - cp.get_base(i, cp.length(i) - 1 - j), globals.kmer_k + 1);
                }

                s_seq = 0;

                for (int j = 0; j < globals.step && j < s_length; ++j) {
                    s_seq |= uint64_t(3 - cp.get_base(i, cp.length(i) - 1 - (globals.kmer_k + 1) - j)) << (31 - j) * 2;
                }

                s_seq |= s_length;
                kmer_p.ann.s_seq = s_seq;
                kmer_p.ann.mul = f_mul[i];
                crusial_kmers.find_or_insert(kmer_p);
            }
        }
    }

    xlog("Number of crusial kmers: %lu\n", crusial_kmers.size());
}

template <uint32_t kNumKmerWord_p, typename kmer_word_p_t>
bool IterateToNextK(IterateGlobalData &globals) {
    if (Kmer<kNumKmerWord_p, kmer_word_p_t>::max_size() >= (unsigned)globals.kmer_k + 1) {
        HashTable<KmerPlus<kNumKmerWord_p, kmer_word_p_t, MarginalKmer>, Kmer<kNumKmerWord_p, kmer_word_p_t> > crusial_kmers;
        ReadContigsAndBuildHash(globals, globals.bubble_file, crusial_kmers);
        ReadContigsAndBuildHash(globals, globals.contig_file, crusial_kmers);
        ReadReadsAndProcess(globals, crusial_kmers);
        return true;
    }

    return false;
}

int main_iterate(int argc, char *argv[]) {
    AutoMaxRssRecorder recorder;

    IterateGlobalData globals;
    ParseIterOptions(argc, argv);
    InitGlobalData(globals);

    while (true) {
        if (IterateToNextK<1, uint64_t>(globals)) break;

        if (IterateToNextK<3, uint32_t>(globals)) break;

        if (IterateToNextK<2, uint64_t>(globals)) break;

        if (IterateToNextK<5, uint32_t>(globals)) break;

        if (IterateToNextK<3, uint64_t>(globals)) break;

        if (IterateToNextK<7, uint32_t>(globals)) break;

        if (IterateToNextK<4, uint64_t>(globals)) break;

        if (IterateToNextK<kUint32PerKmerMaxK, uint32_t>(globals)) break;

        assert(false);
    }

    return 0;
}