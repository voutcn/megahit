/*
 *  MEGAHIT
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

#include "iterate_edges.h"
#include <stdio.h>
#include <pthread.h>
#include <omp.h>
#include <stdlib.h>
#include <assert.h>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include "definitions.h"
#include "io-utility.h"
#include "options_description.h"
#include "atomic_bit_vector.h"
#include "timer.h"
#include "kmer_plus.h"
#include "hash_table.h"

using std::string;
using std::vector;

struct Options {
    string contigs_file;
    string contigs_multi_file;
    string addi_contig_file;
    string addi_multi_file;
    string read_file;
    string read_format;
    int num_cpu_threads;
    int kmer_k;
    int step;
    int max_read_len;
    string output_prefix;

    Options() {
        read_format = "";
        num_cpu_threads = 0;
        kmer_k = 0;
        step = 0;
        max_read_len = 0;
    }

    string output_edges_file() {
        return output_prefix + "edges.0";
    }

    string output_read_file() {
        return output_prefix + "rr.pb";
    }
} options;

static void ParseOptions(int argc, char *argv[]) {
    OptionsDescription desc;

    desc.AddOption("contigs_file", "c", options.contigs_file, "(*) contigs file, fasta/fastq format, output by assembler");
    desc.AddOption("multi_file", "m", options.contigs_multi_file, "(*) contigs's multiplicity file output by assembler");
    desc.AddOption("addi_contig_file", "", options.addi_contig_file, "additional contigs file, fasta/fastq format, output by assembler if remove low local");
    desc.AddOption("addi_multi_file", "", options.addi_multi_file, "contigs's multiplicity file, output by assembler if remove low local");
    desc.AddOption("read_file", "r", options.read_file, "(*) reads to be aligned. \"-\" for stdin. Can be gzip'ed.");
    desc.AddOption("read_format", "f", options.read_format, "(*) reads' format. fasta, fastq or binary.");
    desc.AddOption("num_cpu_threads", "t", options.num_cpu_threads, "number of cpu threads, at least 2. 0 for auto detect.");
    desc.AddOption("kmer_k", "k", options.kmer_k, "(*) current kmer size.");
    desc.AddOption("step", "s", options.step, "(*) step for iteration (<= 29). i.e. this iteration is from kmer_k to (kmer_k + step)");
    desc.AddOption("output_prefix", "o", options.output_prefix, "(*) output_prefix.edges.0 and output_prefix.rr.pb will be created.");
    desc.AddOption("max_read_len", "l", options.max_read_len, "(*) max read length of all reads.");

    try {
        desc.Parse(argc, argv);
        if (options.step + options.kmer_k >= (int)Kmer<4>::max_size()) {
            std::ostringstream os;
            os << "kmer_k + step must less than " << Kmer<4>::max_size();
            throw std::logic_error(os.str());
        } else if (options.contigs_file == "") {
            throw std::logic_error("No contig file!");
        } else if (options.contigs_multi_file == "") {
            throw std::logic_error("No contig's multiplicity file!");
        } else if (options.read_file == "") {
            throw std::logic_error("No reads file!");
        } else if (options.kmer_k <= 0) {
            throw std::logic_error("Invalid kmer size!");
        } else if (options.step <= 0) {
            throw std::logic_error("Invalid step size!");
        } else if (options.output_prefix == "") {
            throw std::logic_error("No output prefix!");
        } else if (options.read_format != "binary" && options.read_format != "fasta" && options.read_format != "fastq") {
            throw std::logic_error("Invalid read format!");
        } else if (options.max_read_len == 0) {
            throw std::logic_error("Invalid max read length!");
        }

        if (options.num_cpu_threads == 0) {
            options.num_cpu_threads = omp_get_max_threads();
        }
        // must set the number of threads before the parallel hash table declared
        if (options.num_cpu_threads > 1) {
            omp_set_num_threads(options.num_cpu_threads - 1);
        } else {
            omp_set_num_threads(1);
        }
    } catch (std::exception &e) {
        std::cerr << e.what() << std::endl;
        std::cerr << "Usage: " << argv[0] << " [options]" << std::endl;
        std::cerr << "options with (*) are must" << std::endl;
        std::cerr << "options:" << std::endl;
        std::cerr << desc << std::endl;
        exit(1);
    }
}

static void InitGlobalData(IterateGlobalData &globals) {
    for (int i = 0; i < 256; ++i) {
        globals.dna_map[i] = 2;
    }
    globals.dna_map['A'] = 0;
    globals.dna_map['C'] = 1;
    globals.dna_map['G'] = 2;
    globals.dna_map['T'] = 3;

    globals.kmer_k = options.kmer_k;
    globals.step = options.step;
    globals.max_read_len = options.max_read_len;
    globals.num_cpu_threads = options.num_cpu_threads;

    if (string(options.read_format) == "fastq") {
        globals.read_format = IterateGlobalData::kFastq;
    } else if (string(options.read_format) == "fasta") {
        globals.read_format = IterateGlobalData::kFasta;
    } else if (string(options.read_format) == "binary") {
        globals.read_format = IterateGlobalData::kBinary;
    } else {
        fprintf(stderr, "Cannot identify read format!\n");
        exit(1);
    }

    globals.contigs_file = gzopen(options.contigs_file.c_str(), "r");
    globals.contigs_multi_file = gzopen(options.contigs_multi_file.c_str(), "r");

    if (string(options.read_file) == "-") {
        globals.read_file = gzdopen(fileno(stdin), "r");
    } else {
        globals.read_file = gzopen(options.read_file.c_str(), "r");
    }
    assert(globals.contigs_file != NULL);
    assert(globals.contigs_multi_file != NULL);
    assert(globals.read_file != NULL);

    if (options.addi_contig_file != "") {
        globals.addi_contig_file = gzopen(options.addi_contig_file.c_str(), "r");
        globals.addi_multi_file = gzopen(options.addi_multi_file.c_str(), "r");
        assert(globals.addi_multi_file != NULL);
        assert(globals.addi_contig_file != NULL);
    } else {
        globals.addi_contig_file = NULL;
        globals.addi_multi_file = NULL;
    }

    globals.output_edge_file = OpenFileAndCheck((string(options.output_prefix) + ".edges.0").c_str(), "wb");
    globals.output_read_file = OpenFileAndCheck((string(options.output_prefix) + ".rr.pb").c_str(), "wb"); // remaining reads packed binary
    assert(globals.output_edge_file != NULL);
    assert(globals.output_read_file != NULL);
}

static void ClearGlobalData(IterateGlobalData &globals) {
    gzclose(globals.contigs_file);
    gzclose(globals.read_file);
    if (globals.addi_contig_file != NULL) {
        gzclose(globals.addi_contig_file);
        gzclose(globals.addi_multi_file);
    }
    fclose(globals.output_edge_file);
    fclose(globals.output_read_file);
}

struct ReadContigsThreadData {
    ContigPackage *contig_package;
    kseq_t *seq;
    IterateGlobalData *globals;
    gzFile *multi_file;
};

static void* ReadContigsThread(void* data) {
    ContigPackage &package = *(((ReadContigsThreadData*)data)->contig_package);
    kseq_t *seq = ((ReadContigsThreadData*)data)->seq;
    IterateGlobalData &globals = *(((ReadContigsThreadData*)data)->globals);
    gzFile &multi_file = *(((ReadContigsThreadData*)data)->multi_file);
    char *dna_map = globals.dna_map;

    printf("Reading contigs...\n");
    package.ReadContigs(seq, dna_map);
    package.ReadMultiplicity(multi_file);
    printf("Read %lu contigs, total length: %lu\n", package.size(), package.seqs.length());
    return NULL;
}

struct ReadReadsThreadData {
    ReadPackage *read_package;
    kseq_t *seq;
    IterateGlobalData *globals;
};

static void* ReadReadsThread(void* data) {
    ReadPackage &package = *(((ReadReadsThreadData*)data)->read_package);
    IterateGlobalData &globals = *(((ReadReadsThreadData*)data)->globals);
    kseq_t *seq = ((ReadReadsThreadData*)data)->seq;
    package.clear();

    if (globals.read_format == IterateGlobalData::kFastq || globals.read_format == IterateGlobalData::kFasta) {
        package.ReadFastxReads(seq, globals.dna_map);
    } else {
        package.ReadBinaryReads(globals.read_file);
    }
    return NULL;
}

template<uint32_t kNumKmerWord_n, typename kmer_word_n_t, uint32_t kNumKmerWord_p, typename kmer_word_p_t>
static bool ReadReadsAndProcessKernel(IterateGlobalData &globals, 
                                      HashTable<KmerPlus<kNumKmerWord_p, kmer_word_p_t, uint64_t>, Kmer<kNumKmerWord_p, kmer_word_p_t> > &crusial_kmers)
{
    if (Kmer<kNumKmerWord_n, kmer_word_n_t>::max_size() < (unsigned)globals.kmer_k + globals.step + 1) {
        return false;
    }
    else 
    {
        HashTable<KmerPlus<kNumKmerWord_n, kmer_word_n_t, uint16_t>, Kmer<kNumKmerWord_n, kmer_word_n_t> > iterative_edges;
        ReadPackage packages[2];
        packages[0].init(globals.max_read_len);
        packages[1].init(globals.max_read_len);
        kseq_t *seq = NULL;
        if (globals.read_format != IterateGlobalData::kBinary) {
            seq = kseq_init(globals.read_file);
        }
        int input_thread_index = 0;
        static const int kWordsPerEdge = ((globals.kmer_k + globals.step + 1) * 2 + kBitsPerMulti_t + 31) / 32;
        uint32_t packed_edge[kWordsPerEdge];

        int64_t num_aligned_reads = 0;
        int64_t num_total_reads = 0;

        pthread_t input_thread;
        ReadReadsThreadData input_thread_data;
        input_thread_data.read_package = &packages[input_thread_index];
        input_thread_data.seq = seq;
        input_thread_data.globals = &globals;

        pthread_create(&input_thread, NULL, ReadReadsThread, &input_thread_data);
        iterative_edges.reserve(crusial_kmers.size() * 4); // tunable
        AtomicBitVector is_aligned;

        if (globals.num_cpu_threads == 1) {
            pthread_join(input_thread, NULL);
        }

        while (true) {
            if (globals.num_cpu_threads > 1) {
                pthread_join(input_thread, NULL);
            }

            if (packages[input_thread_index].num_of_reads == 0) {
                break;
            }

            input_thread_index ^= 1;
            input_thread_data.read_package = &packages[input_thread_index];
            pthread_create(&input_thread, NULL, ReadReadsThread, &input_thread_data);
            ReadPackage &cur_package = packages[input_thread_index ^ 1];
            is_aligned.reset(cur_package.num_of_reads);

            if (globals.num_cpu_threads == 1) {
                pthread_join(input_thread, NULL);
            }

     #pragma omp parallel for
            for (unsigned i = 0; i < (unsigned)cur_package.num_of_reads; ++i) {
                int length = cur_package.length(i);
                assert(length <= cur_package.max_read_len);
                if (length < globals.kmer_k + globals.step + 1) {
                    continue;
                }

                vector<bool> kmer_exist(length, false);
                int cur_pos = 0;
                int last_marked_pos = -1;
                Kmer<kNumKmerWord_p, kmer_word_p_t> kmer(globals.kmer_k);
                for (int j = 0; j < globals.kmer_k; ++j) {
                    kmer.ShiftAppend(cur_package.CharAt(i, j));
                }
                
                Kmer<kNumKmerWord_p, kmer_word_p_t> rev_kmer(kmer);
                rev_kmer.ReverseComplement();

                while (cur_pos + globals.kmer_k <= length) {
                    int next_pos = cur_pos + 1;
                    if (!kmer_exist[cur_pos]) {
                        auto iter = crusial_kmers.find(kmer);
                        if (iter != crusial_kmers.end()) {
                            kmer_exist[cur_pos] = true;
                            int64_t s_seq = iter->ann;
                            int s_seq_length = s_seq & 63;
                            int j;
                            for (j = 0; j < s_seq_length && cur_pos + globals.kmer_k + j < length; ++j) {
                                if (cur_package.CharAt(i, cur_pos + globals.kmer_k + j) == ((s_seq >> (31 - j) * 2) & 3)) {
                                    kmer_exist[cur_pos + j + 1] = true;
                                } else {
                                    break;
                                }
                            }

                            last_marked_pos = cur_pos + j;
                            next_pos = last_marked_pos + 1;
                        } else if ((iter = crusial_kmers.find(rev_kmer)) != crusial_kmers.end()) {
                            kmer_exist[cur_pos] = true;
                            int64_t s_seq = iter->ann;
                            int s_seq_length = s_seq & 63;
                            int j;
                            for (j = 0; j < s_seq_length && cur_pos - 1 - j > last_marked_pos; ++j) {
                                if (3 - cur_package.CharAt(i, cur_pos - 1 - j) == ((s_seq >> (31 - j) * 2) & 3)) {
                                    kmer_exist[cur_pos - 1 - j] = true;
                                } else {
                                    break;
                                }
                            }
                        }
                    }

                    if (next_pos + globals.kmer_k <= length) {
                        while (cur_pos < next_pos) {
                            ++cur_pos;
                            uint8_t c = cur_package.CharAt(i, cur_pos + globals.kmer_k - 1);
                            kmer.ShiftAppend(c);
                            rev_kmer.ShiftPreappend(3 - c);
                        }
                    } else {
                        break;
                    }
                }

                bool aligned = false;
                KmerPlus<kNumKmerWord_n, kmer_word_n_t, uint16_t> kmer_p(globals.kmer_k + 1 + globals.step);
                KmerPlus<kNumKmerWord_n, kmer_word_n_t, uint16_t> rev_kmer_p(globals.kmer_k + 1 + globals.step);

                for (int j = 0, last_j = -globals.kmer_k, acc_exist = 0; j + globals.kmer_k <= length; ++j) {
                    acc_exist = kmer_exist[j] ? acc_exist + 1 : 0;

                    if (acc_exist >= globals.step + 2) {
                        if (j - last_j < 8) { // tunable
                            for (int x = last_j + 1; x <= j; ++x) {
                                uint8_t c = cur_package.CharAt(i, x + globals.kmer_k - 1);
                                kmer_p.kmer.ShiftAppend(c);
                                rev_kmer_p.kmer.ShiftPreappend(3 - c);
                            }
                        } else if (j - last_j < globals.kmer_k + globals.step + 1) {
                            for (int x = last_j + 1; x <= j; ++x) {
                                kmer_p.kmer.ShiftAppend(cur_package.CharAt(i, x + globals.kmer_k - 1));
                            }
                            rev_kmer_p.kmer = kmer_p.kmer;
                            rev_kmer_p.kmer.ReverseComplement();
                        } else {
                            for (int k = j - globals.step - 1; k < j + globals.kmer_k; ++k) {
                                kmer_p.kmer.ShiftAppend(cur_package.CharAt(i, k));
                            }
                            rev_kmer_p.kmer = kmer_p.kmer;
                            rev_kmer_p.kmer.ReverseComplement();
                        }

                        if (kmer_p.kmer < rev_kmer_p.kmer) {
                            KmerPlus<kNumKmerWord_n, kmer_word_n_t, uint16_t> &kp = iterative_edges.find_or_insert_with_lock(kmer_p);
                            if (kp.ann < kMaxMulti_t) { ++kp.ann; }
                            iterative_edges.unlock(kmer_p);
                        } else {
                            KmerPlus<kNumKmerWord_n, kmer_word_n_t, uint16_t> &kp = iterative_edges.find_or_insert_with_lock(rev_kmer_p);
                            if (kp.ann < kMaxMulti_t) { ++kp.ann; }
                            iterative_edges.unlock(rev_kmer_p);
                        }
                        last_j = j;
                        aligned = true;
                    }
                }
                if (aligned) {
                    is_aligned.set(i);
     #pragma omp atomic
                    ++num_aligned_reads;
                }
            }

            num_total_reads += cur_package.num_of_reads;

            for (int64_t i = 0; i < cur_package.num_of_reads; ++i) {
                if (is_aligned.get(i)) {
                    fwrite(cur_package.packed_reads + i * cur_package.words_per_read,
                           sizeof(uint32_t), cur_package.words_per_read, globals.output_read_file);
                }
            }

            if (num_total_reads % (16 * cur_package.kMaxNumReads) == 0) {
                printf("Total: %lld, aligned: %lld. Iterative edges: %llu\n", (long long)num_total_reads, (long long)num_aligned_reads, (unsigned long long)iterative_edges.size());
            }
        }
        printf("Total: %lld, aligned: %lld. Iterative edges: %llu\n", (long long)num_total_reads, (long long)num_aligned_reads, (unsigned long long)iterative_edges.size());

        kseq_destroy(seq);

        printf("Writing iterative edges...\n");
        int next_k = globals.step + globals.kmer_k;
        int last_shift = (next_k + 1) % 16;
        last_shift = (last_shift == 0 ? 0 : 16 - last_shift) * 2;
        for (auto iter = iterative_edges.begin(); iter != iterative_edges.end(); ++iter) {
            memset(packed_edge, 0, sizeof(uint32_t) * kWordsPerEdge);
            int w = 0;
            int end_word = 0;
            for (int j = 0; j < next_k + 1; ) {
                w = (w << 2) | iter->kmer.get_base(next_k - j);
                ++j;
                if (j % 16 == 0) {
                    packed_edge[end_word] = w;
                    w = 0;
                    end_word++;
                }
            }
            packed_edge[end_word] = (w << last_shift);
            assert((packed_edge[kWordsPerEdge - 1] & kMaxMulti_t) == 0);
            packed_edge[kWordsPerEdge - 1] |= iter->ann;
            fwrite(packed_edge, sizeof(uint32_t), kWordsPerEdge, globals.output_edge_file);
        }
    }
    return true;
}

template<uint32_t kNumKmerWord_p, typename kmer_word_p_t>
static void ReadReadsAndProcess(IterateGlobalData &globals, 
                                HashTable<KmerPlus<kNumKmerWord_p, kmer_word_p_t, uint64_t>, Kmer<kNumKmerWord_p, kmer_word_p_t> > &crusial_kmers)
{
    if (ReadReadsAndProcessKernel<1, uint64_t>(globals, crusial_kmers)) return;
    if (ReadReadsAndProcessKernel<3, uint32_t>(globals, crusial_kmers)) return;
    if (ReadReadsAndProcessKernel<2, uint64_t>(globals, crusial_kmers)) return;
    if (ReadReadsAndProcessKernel<5, uint32_t>(globals, crusial_kmers)) return;
    if (ReadReadsAndProcessKernel<3, uint64_t>(globals, crusial_kmers)) return;
    if (ReadReadsAndProcessKernel<7, uint32_t>(globals, crusial_kmers)) return;
    if (ReadReadsAndProcessKernel<4, uint64_t>(globals, crusial_kmers)) return;
    assert (false);
}

template<uint32_t kNumKmerWord_p, typename kmer_word_p_t>
static void ReadContigsAndBuildHash(IterateGlobalData &globals, 
                                    HashTable<KmerPlus<kNumKmerWord_p, kmer_word_p_t, uint64_t>, Kmer<kNumKmerWord_p, kmer_word_p_t> > &crusial_kmers, 
                                    bool is_addi_contigs)
{
    ContigPackage packages[2];
    kseq_t *seq;
    if (is_addi_contigs) {
        seq = kseq_init(globals.addi_contig_file);
    } else {
        seq = kseq_init(globals.contigs_file);
    }
    int input_thread_index = 0;
    bool is_first_round = !is_addi_contigs;
    static const int kWordsPerEdge = ((globals.kmer_k + globals.step + 1) * kBitsPerEdgeChar + kBitsPerMulti_t + 31) / 32;
    uint32_t packed_edge[kWordsPerEdge];

    pthread_t input_thread;
    ReadContigsThreadData input_thread_data;
    input_thread_data.contig_package = &packages[input_thread_index];
    input_thread_data.seq = seq;
    input_thread_data.globals = &globals;
    if (!is_addi_contigs) {
        input_thread_data.multi_file = &globals.contigs_multi_file;
    } else {
        input_thread_data.multi_file = &globals.addi_multi_file;
    }

    pthread_create(&input_thread, NULL, ReadContigsThread, &input_thread_data);
    if (globals.num_cpu_threads == 1) {
        pthread_join(input_thread, NULL);
    }

    while (true) {
        if (globals.num_cpu_threads > 1) {
            pthread_join(input_thread, NULL);
        }

        if (packages[input_thread_index].size() == 0) {
            break;
        }

        input_thread_index ^= 1;
        input_thread_data.contig_package = &packages[input_thread_index];
        pthread_create(&input_thread, NULL, ReadContigsThread, &input_thread_data);
        ContigPackage &cur_package = packages[input_thread_index ^ 1];

        if (globals.num_cpu_threads == 1) {
            pthread_join(input_thread, NULL);
        }

        if (!is_addi_contigs) {
 #pragma omp parallel for
            for (unsigned i = 0; i < cur_package.size(); ++i) {
                if (cur_package.seq_lengths[i] < globals.kmer_k) {
                    continue;
                }

                KmerPlus<kNumKmerWord_p, kmer_word_p_t, uint64_t> kmer_p(globals.kmer_k);
                Kmer<kNumKmerWord_p, kmer_word_p_t> &kmer = kmer_p.kmer;
                for (int j = 0; j < globals.kmer_k; ++j) {
                    kmer.ShiftAppend(cur_package.CharAt(i, j));
                }
                uint64_t s_seq = 0;
                int s_length = std::min(globals.step, cur_package.seq_lengths[i] - globals.kmer_k);
                for (int j = 0; j < globals.step && j < s_length; ++j) {
                    s_seq |= uint64_t(cur_package.CharAt(i, j + globals.kmer_k)) << (31 - j) * 2;
                }
                s_seq |= s_length;
                KmerPlus<kNumKmerWord_p, kmer_word_p_t, uint64_t> &kp = crusial_kmers.find_or_insert(kmer_p);
                kp.ann = s_seq;

                if (cur_package.seq_lengths[i] > globals.kmer_k) {
                    for (int j = 0; j < globals.kmer_k; ++j) {
                        kmer.ShiftAppend(3 - cur_package.CharAt(i, cur_package.seq_lengths[i] - 1 - j));
                    }

                    s_seq = 0;
                    for (int j = 0; j < globals.step && j < s_length; ++j) {
                        s_seq |= uint64_t(3 - cur_package.CharAt(i, cur_package.seq_lengths[i] - globals.kmer_k - 1 - j)) << (31 - j) * 2;
                    }
                    s_seq |= s_length;
                    KmerPlus<kNumKmerWord_p, kmer_word_p_t, uint64_t> &kp = crusial_kmers.find_or_insert(kmer_p);
                    kp.ann = s_seq;
                }
            }
        }

        if (is_first_round) {
            uint32_t next_k = globals.kmer_k + globals.step;
            fwrite(&next_k, sizeof(uint32_t), 1, globals.output_edge_file);
            fwrite(&kWordsPerEdge, sizeof(uint32_t), 1, globals.output_edge_file);
            is_first_round = false;
        }

        int next_k = globals.kmer_k + globals.step;
        int last_shift = (next_k + 1) % 16;
        last_shift = (last_shift == 0 ? 0 : 16 - last_shift) * 2;
        for (unsigned i = 0; i < cur_package.size(); ++i) {
            if (cur_package.seq_lengths[i] < next_k + 1) {
                continue;
            }

            double multiplicity_prev = cur_package.multiplicity[i];
            uint16_t multiplicity;
            // convert the multiplicity from k to k+s+1
            {
                int num_kmer = cur_package.seq_lengths[i] - globals.kmer_k + 1;
                int num_nextk1 = cur_package.seq_lengths[i] - (next_k + 1) + 1;
                int internal_max = std::min(next_k + 1 - globals.kmer_k + 1, num_nextk1);
                int num_external = internal_max - 1;
                int num_internal = num_kmer - num_external * 2;

                double exp_num_kmer = (double)num_external * (num_external + 1) / (next_k + 1 - globals.kmer_k + 1)
                                      + (double)internal_max / (next_k + 1 - globals.kmer_k + 1) * num_internal;
                exp_num_kmer *= multiplicity_prev;
                multiplicity = std::min(int(exp_num_kmer * globals.kmer_k / (next_k + 1) / num_nextk1 + 0.5), kMaxMulti_t);
            }

            memset(packed_edge, 0, sizeof(uint32_t) * kWordsPerEdge);

            int w = 0;
            int end_word = 0;
            for (int j = 0; j < next_k + 1; ) {
                w = (w << 2) | cur_package.CharAt(i, next_k - j);
                ++j;
                if (j % 16 == 0) {
                    packed_edge[end_word] = w;
                    w = 0;
                    end_word++;
                }
            }
            packed_edge[end_word] = (w << last_shift);
            packed_edge[kWordsPerEdge - 1] |= multiplicity;

            fwrite(packed_edge, sizeof(uint32_t), kWordsPerEdge, globals.output_edge_file);

            for (int j = next_k + 1; j < cur_package.seq_lengths[i]; ++j) {
                packed_edge[kWordsPerEdge - 1] ^= multiplicity;
                packed_edge[next_k / 16] &= ~(3 << (15 - next_k % 16) * 2);
                for (int k = kWordsPerEdge - 1; k > 0; --k) {
                    packed_edge[k] >>= 2;
                    packed_edge[k] |= packed_edge[k - 1] << 30;
                }
                packed_edge[0] >>= 2;
                packed_edge[0] |= cur_package.CharAt(i, j) << 30;
                assert((packed_edge[kWordsPerEdge - 1] & kMaxMulti_t) == 0);
                packed_edge[kWordsPerEdge - 1] |= multiplicity;
                fwrite(packed_edge, sizeof(uint32_t), kWordsPerEdge, globals.output_edge_file);
            }
        }
    }
    printf("Number of crusial kmers: %lu\n", crusial_kmers.size());

    kseq_destroy(seq);
}

template <uint32_t kNumKmerWord_p, typename kmer_word_p_t>
bool IterateToNextK(IterateGlobalData &globals)
{
    if (Kmer<kNumKmerWord_p, kmer_word_p_t>::max_size() >= (unsigned)globals.kmer_k)
    {
        HashTable<KmerPlus<kNumKmerWord_p, kmer_word_p_t, uint64_t>, Kmer<kNumKmerWord_p, kmer_word_p_t> > crusial_kmers;
        ReadContigsAndBuildHash(globals, crusial_kmers, false);

        if (options.addi_multi_file != "") {
            ReadContigsAndBuildHash(globals, crusial_kmers, true);
        }

        ReadReadsAndProcess(globals, crusial_kmers);
        return true;
    }
    return false;
}

static AutoMaxRssRecorder recorder;

int main(int argc, char *argv[]) {
    // set stdout line buffered
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    IterateGlobalData globals;
    ParseOptions(argc, argv);
    InitGlobalData(globals);

    while (true) {
        if (IterateToNextK<1, uint64_t>(globals)) break;
        if (IterateToNextK<3, uint32_t>(globals)) break;
        if (IterateToNextK<2, uint64_t>(globals)) break;
        if (IterateToNextK<5, uint32_t>(globals)) break;
        if (IterateToNextK<3, uint64_t>(globals)) break;
        if (IterateToNextK<7, uint32_t>(globals)) break;
        if (IterateToNextK<4, uint64_t>(globals)) break;
        assert(false);
    }

    ClearGlobalData(globals);
    return 0;
}