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

#include "sequence/kmer_plus.h"
#include "iterate/juncion_index.h"
#include "definitions.h"
#include "options_description.h"
#include "kmlib/kmbitvector.h"
#include "utils.h"
#include "sparsepp/sparsepp/spp.h"
#include "sequence_manager.h"
#include "sequence_package.h"
#include "edge_io.h"

using std::string;
using std::vector;

struct IterateGlobalData {
  std::string contig_file;
  std::string bubble_file;
  std::string read_file;
  std::string read_format;
  std::string output_prefix;

  unsigned kmer_k;
  unsigned step;
  unsigned next_k1; // = next_k + 1
  unsigned num_cpu_threads;
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

static iter_opt_t opt;

static void ParseIterOptions(int argc, char *argv[]) {
  OptionsDescription desc;

  desc.AddOption("contig_file", "c", opt.contig_file, "(*) contigs file, fasta/fastq format, output by assembler");
  desc.AddOption("bubble_file", "b", opt.bubble_file, "(*) bubble file, fasta/fastq format, output by assembler");
  desc.AddOption("read_file", "r", opt.read_file, "(*) reads to be aligned. \"-\" for stdin. Can be gzip'ed.");
  desc.AddOption("read_format", "f", opt.read_format, "(*) reads' format. fasta, fastq or binary.");
  desc.AddOption("num_cpu_threads", "t", opt.num_cpu_threads, "number of cpu threads, at least 2. 0 for auto detect.");
  desc.AddOption("kmer_k", "k", opt.kmer_k, "(*) current kmer size.");
  desc.AddOption("step",
                 "s",
                 opt.step,
                 "(*) step for iteration (<= 29). i.e. this iteration is from kmer_k to (kmer_k + step)");
  desc.AddOption("output_prefix",
                 "o",
                 opt.output_prefix,
                 "(*) output_prefix.edges.0 and output_prefix.rr.pb will be created.");

  try {
    desc.Parse(argc, argv);

    if (opt.step + opt.kmer_k >= std::max((int) Kmer<4>::max_size(), (int) GenericKmer::max_size())) {
      std::ostringstream os;
      os << "kmer_k + step must less than " << std::max((int) Kmer<4>::max_size(), (int) GenericKmer::max_size());
      throw std::logic_error(os.str());
    } else if (opt.contig_file == "") {
      throw std::logic_error("No contig file!");
    } else if (opt.bubble_file == "") {
      throw std::logic_error("No bubble file!");
    } else if (opt.read_file == "") {
      throw std::logic_error("No reads file!");
    } else if (opt.kmer_k <= 0) {
      throw std::logic_error("Invalid kmer size!");
    } else if (opt.step <= 0 || opt.step > 28 || opt.step % 2 == 1) {
      throw std::logic_error("Invalid step size!");
    } else if (opt.output_prefix == "") {
      throw std::logic_error("No output prefix!");
    } else if (opt.read_format != "binary" && opt.read_format != "fasta" && opt.read_format != "fastq") {
      throw std::logic_error("Invalid read format!");
    }

    if (opt.num_cpu_threads == 0) {
      opt.num_cpu_threads = omp_get_max_threads();
    }

    // must set the number of threads before the parallel hash table declared
    if (opt.num_cpu_threads > 1) {
      omp_set_num_threads(opt.num_cpu_threads - 1);
    } else {
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
  auto *sm = (SequenceManager *) seq_manager;

  int64_t kMaxNumContigs = 1 << 22;
  int64_t kMaxNumBases = 1 << 28;
  bool append = false;
  bool reverse = false;
  int discard_flag = contig_flag::kLoop | contig_flag::kStandalone;
  bool extend_loop = false;
  bool calc_depth = false;
  sm->ReadMegahitContigs(kMaxNumContigs, kMaxNumBases, append, reverse, discard_flag, extend_loop, calc_depth);

  return nullptr;
}

static void *ReadReadsThread(void *seq_manager) {
  auto *sm = (SequenceManager *) seq_manager;

  int64_t kMaxNumReads = 1 << 22;
  int64_t kMaxNumBases = 1 << 28;
  bool append = false;
  bool reverse = false;
  sm->ReadShortReads(kMaxNumReads, kMaxNumBases, append, reverse);

  return nullptr;
}

template <class KmerType>
class WriteEdgeFunc {
 public:
  using kp_t = KmerPlus<KmerType::kNumWords, typename KmerType::word_type, mul_t>;

  void set_edge_writer(EdgeWriter *edge_writer) {
    edge_writer_ = edge_writer;
  }
  void set_nextk1_size(unsigned nextk1) {
    next_k1_ = nextk1;
    next_k_ = nextk1 - 1;
    last_shift_ = nextk1 % 16;
    last_shift_ = (last_shift_ == 0 ? 0 : 16 - last_shift_) * 2;
    words_per_edge_ = DivCeiling(nextk1 * 2 + kBitsPerMul, 32);
    packed_edges_.resize(omp_get_max_threads() * words_per_edge_);
  }

  void operator()(const kp_t &kp) {
    int tid = omp_get_thread_num();
    uint32_t *packed_edge = &packed_edges_[0] + tid * words_per_edge_;
    memset(packed_edge, 0, sizeof(uint32_t) * words_per_edge_);

    int w = 0;
    int end_word = 0;

    for (unsigned j = 0; j < next_k1_;) {
      w = (w << 2) | kp.kmer.GetBase(next_k_ - j);
      ++j;

      if (j % 16 == 0) {
        packed_edge[end_word] = w;
        w = 0;
        end_word++;
      }
    }

    packed_edge[end_word] = (w << last_shift_);
    assert((packed_edge[words_per_edge_ - 1] & kMaxMul) == 0);
    packed_edge[words_per_edge_ - 1] |= kp.aux;
    edge_writer_->write_unsorted(packed_edge, tid);
  }

 private:
  EdgeWriter *edge_writer_;
  unsigned last_shift_;
  unsigned words_per_edge_;
  unsigned next_k1_, next_k_;
  std::vector<uint32_t> packed_edges_;
};

template <unsigned NumWords, class WordType, class IndexType>
inline bool ReadReadsAndProcessKernel(
    IterateGlobalData &globals,
    IndexType &index) {
  using KmerType = Kmer<NumWords, WordType>;
  if (KmerType::max_size() < (unsigned) globals.kmer_k + globals.step + 1) {
    return false;
  }

  spp::sparse_hash_set<KmerPlus<NumWords, WordType, uint16_t>, KmerHash> iterative_edges;
  SequencePackage packages[2];
  SequenceManager seq_manager;
  pthread_t input_thread;
  int input_thread_index = 0;
  int64_t num_aligned_reads = 0;
  int64_t num_total_reads = 0;

  std::mutex lock;

  if (globals.read_format == "binary") {
    seq_manager.set_file_type(SequenceManager::kBinaryReads);
  } else {
    seq_manager.set_file_type(SequenceManager::kFastxReads);
  }

  seq_manager.set_file(globals.read_file);
  seq_manager.set_readlib_type(SequenceManager::kSingle); // PE info not used
  seq_manager.set_package(&packages[input_thread_index]);

  pthread_create(&input_thread, nullptr, ReadReadsThread, &seq_manager);
  iterative_edges.reserve(index.size() * 4); // tunable

  while (true) {
    pthread_join(input_thread, nullptr);
    SequencePackage &rp = packages[input_thread_index];

    if (rp.size() == 0) {
      break;
    }

    input_thread_index ^= 1;
    seq_manager.set_package(&packages[input_thread_index]);
    pthread_create(&input_thread, nullptr, ReadReadsThread, &seq_manager);
    vector<KmerPlus<NumWords, WordType, mul_t>> iter_kmers;

#pragma omp parallel for reduction(+: num_aligned_reads) private(iter_kmers)
    for (unsigned i = 0; i < rp.size(); ++i) {
      iter_kmers.clear();
      num_aligned_reads += index.FindNextKmersFromRead(rp, i, &iter_kmers) > 0;
      std::lock_guard<std::mutex> lk(lock);
      for (auto &item : iter_kmers) {
        iterative_edges.emplace(item);
      }
    }

    num_total_reads += rp.size();
    xinfo("Processed: %lld, aligned: %lld. Iterative edges: %llu\n",
          (long long) num_total_reads,
          (long long) num_aligned_reads,
          (unsigned long long) iterative_edges.size());
  }

  xinfo("Total: %lld, aligned: %lld. Iterative edges: %llu\n",
        (long long) num_total_reads,
        (long long) num_aligned_reads,
        (unsigned long long) iterative_edges.size());

  // write iterative edges
  if (iterative_edges.size() > 0) {
    xinfo("Writing iterative edges...\n");
    EdgeWriter edge_writer;

    uint32_t next_k = globals.kmer_k + globals.step;
    omp_set_num_threads(globals.num_cpu_threads);

    edge_writer.set_num_threads(globals.num_cpu_threads);
    edge_writer.set_file_prefix(globals.output_prefix);
    edge_writer.set_unsorted();
    edge_writer.set_kmer_size(next_k);
    edge_writer.init_files();

    WriteEdgeFunc<KmerType> write_edge_func;
    write_edge_func.set_nextk1_size(next_k + 1);
    write_edge_func.set_edge_writer(&edge_writer);
#pragma omp parallel for
    for (size_t bucket_i = 0; bucket_i < iterative_edges.bucket_count(); ++bucket_i) {
      for (auto it = iterative_edges.begin(bucket_i), end = iterative_edges.end(bucket_i); it != end; ++it) {
        write_edge_func(*it);
      }
    }
  }

  return true;
}

template <class IndexType>
static void ReadReadsAndProcess(
    IterateGlobalData &globals,
    IndexType &index) {
  if (ReadReadsAndProcessKernel<1, uint64_t>(globals, index)) return;
  if (ReadReadsAndProcessKernel<3, uint32_t>(globals, index)) return;
  if (ReadReadsAndProcessKernel<2, uint64_t>(globals, index)) return;
  if (ReadReadsAndProcessKernel<5, uint32_t>(globals, index)) return;
  if (ReadReadsAndProcessKernel<3, uint64_t>(globals, index)) return;
  if (ReadReadsAndProcessKernel<7, uint32_t>(globals, index)) return;
  if (ReadReadsAndProcessKernel<4, uint64_t>(globals, index)) return;
  if (ReadReadsAndProcessKernel<kUint32PerKmerMaxK, uint32_t>(globals, index)) return;
  assert (false);
}

template<class IndexType>
static void ReadContigsAndBuildHash(
    IterateGlobalData &globals, std::string file_name,
    IndexType *index) {
  SequencePackage packages[2];
  std::vector<float> f_muls[2];
  SequenceManager seq_manager;
  int input_thread_index = 0;
  pthread_t input_thread;

  seq_manager.set_file_type(SequenceManager::kMegahitContigs);
  seq_manager.set_package(&packages[input_thread_index]);
  seq_manager.set_float_multiplicity_vector(&f_muls[input_thread_index]);
  seq_manager.set_file(file_name);

  pthread_create(&input_thread, nullptr, ReadContigsThread, &seq_manager);

  while (true) {
    pthread_join(input_thread, nullptr);
    SequencePackage &cp = packages[input_thread_index];
    std::vector<float> &f_mul = f_muls[input_thread_index];

    if (cp.size() == 0) {
      break;
    }

    input_thread_index ^= 1;
    seq_manager.set_float_multiplicity_vector(&f_muls[input_thread_index]);
    seq_manager.set_package(&packages[input_thread_index]);
    pthread_create(&input_thread, nullptr, ReadContigsThread, &seq_manager);
    
    index->FeedBatchContigs(cp, f_mul);
  }
  xinfo("Number of junction kmers: %lu\n", index->size());
}

template<uint32_t NumWords, typename WordType>
bool IterateToNextK(IterateGlobalData &globals) {
  if (Kmer<NumWords, WordType>::max_size() >= globals.kmer_k + 1) {
    JunctionIndex<Kmer<NumWords, WordType>> index(globals.kmer_k, globals.step);
    ReadContigsAndBuildHash(globals, globals.bubble_file, &index);
    ReadContigsAndBuildHash(globals, globals.contig_file, &index);
    ReadReadsAndProcess(globals, index);
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