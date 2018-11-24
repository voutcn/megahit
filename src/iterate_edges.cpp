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
#include <omp.h>
#include <assert.h>

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <stdexcept>

#include "sequence/kmer_plus.h"
#include "iterate/async_sequence_reader.h"
#include "iterate/kmer_collector.h"
#include "iterate/juncion_index.h"
#include "definitions.h"
#include "options_description.h"
#include "edge_io.h"

using std::string;
using std::vector;

struct IterParam {
  std::string contig_file;
  std::string bubble_file;
  std::string read_file;
  std::string output_prefix;

  unsigned kmer_k;
  unsigned step;
  unsigned num_cpu_threads;
};

struct iter_opt_t {
  string contig_file;
  string bubble_file;
  string read_file;
  int num_cpu_threads;
  int kmer_k;
  int step;
  string output_prefix;

  iter_opt_t() {
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
  desc.AddOption("num_cpu_threads", "t", opt.num_cpu_threads, "number of cpu threads, at least 2. 0 for auto detect.");
  desc.AddOption("kmer_k", "k", opt.kmer_k, "(*) current kmer size.");
  desc.AddOption("step", "s", opt.step,
                 "(*) step for iteration (<= 29). i.e. this iteration is from kmer_k to (kmer_k + step)");
  desc.AddOption("output_prefix", "o", opt.output_prefix,
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

static void InitGlobalData(IterParam &globals) {
  globals.kmer_k = opt.kmer_k;
  globals.step = opt.step;
  globals.num_cpu_threads = opt.num_cpu_threads;
  globals.contig_file = opt.contig_file;
  globals.bubble_file = opt.bubble_file;
  globals.read_file = opt.read_file;
  globals.output_prefix = opt.output_prefix;
}

template<unsigned NumWords, class WordType, class IndexType>
inline bool ReadReadsAndProcessKernel(IterParam &globals, IndexType &index) {
  using KmerType = Kmer<NumWords, WordType>;
  if (KmerType::max_size() < globals.kmer_k + globals.step + 1) {
    return false;
  }
  AsyncReadReader reader(globals.read_file);
  KmerCollector<KmerType> collector(globals.kmer_k + globals.step + 1, globals.output_prefix,
      globals.num_cpu_threads);
  int64_t num_aligned_reads = 0;
  int64_t num_total_reads = 0;

  while (true) {
    auto read_pkg = reader.Next();
    if (read_pkg.size() == 0) {
      break;
    }

#pragma omp parallel for reduction(+: num_aligned_reads)
    for (unsigned i = 0; i < read_pkg.size(); ++i) {
      num_aligned_reads += index.FindNextKmersFromRead(read_pkg, i, &collector) > 0;
    }

    num_total_reads += read_pkg.size();
    xinfo("Processed: %lld, aligned: %lld. Iterative edges: %llu\n",
          num_total_reads, num_aligned_reads, collector.collection().size());
  }

  collector.FlushToFile();
  xinfo("Total: %lld, aligned: %lld. Iterative edges: %llu\n",
        num_total_reads, num_aligned_reads, collector.collection().size());
  return true;
}

template<class IndexType>
static void ReadReadsAndProcess(
    IterParam &globals,
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
    IterParam &globals, const std::string &file_name,
    IndexType *index) {
  AsyncContigReader reader(file_name);
  while (true) {
    auto &pkg = reader.Next();
    auto &contig_pkg = pkg.first;
    auto &mul = pkg.second;

    if (contig_pkg.size() == 0) {
      break;
    }
    index->FeedBatchContigs(contig_pkg, mul);
  }
  xinfo("Number of junction kmers: %lu\n", index->size());
}

template<uint32_t NumWords, typename WordType>
bool IterateToNextK(IterParam &globals) {
  if (Kmer<NumWords, WordType>::max_size() >= globals.kmer_k + 1) {
    JunctionIndex<Kmer<NumWords, WordType>> index(globals.kmer_k, globals.step);
    ReadContigsAndBuildHash(globals, globals.contig_file, &index);
    ReadContigsAndBuildHash(globals, globals.bubble_file, &index);
    ReadReadsAndProcess(globals, index);
    return true;
  }
  return false;
}

int main_iterate(int argc, char *argv[]) {
  AutoMaxRssRecorder recorder;

  IterParam globals;
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