/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong & L3 Bioinformatics
 * Limited
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

#include <assert.h>
#include <omp.h>
#include <stdio.h>

#include <algorithm>
#include <functional>
#include <iostream>
#include <list>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "definitions.h"
#include "iterate/contig_flank_index.h"
#include "iterate/kmer_collector.h"
#include "sequence/io/async_sequence_reader.h"
#include "utils/options_description.h"

using std::string;
using std::vector;

namespace {

struct LocalAsmOption {
  string contig_file;
  string bubble_file;
  string read_file;
  int num_cpu_threads{0};
  int kmer_k{0};
  int step{0};
  string output_prefix;
} opt;

static void ParseIterOptions(int argc, char *argv[]) {
  OptionsDescription desc;

  desc.AddOption("contig_file", "c", opt.contig_file,
                 "(*) contigs file, fasta/fastq format, output by assembler");
  desc.AddOption("bubble_file", "b", opt.bubble_file,
                 "(*) bubble file, fasta/fastq format, output by assembler");
  desc.AddOption("read_file", "r", opt.read_file,
                 "(*) reads to be aligned. \"-\" for stdin. Can be gzip'ed.");
  desc.AddOption("num_cpu_threads", "t", opt.num_cpu_threads,
                 "number of cpu threads, at least 2. 0 for auto detect.");
  desc.AddOption("kmer_k", "k", opt.kmer_k, "(*) current kmer size.");
  desc.AddOption("step", "s", opt.step,
                 "(*) step for iteration (<= 28). i.e. this iteration is from "
                 "kmer_k to (kmer_k + step)");
  desc.AddOption("output_prefix", "o", opt.output_prefix,
                 "(*) output_prefix.edges.0 will be created.");

  try {
    desc.Parse(argc, argv);

    if (opt.step + opt.kmer_k >=
        static_cast<int>(
            std::max(Kmer<4>::max_size(), GenericKmer::max_size()))) {
      std::ostringstream os;
      os << "kmer_k + step must less than "
         << std::max(Kmer<4>::max_size(), GenericKmer::max_size());
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
    if (opt.num_cpu_threads > 1) {
      omp_set_num_threads(opt.num_cpu_threads - 1);
    } else {
      omp_set_num_threads(1);
    }
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << "Usage: " << argv[0] << " [opt]" << std::endl;
    std::cerr << "opt with (*) are must" << std::endl;
    std::cerr << "opt:" << std::endl;
    std::cerr << desc << std::endl;
    exit(1);
  }
}

}  // namespace

template <class KmerType, class IndexType>
static bool ReadReadsAndProcessKernel(const LocalAsmOption &opt,
                                      const IndexType &index) {
  if (KmerType::max_size() < static_cast<unsigned>(opt.kmer_k + opt.step + 1)) {
    return false;
  }
  xinfo("Selected kmer type size for next k: {}\n", sizeof(KmerType));
  BinaryReader binary_reader(opt.read_file);
  AsyncSequenceReader reader(&binary_reader);
  KmerCollector<KmerType> collector(opt.kmer_k + opt.step + 1,
                                    opt.output_prefix);
  int64_t num_aligned_reads = 0;
  int64_t num_total_reads = 0;

  while (true) {
    const auto &read_pkg = reader.Next();
    if (read_pkg.seq_count() == 0) {
      break;
    }
    num_aligned_reads += index.FindNextKmersFromReads(read_pkg, &collector);
    num_total_reads += read_pkg.seq_count();
    xinfo("Processed: {}, aligned: {}. Iterative edges: {}\n", num_total_reads,
          num_aligned_reads, collector.collection().size());
  }
  collector.FlushToFile();
  xinfo("Total: {}, aligned: {}. Iterative edges: {}\n", num_total_reads,
        num_aligned_reads, collector.collection().size());
  return true;
}

template <class IndexType>
static void ReadReadsAndProcess(const LocalAsmOption &opt,
                                const IndexType &index) {
  if (ReadReadsAndProcessKernel<Kmer<1, uint64_t>>(opt, index)) return;
  if (ReadReadsAndProcessKernel<Kmer<3, uint32_t>>(opt, index)) return;
  if (ReadReadsAndProcessKernel<Kmer<2, uint64_t>>(opt, index)) return;
  if (ReadReadsAndProcessKernel<Kmer<5, uint32_t>>(opt, index)) return;
  if (ReadReadsAndProcessKernel<Kmer<3, uint64_t>>(opt, index)) return;
  if (ReadReadsAndProcessKernel<Kmer<7, uint32_t>>(opt, index)) return;
  if (ReadReadsAndProcessKernel<Kmer<4, uint64_t>>(opt, index)) return;
  if (ReadReadsAndProcessKernel<Kmer<kUint32PerKmerMaxK, uint32_t>>(opt, index))
    return;
  xfatal("k is too large!\n");
}

template <class IndexType>
static void ReadContigsAndBuildIndex(const LocalAsmOption &opt,
                                     const std::string &file_name,
                                     IndexType *index) {
  AsyncContigReader reader(file_name);
  while (true) {
    auto &pkg = reader.Next();
    auto &contig_pkg = pkg.first;
    auto &mul = pkg.second;
    if (contig_pkg.seq_count() == 0) {
      break;
    }
    xinfo("Read {} contigs\n", contig_pkg.seq_count());
    index->FeedBatchContigs(contig_pkg, mul);
    xinfo("Number of flank kmers: {}\n", index->size());
  }
}

struct BaseRunner {
  virtual ~BaseRunner() = default;
  virtual void Run(const LocalAsmOption &opt) = 0;
  virtual uint32_t max_k() const = 0;
};

template <class KmerType>
struct Runner : public BaseRunner {
  ~Runner() override = default;
  void Run(const LocalAsmOption &opt) override {
    xinfo("Selected kmer type size for k: {}\n", sizeof(KmerType));
    ContigFlankIndex<KmerType> index(opt.kmer_k, opt.step);
    ReadContigsAndBuildIndex(opt, opt.contig_file, &index);
    ReadContigsAndBuildIndex(opt, opt.bubble_file, &index);
    ReadReadsAndProcess(opt, index);
  }
  uint32_t max_k() const override { return KmerType::max_size(); }
};

int main_iterate(int argc, char *argv[]) {
  AutoMaxRssRecorder recorder;
  ParseIterOptions(argc, argv);

  std::list<std::unique_ptr<BaseRunner>> runners;
  runners.emplace_back(new Runner<Kmer<1, uint64_t>>());
  runners.emplace_back(new Runner<Kmer<3, uint32_t>>());
  runners.emplace_back(new Runner<Kmer<2, uint64_t>>());
  runners.emplace_back(new Runner<Kmer<5, uint32_t>>());
  runners.emplace_back(new Runner<Kmer<3, uint64_t>>());
  runners.emplace_back(new Runner<Kmer<7, uint32_t>>());
  runners.emplace_back(new Runner<Kmer<4, uint64_t>>());
  runners.emplace_back(new Runner<Kmer<kUint32PerKmerMaxK, uint32_t>>());

  for (auto &runner : runners) {
    if (runner->max_k() >= static_cast<uint32_t>(opt.kmer_k + 1)) {
      runner->Run(opt);
      return 0;
    }
  }

  xfatal("k is too large!\n");
  return 1;
}