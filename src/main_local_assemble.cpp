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

#include <omp.h>
#include <iostream>
#include <string>

#include "localasm/local_assembler.h"
#include "utils/options_description.h"
#include "utils/utils.h"

namespace {
struct Option {
  std::string contig_file;
  std::string lib_file_prefix;

  int kmin{11};
  int kmax{41};
  int step{6};
  int seed_kmer{31};

  int min_contig_len{200};
  int sparsity{8};
  double similarity{0.8};
  int min_mapping_len{75};

  int num_threads{0};
  std::string output_file;
} opt;

void ParseLocalAsmOptions(int argc, char *argv[]) {
  OptionsDescription desc;

  desc.AddOption("contig_file", "c", opt.contig_file, "contig file");
  desc.AddOption("lib_file_prefix", "l", opt.lib_file_prefix, "lib file prefix");
  desc.AddOption("kmin", "", opt.kmin, "");
  desc.AddOption("kmax", "", opt.kmax, "");
  desc.AddOption("step", "", opt.step, "");
  desc.AddOption("seed_kmer", "", opt.seed_kmer, "kmer size for seeding alignments");
  desc.AddOption("min_contig_len", "", opt.min_contig_len, "");
  desc.AddOption("min_mapping_len", "", opt.min_mapping_len, "");
  desc.AddOption("sparsity", "", opt.sparsity, "sparsity of hash mapper");
  desc.AddOption("similarity", "", opt.similarity, "alignment similarity threshold");
  desc.AddOption("num_threads", "t", opt.num_threads, "");
  desc.AddOption("output_file", "o", opt.output_file, "");

  try {
    desc.Parse(argc, argv);
    if (opt.contig_file == "") {
      throw std::logic_error("no contig file!");
    }
    if (opt.lib_file_prefix == "") {
      throw std::logic_error("no read file!");
    }
    if (opt.output_file == "") {
      throw std::logic_error("no output file!");
    }
    if (opt.num_threads == 0) {
      opt.num_threads = omp_get_max_threads();
    }
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << "Usage: " << argv[0] << " -c contigs.fa -r reads.fq -o out.local_contig.fa" << std::endl;
    std::cerr << "options:" << std::endl;
    std::cerr << desc << std::endl;
    exit(1);
  }
}
}  // namespace

int main_local(int argc, char **argv) {
  AutoMaxRssRecorder recorder;

  ParseLocalAsmOptions(argc, argv);

  omp_set_num_threads(opt.num_threads);

  LocalAssembler la(opt.min_contig_len, opt.seed_kmer, opt.sparsity);
  la.set_kmer(opt.kmin, opt.kmax, opt.step);
  la.set_mapping_threshold(opt.similarity, opt.min_mapping_len);
  la.set_local_file(opt.output_file);

  la.ReadContigs(opt.contig_file);
  la.BuildHashMapper();
  la.AddReadLib(opt.lib_file_prefix);
  la.EstimateInsertSize();
  la.MapToContigs();
  la.LocalAssemble();

  return 0;
}