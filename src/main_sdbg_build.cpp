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

/* contact: Dinghua Li <dhli@cs.hku.hk> */

#include <omp.h>
#include <stdio.h>

#include <iostream>
#include <stdexcept>
#include <string>

#include "definitions.h"
#include "sorting/kmer_counter.h"
#include "sorting/read_to_sdbg.h"
#include "sorting/seq_to_sdbg.h"
#include "utils/options_description.h"
#include "utils/utils.h"

int main_kmer_count(int argc, char **argv) {
  AutoMaxRssRecorder recorder;

  // parse option
  OptionsDescription desc;
  KmerCounterOption opt;

  desc.AddOption("kmer_k", "k", opt.k, "kmer size");
  desc.AddOption("min_kmer_frequency", "m", opt.solid_threshold,
                 "min frequency to output an edge");
  desc.AddOption(
      "host_mem", "", opt.host_mem,
      "Max memory to be used. 90% of the free memory is recommended.");
  desc.AddOption("num_cpu_threads", "", opt.n_threads,
                 "number of CPU threads. At least 2.");
  desc.AddOption("read_lib_file", "", opt.read_lib_file,
                 "read library configuration file.");
  desc.AddOption("output_prefix", "", opt.output_prefix, "output prefix");
  desc.AddOption("mem_flag", "", opt.mem_flag,
                 "memory options. 0: minimize memory usage; 1: automatically "
                 "use moderate memory; "
                 "other: use all "
                 "available mem specified by '--host_mem'");

  try {
    desc.Parse(argc, argv);

    if (opt.read_lib_file.empty()) {
      throw std::logic_error("No read library configuration file!");
    }

    if (opt.n_threads == 0) {
      opt.n_threads = omp_get_max_threads();
    }

    if (opt.host_mem == 0) {
      throw std::logic_error("Please specify the host memory!");
    }
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << "Usage: sdbg_builder count --input_file fastx_file -o out"
              << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << desc << std::endl;
    exit(1);
  }

  KmerCounter runner(opt);
  runner.Run();

  return 0;
}

int main_read2sdbg(int argc, char **argv) {
  AutoMaxRssRecorder recorder;

  // parse option the same as kmer_count
  OptionsDescription desc;
  Read2SdbgOption opt;

  desc.AddOption("kmer_k", "k", opt.k, "kmer size");
  desc.AddOption("min_kmer_frequency", "m", opt.solid_threshold,
                 "min frequency to output an edge");
  desc.AddOption(
      "host_mem", "", opt.host_mem,
      "Max memory to be used. 90% of the free memory is recommended.");
  desc.AddOption("num_cpu_threads", "", opt.n_threads,
                 "number of CPU threads. At least 2.");
  desc.AddOption("read_lib_file", "", opt.read_lib_file,
                 "input fast[aq] file, can be gzip'ed. \"-\" for stdin.");
  desc.AddOption("output_prefix", "", opt.output_prefix, "output prefix");
  desc.AddOption("mem_flag", "", opt.mem_flag,
                 "memory options. 0: minimize memory usage; 1: automatically "
                 "use moderate memory; "
                 "other: use all "
                 "available mem specified by '--host_mem'");
  desc.AddOption("need_mercy", "", opt.need_mercy, "to add mercy edges.");

  try {
    desc.Parse(argc, argv);

    if (opt.read_lib_file.empty()) {
      throw std::logic_error("No input file!");
    }

    if (opt.n_threads == 0) {
      opt.n_threads = omp_get_max_threads();
    }

    if (opt.host_mem == 0) {
      throw std::logic_error("Please specify the host memory!");
    }
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cerr
        << "Usage: sdbg_builder read2sdbg --read_lib_file fastx_file -o out"
        << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << desc << std::endl;
    exit(1);
  }

  SeqPkgWithSolidMarker pkg;

  {
    // stage 1
    Read2SdbgS1 runner(opt, &pkg);
    if (opt.solid_threshold > 1) {
      runner.Run();
    } else {
      runner.Initialize();
    }
  }

  {
    // stage 2
    Read2SdbgS2 runner(opt, &pkg);
    runner.Run();
  }

  return 0;
}

int main_seq2sdbg(int argc, char **argv) {
  AutoMaxRssRecorder recorder;

  OptionsDescription desc;
  Seq2SdbgOption opt;

  desc.AddOption("host_mem", "", opt.host_mem,
                 "memory to be used. No more than 95% of the free memory is "
                 "recommended. 0 for auto detect.");
  desc.AddOption("kmer_size", "k", opt.k, "kmer size");
  desc.AddOption("kmer_from", "", opt.k_from, "previous k");
  desc.AddOption("num_cpu_threads", "t", opt.n_threads,
                 "number of CPU threads. At least 2.");
  desc.AddOption("contig", "", opt.contig, "contigs from previous k");
  desc.AddOption("bubble", "", opt.bubble_seq,
                 "bubble sequence from previous k");
  desc.AddOption("addi_contig", "", opt.addi_contig,
                 "additional contigs from previous k");
  desc.AddOption("local_contig", "", opt.local_contig,
                 "local contigs from previous k");
  desc.AddOption(
      "input_prefix", "", opt.input_prefix,
      "files input_prefix.edges.* output by count module, can be gzip'ed.");
  desc.AddOption("output_prefix", "o", opt.output_prefix, "output prefix");
  desc.AddOption("need_mercy", "", opt.need_mercy,
                 "to add mercy edges. The file input_prefix.cand output by "
                 "count module should exist.");
  desc.AddOption("mem_flag", "", opt.mem_flag,
                 "memory options. 0: minimize memory usage; 1: automatically "
                 "use moderate memory; "
                 "other: use all "
                 "available mem specified by '--host_mem'");

  try {
    desc.Parse(argc, argv);

    if (opt.input_prefix.empty() && opt.contig.empty() &&
        opt.addi_contig.empty()) {
      throw std::logic_error("No input files!");
    }

    if (opt.n_threads == 0) {
      opt.n_threads = omp_get_max_threads();
    }

    if (opt.k < 9) {
      throw std::logic_error("kmer size must be >= 9!");
    }

    if (opt.host_mem == 0) {
      throw std::logic_error("Please specify the host memory!");
    }
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << "Usage: sdbg_builder seq2sdbg -k kmer_size --contig "
                 "contigs.fa [--addi_contig "
                 "add.fa] [--input_prefix input] -o out"
              << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << desc << std::endl;
    exit(1);
  }

  SeqToSdbg runner(opt);
  runner.Run();
  return 0;
}
