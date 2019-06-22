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

#include "cx1/cx1_kmer_count.h"
#include "cx1/cx1_read2sdbg.h"
#include "cx1/cx1_seq2sdbg.h"
#include "definitions.h"
#include "utils/options_description.h"
#include "utils/utils.h"

int main_kmer_count(int argc, char **argv) {
  AutoMaxRssRecorder recorder;

  // parse option
  OptionsDescription desc;
  count_opt_t opt;

  desc.AddOption("kmer_k", "k", opt.kmer_k, "kmer size");
  desc.AddOption("min_kmer_frequency", "m", opt.kmer_freq_threshold, "min frequency to output an edge");
  desc.AddOption("host_mem", "", opt.host_mem, "Max memory to be used. 90% of the free memory is recommended.");
  desc.AddOption("num_cpu_threads", "", opt.num_cpu_threads, "number of CPU threads. At least 2.");
  desc.AddOption("read_lib_file", "", opt.read_lib_file, "read library configuration file.");
  desc.AddOption("assist_seq", "", opt.assist_seq_file,
                 "input assisting fast[aq] file (FILE_NAME.info should exist), can be gzip'ed.");
  desc.AddOption("output_prefix", "", opt.output_prefix, "output prefix");
  desc.AddOption("mem_flag", "", opt.mem_flag,
                 "memory options. 0: minimize memory usage; 1: automatically use moderate memory; "
                 "other: use all "
                 "available mem specified by '--host_mem'");

  try {
    desc.Parse(argc, argv);

    if (opt.read_lib_file == "") {
      throw std::logic_error("No read library configuration file!");
    }

    if (opt.num_cpu_threads == 0) {
      opt.num_cpu_threads = omp_get_max_threads();
    }

    if (opt.host_mem == 0) {
      throw std::logic_error("Please specify the host memory!");
    }
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << "Usage: sdbg_builder count --input_file fastx_file -o out" << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << desc << std::endl;
    exit(1);
  }

  CX1KmerCount runner(opt.kmer_k, opt.kmer_freq_threshold, opt.host_mem, opt.mem_flag, opt.num_cpu_threads,
      opt.read_lib_file, opt.assist_seq_file, opt.output_prefix);
  xinfo("Host memory to be used: %lld\n", (long long)opt.host_mem);
  xinfo("Number CPU threads: %d\n", opt.num_cpu_threads);
  runner.run();

  return 0;
}

int main_read2sdbg(int argc, char **argv) {
  AutoMaxRssRecorder recorder;

  // parse option the same as kmer_count
  OptionsDescription desc;
  read2sdbg_opt_t opt;

  desc.AddOption("kmer_k", "k", opt.kmer_k, "kmer size");
  desc.AddOption("min_kmer_frequency", "m", opt.kmer_freq_threshold, "min frequency to output an edge");
  desc.AddOption("host_mem", "", opt.host_mem, "Max memory to be used. 90% of the free memory is recommended.");
  desc.AddOption("num_cpu_threads", "", opt.num_cpu_threads, "number of CPU threads. At least 2.");
  desc.AddOption("read_lib_file", "", opt.read_lib_file, "input fast[aq] file, can be gzip'ed. \"-\" for stdin.");
  desc.AddOption("assist_seq", "", opt.assist_seq_file,
                 "input assisting fast[aq] file (FILE_NAME.info should exist), can be gzip'ed.");
  desc.AddOption("output_prefix", "", opt.output_prefix, "output prefix");
  desc.AddOption("mem_flag", "", opt.mem_flag,
                 "memory options. 0: minimize memory usage; 1: automatically use moderate memory; "
                 "other: use all "
                 "available mem specified by '--host_mem'");
  desc.AddOption("need_mercy", "", opt.need_mercy, "to add mercy edges.");

  try {
    desc.Parse(argc, argv);

    if (opt.read_lib_file == "") {
      throw std::logic_error("No input file!");
    }

    if (opt.num_cpu_threads == 0) {
      opt.num_cpu_threads = omp_get_max_threads();
    }

    if (opt.host_mem == 0) {
      throw std::logic_error("Please specify the host memory!");
    }
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << "Usage: sdbg_builder read2sdbg --read_lib_file fastx_file -o out" << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << desc << std::endl;
    exit(1);
  }

  cx1_read2sdbg::read2sdbg_global_t globals;
  globals.kmer_k = opt.kmer_k;
  globals.kmer_freq_threshold = opt.kmer_freq_threshold;
  globals.host_mem = opt.host_mem;
  globals.num_cpu_threads = opt.num_cpu_threads;
  globals.read_lib_file = opt.read_lib_file;
  globals.assist_seq_file = opt.assist_seq_file;
  globals.output_prefix = opt.output_prefix;
  globals.mem_flag = opt.mem_flag;
  globals.need_mercy = opt.need_mercy;
  globals.cx1.reset(new cx1_read2sdbg::CX1Read2SdbgS1());
  globals.cx1->g_ = &globals;

  // stage1
  if (opt.kmer_freq_threshold > 1) {
    globals.cx1->run();
  } else {
    globals.cx1->prepare_func_(globals);
  }

  // stage2
  globals.cx1.reset(new cx1_read2sdbg::CX1Read2SdbgS2());
  globals.cx1->g_ = &globals;
  globals.cx1->run();

  return 0;
}

int main_seq2sdbg(int argc, char **argv) {
  AutoMaxRssRecorder recorder;

  OptionsDescription desc;
  seq2sdbg_opt_t opt;

  desc.AddOption("host_mem", "", opt.host_mem,
                 "memory to be used. No more than 95% of the free memory is recommended. 0 for auto detect.");
  desc.AddOption("kmer_size", "k", opt.kmer_k, "kmer size");
  desc.AddOption("kmer_from", "", opt.kmer_from, "previous k");
  desc.AddOption("num_cpu_threads", "t", opt.num_cpu_threads, "number of CPU threads. At least 2.");
  desc.AddOption("contig", "", opt.contig, "contigs from previous k");
  desc.AddOption("bubble", "", opt.bubble_seq, "bubble sequence from previous k");
  desc.AddOption("addi_contig", "", opt.addi_contig, "additional contigs from previous k");
  desc.AddOption("local_contig", "", opt.local_contig, "local contigs from previous k");
  desc.AddOption("input_prefix", "", opt.input_prefix,
                 "files input_prefix.edges.* output by count module, can be gzip'ed.");
  desc.AddOption("output_prefix", "o", opt.output_prefix, "output prefix");
  desc.AddOption("need_mercy", "", opt.need_mercy,
                 "to add mercy edges. The file input_prefix.cand output by count module should exist.");
  desc.AddOption("mem_flag", "", opt.mem_flag,
                 "memory options. 0: minimize memory usage; 1: automatically use moderate memory; "
                 "other: use all "
                 "available mem specified by '--host_mem'");

  try {
    desc.Parse(argc, argv);

    if (opt.input_prefix == "" && opt.contig == "" && opt.addi_contig == "") {
      throw std::logic_error("No input files!");
    }

    if (opt.num_cpu_threads == 0) {
      opt.num_cpu_threads = omp_get_max_threads();
    }

    if (opt.kmer_k < 9) {
      throw std::logic_error("kmer size must be >= 9!");
    }

    if (opt.host_mem == 0) {
      throw std::logic_error("Please specify the host memory!");
    }
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << "Usage: sdbg_builder seq2sdbg -k kmer_size --contig contigs.fa [--addi_contig "
                 "add.fa] [--input_prefix "
                 "input] -o out"
              << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << desc << std::endl;
    exit(1);
  }

  cx1_seq2sdbg::seq2sdbg_global_t globals;
  globals.host_mem = opt.host_mem;
  globals.num_cpu_threads = opt.num_cpu_threads;
  globals.input_prefix = opt.input_prefix;
  globals.output_prefix = opt.output_prefix;
  globals.contig = opt.contig;
  globals.bubble_seq = opt.bubble_seq;
  globals.addi_contig = opt.addi_contig;
  globals.local_contig = opt.local_contig;
  globals.mem_flag = opt.mem_flag;
  globals.kmer_k = opt.kmer_k;
  globals.kmer_from = opt.kmer_from;
  globals.need_mercy = opt.need_mercy;

  xinfo("Host memory to be used: %lld\n", (long long)globals.host_mem);
  xinfo("Number CPU threads: %d\n", globals.num_cpu_threads);

  // set & run cx1
  globals.cx1.reset(new cx1_seq2sdbg::CX1Seq2Sdbg());
  globals.cx1->g_ = &globals;
  globals.cx1->run();
  return 0;
}
