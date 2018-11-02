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

#include <stdio.h>
#include <omp.h>

#include <iostream>
#include <string>
#include <stdexcept>

#include "cx1_kmer_count.h"
#include "cx1_read2sdbg.h"
#include "cx1_seq2sdbg.h"
#include "lv2_gpu_functions.h"
#include "options_description.h"
#include "utils.h"
#include "definitions.h"

int main_kmer_count(int argc, char **argv) {
    AutoMaxRssRecorder recorder;

    // parse option
    OptionsDescription desc;
    count_opt_t opt;

    desc.AddOption("kmer_k", "k", opt.kmer_k, "kmer size");
    desc.AddOption("min_kmer_frequency", "m", opt.kmer_freq_threshold, "min frequency to output an edge");
    desc.AddOption("host_mem", "", opt.host_mem, "Max memory to be used. 90% of the free memory is recommended.");
    desc.AddOption("gpu_mem", "", opt.gpu_mem, "gpu memory to be used. 0 for auto detect.");
    desc.AddOption("num_cpu_threads", "", opt.num_cpu_threads, "number of CPU threads. At least 2.");
    desc.AddOption("num_output_threads", "", opt.num_output_threads, "number of threads for output. Must be less than num_cpu_threads");
    desc.AddOption("read_lib_file", "", opt.read_lib_file, "read library configuration file.");
    desc.AddOption("assist_seq", "", opt.assist_seq_file, "input assisting fast[aq] file (FILE_NAME.info should exist), can be gzip'ed.");
    desc.AddOption("output_prefix", "", opt.output_prefix, "output prefix");
    desc.AddOption("mem_flag", "", opt.mem_flag, "memory options. 0: minimize memory usage; 1: automatically use moderate memory; other: use all available mem specified by '--host_mem'");

    try {
        desc.Parse(argc, argv);

        if (opt.read_lib_file == "") {
            throw std::logic_error("No read library configuration file!");
        }

        if (opt.num_cpu_threads == 0) {
            opt.num_cpu_threads = omp_get_max_threads();
        }

        if (opt.num_output_threads == 0) {
            opt.num_output_threads = std::max(1, opt.num_cpu_threads / 3);
        }

        if (opt.host_mem == 0) {
            throw std::logic_error("Please specify the host memory!");
        }

        if (opt.gpu_mem == 0) {
#ifdef USE_GPU
            size_t free_gpu_mem, total_gpu_mem;
            get_cuda_memory(free_gpu_mem, total_gpu_mem);
            opt.gpu_mem = free_gpu_mem;
#else
            opt.gpu_mem = 0;
#endif
        }

#ifdef USE_GPU

        if (opt.num_output_threads >= opt.num_cpu_threads) {
            throw std::logic_error("Number of output threads must be less than number of CPU threads!");
        }

#endif

    }
    catch (std::exception &e) {
        std::cerr << e.what() << std::endl;
        std::cerr << "Usage: sdbg_builder count --input_file fastx_file -o out" << std::endl;
        std::cerr << "Options:" << std::endl;
        std::cerr << desc << std::endl;
        exit(1);
    }

    cx1_kmer_count::count_global_t globals;
    globals.kmer_k = opt.kmer_k;
    globals.kmer_freq_threshold = opt.kmer_freq_threshold;
    globals.host_mem = opt.host_mem;
    globals.gpu_mem = opt.gpu_mem;
    globals.num_cpu_threads = opt.num_cpu_threads;
    globals.num_output_threads = opt.num_output_threads;
    globals.read_lib_file = opt.read_lib_file.c_str();
    globals.assist_seq_file = opt.assist_seq_file;
    globals.output_prefix = opt.output_prefix.c_str();
    globals.mem_flag = opt.mem_flag;

    xlog("Host memory to be used: %lld\n", (long long)globals.host_mem);
    xlog("Number CPU threads: %d\n", globals.num_cpu_threads);

    // set & run cx1
    globals.cx1.g_ = &globals;
    globals.cx1.encode_lv1_diff_base_func_ = cx1_kmer_count::encode_lv1_diff_base;
    globals.cx1.prepare_func_ = cx1_kmer_count::read_input_prepare;
    globals.cx1.lv0_calc_bucket_size_func_ = cx1_kmer_count::lv0_calc_bucket_size;
    globals.cx1.init_global_and_set_cx1_func_ = cx1_kmer_count::init_global_and_set_cx1;
    globals.cx1.lv1_fill_offset_func_ = cx1_kmer_count::lv1_fill_offset;
    globals.cx1.lv1_sort_and_proc = cx1_kmer_count::lv1_direct_sort_and_count;
    globals.cx1.lv2_extract_substr_func_ = cx1_kmer_count::lv2_extract_substr;
    globals.cx1.lv2_sort_func_ = cx1_kmer_count::lv2_sort;
    globals.cx1.lv2_pre_output_partition_func_ = cx1_kmer_count::lv2_pre_output_partition;
    globals.cx1.lv2_output_func_ = cx1_kmer_count::lv2_output;
    globals.cx1.lv2_post_output_func_ = cx1_kmer_count::lv2_post_output;
    globals.cx1.post_proc_func_ = cx1_kmer_count::post_proc;

    globals.cx1.run();

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
    desc.AddOption("gpu_mem", "", opt.gpu_mem, "gpu memory to be used. 0 for auto detect.");
    desc.AddOption("num_cpu_threads", "", opt.num_cpu_threads, "number of CPU threads. At least 2.");
    desc.AddOption("num_output_threads", "", opt.num_output_threads, "number of threads for output. Must be less than num_cpu_threads");
    desc.AddOption("read_lib_file", "", opt.read_lib_file, "input fast[aq] file, can be gzip'ed. \"-\" for stdin.");
    desc.AddOption("assist_seq", "", opt.assist_seq_file, "input assisting fast[aq] file (FILE_NAME.info should exist), can be gzip'ed.");
    desc.AddOption("output_prefix", "", opt.output_prefix, "output prefix");
    desc.AddOption("mem_flag", "", opt.mem_flag, "memory options. 0: minimize memory usage; 1: automatically use moderate memory; other: use all available mem specified by '--host_mem'");
    desc.AddOption("need_mercy", "", opt.need_mercy, "to add mercy edges.");

    try {
        desc.Parse(argc, argv);

        if (opt.read_lib_file == "") {
            throw std::logic_error("No input file!");
        }

        if (opt.num_cpu_threads == 0) {
            opt.num_cpu_threads = omp_get_max_threads();
        }

        if (opt.num_output_threads == 0) {
            opt.num_output_threads = std::max(1, opt.num_cpu_threads / 3);
        }

        if (opt.host_mem == 0) {
            throw std::logic_error("Please specify the host memory!");
        }

        if (opt.gpu_mem == 0) {
#ifdef USE_GPU
            size_t free_gpu_mem, total_gpu_mem;
            get_cuda_memory(free_gpu_mem, total_gpu_mem);
            opt.gpu_mem = free_gpu_mem;
#else
            opt.gpu_mem = 0;
#endif
        }

#ifdef USE_GPU
        if (opt.num_output_threads >= opt.num_cpu_threads) {
            throw std::logic_error("Number of output threads must be less than number of CPU threads!");
        }
#endif
    }
    catch (std::exception &e) {
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
    globals.gpu_mem = opt.gpu_mem;
    globals.num_cpu_threads = opt.num_cpu_threads;
    globals.num_output_threads = opt.num_output_threads;
    globals.read_lib_file = opt.read_lib_file;
    globals.assist_seq_file = opt.assist_seq_file;
    globals.output_prefix = opt.output_prefix;
    globals.mem_flag = opt.mem_flag;
    globals.need_mercy = opt.need_mercy;
    globals.cx1.g_ = &globals;

    // stage1
    if (opt.kmer_freq_threshold > 1) {
        globals.cx1.encode_lv1_diff_base_func_ = cx1_read2sdbg::s1::s1_encode_lv1_diff_base;
        globals.cx1.prepare_func_ = cx1_read2sdbg::s1::s1_read_input_prepare;
        globals.cx1.lv0_calc_bucket_size_func_ = cx1_read2sdbg::s1::s1_lv0_calc_bucket_size;
        globals.cx1.init_global_and_set_cx1_func_ = cx1_read2sdbg::s1::s1_init_global_and_set_cx1;
        globals.cx1.lv1_fill_offset_func_ = cx1_read2sdbg::s1::s1_lv1_fill_offset;
        globals.cx1.lv1_sort_and_proc = cx1_read2sdbg::s1::s1_lv1_direct_sort_and_count;
        globals.cx1.lv2_extract_substr_func_ = cx1_read2sdbg::s1::s1_lv2_extract_substr;
        globals.cx1.lv2_sort_func_ = cx1_read2sdbg::s1::s1_lv2_sort;
        globals.cx1.lv2_pre_output_partition_func_ = cx1_read2sdbg::s1::s1_lv2_pre_output_partition;
        globals.cx1.lv2_output_func_ = cx1_read2sdbg::s1::s1_lv2_output;
        globals.cx1.lv2_post_output_func_ = cx1_read2sdbg::s1::s1_lv2_post_output;
        globals.cx1.post_proc_func_ = cx1_read2sdbg::s1::s1_post_proc;
        globals.cx1.run();
    }
    else {
        cx1_read2sdbg::s1::s1_read_input_prepare(globals);
    }

    // stage2
    globals.cx1.encode_lv1_diff_base_func_ = cx1_read2sdbg::s2::s2_encode_lv1_diff_base;
    globals.cx1.prepare_func_ = cx1_read2sdbg::s2::s2_read_mercy_prepare;
    globals.cx1.lv0_calc_bucket_size_func_ = cx1_read2sdbg::s2::s2_lv0_calc_bucket_size;
    globals.cx1.init_global_and_set_cx1_func_ = cx1_read2sdbg::s2::s2_init_global_and_set_cx1;
    globals.cx1.lv1_fill_offset_func_ = cx1_read2sdbg::s2::s2_lv1_fill_offset;
    globals.cx1.lv1_sort_and_proc = cx1_read2sdbg::s2::s2_lv1_direct_sort_and_proc;
    globals.cx1.lv2_extract_substr_func_ = cx1_read2sdbg::s2::s2_lv2_extract_substr;
    globals.cx1.lv2_sort_func_ = cx1_read2sdbg::s2::s2_lv2_sort;
    globals.cx1.lv2_pre_output_partition_func_ = cx1_read2sdbg::s2::s2_lv2_pre_output_partition;
    globals.cx1.lv2_output_func_ = cx1_read2sdbg::s2::s2_lv2_output;
    globals.cx1.lv2_post_output_func_ = cx1_read2sdbg::s2::s2_lv2_post_output;
    globals.cx1.post_proc_func_ = cx1_read2sdbg::s2::s2_post_proc;
    globals.cx1.run();

    return 0;
}

int main_seq2sdbg(int argc, char **argv) {
    AutoMaxRssRecorder recorder;

    OptionsDescription desc;
    seq2sdbg_opt_t opt;

    desc.AddOption("host_mem", "", opt.host_mem, "memory to be used. No more than 95% of the free memory is recommended. 0 for auto detect.");
    desc.AddOption("gpu_mem", "", opt.gpu_mem, "gpu memory to be used. 0 for auto detect.");
    desc.AddOption("kmer_size", "k", opt.kmer_k, "kmer size");
    desc.AddOption("kmer_from", "", opt.kmer_from, "previous k");
    desc.AddOption("num_cpu_threads", "t", opt.num_cpu_threads, "number of CPU threads. At least 2.");
    desc.AddOption("num_output_threads", "", opt.num_output_threads, "number of threads for output. Must be less than num_cpu_threads");
    desc.AddOption("contig", "", opt.contig, "contigs from previous k");
    desc.AddOption("bubble", "", opt.bubble_seq, "bubble sequence from previous k");
    desc.AddOption("addi_contig", "", opt.addi_contig, "additional contigs from previous k");
    desc.AddOption("local_contig", "", opt.local_contig, "local contigs from previous k");
    desc.AddOption("input_prefix", "", opt.input_prefix, "files input_prefix.edges.* output by count module, can be gzip'ed.");
    desc.AddOption("num_edge_files", "", opt.num_edge_files, "the number of files with name input_prefix.edges.*");
    desc.AddOption("output_prefix", "o", opt.output_prefix, "output prefix");
    desc.AddOption("need_mercy", "", opt.need_mercy, "to add mercy edges. The file input_prefix.cand output by count module should exist.");
    desc.AddOption("mem_flag", "", opt.mem_flag, "memory options. 0: minimize memory usage; 1: automatically use moderate memory; other: use all available mem specified by '--host_mem'");

    try {
        desc.Parse(argc, argv);

        if (opt.input_prefix == "" && opt.contig == "" && opt.addi_contig == "") {
            throw std::logic_error("No input files!");
        }

        if (opt.num_cpu_threads == 0) {
            opt.num_cpu_threads = omp_get_max_threads();
        }

        if (opt.num_output_threads == 0) {
            opt.num_output_threads = std::max(1, opt.num_cpu_threads / 3);
        }

        if (opt.num_edge_files == 0) {
            throw std::logic_error("num edges files cannot be 0!");
        }

        if (opt.kmer_k < 9) {
            throw std::logic_error("kmer size must be >= 9!");
        }

        if (opt.host_mem == 0) {
            throw std::logic_error("Please specify the host memory!");
        }

        if (opt.gpu_mem == 0) {
#ifdef USE_GPU
            size_t free_gpu_mem, total_gpu_mem;
            get_cuda_memory(free_gpu_mem, total_gpu_mem);
            opt.gpu_mem = free_gpu_mem;
#else
            opt.gpu_mem = 0;
#endif
        }

#ifdef USE_GPU
        if (opt.num_output_threads >= opt.num_cpu_threads) {
            throw std::logic_error("Number of output threads must be less than number of CPU threads!");
        }
#endif
    }
    catch (std::exception &e) {
        std::cerr << e.what() << std::endl;
        std::cerr << "Usage: sdbg_builder seq2sdbg -k kmer_size --contig contigs.fa [--addi_contig add.fa] [--input_prefix input] -o out" << std::endl;
        std::cerr << "Options:" << std::endl;
        std::cerr << desc << std::endl;
        exit(1);
    }

    cx1_seq2sdbg::seq2sdbg_global_t globals;
    globals.host_mem = opt.host_mem;
    globals.gpu_mem = opt.gpu_mem;
    globals.num_cpu_threads = opt.num_cpu_threads;
    globals.num_output_threads = opt.num_output_threads;
    globals.input_prefix = opt.input_prefix;
    globals.num_edge_files = opt.num_edge_files;
    globals.output_prefix = opt.output_prefix;
    globals.contig = opt.contig;
    globals.bubble_seq = opt.bubble_seq;
    globals.addi_contig = opt.addi_contig;
    globals.local_contig = opt.local_contig;
    globals.mem_flag = opt.mem_flag;
    globals.kmer_k = opt.kmer_k;
    globals.kmer_from = opt.kmer_from;
    globals.need_mercy = opt.need_mercy;

    xlog("Host memory to be used: %lld\n", (long long)globals.host_mem);
    xlog("Number CPU threads: %d\n", globals.num_cpu_threads);

    // set & run cx1
    globals.cx1.g_ = &globals;
    globals.cx1.encode_lv1_diff_base_func_ = cx1_seq2sdbg::encode_lv1_diff_base;
    globals.cx1.prepare_func_ = cx1_seq2sdbg::read_seq_and_prepare;
    globals.cx1.lv0_calc_bucket_size_func_ = cx1_seq2sdbg::lv0_calc_bucket_size;
    globals.cx1.init_global_and_set_cx1_func_ = cx1_seq2sdbg::init_global_and_set_cx1;
    globals.cx1.lv1_fill_offset_func_ = cx1_seq2sdbg::lv1_fill_offset;
    globals.cx1.lv1_sort_and_proc = cx1_seq2sdbg::lv1_direct_sort_and_proc;
    globals.cx1.lv2_extract_substr_func_ = cx1_seq2sdbg::lv2_extract_substr;
    globals.cx1.lv2_sort_func_ = cx1_seq2sdbg::lv2_sort;
    globals.cx1.lv2_pre_output_partition_func_ = cx1_seq2sdbg::lv2_pre_output_partition;
    globals.cx1.lv2_output_func_ = cx1_seq2sdbg::lv2_output;
    globals.cx1.lv2_post_output_func_ = cx1_seq2sdbg::lv2_post_output;
    globals.cx1.post_proc_func_ = cx1_seq2sdbg::post_proc;

    globals.cx1.run();
    return 0;
}

void DisplayHelp(char *program_name) {
    fprintf(stderr, "Usage: \n");
    fprintf(stderr, "       %s count          kmer counting\n", program_name);
    fprintf(stderr, "       %s read2sdbg      build sdbg from reads\n", program_name);
    fprintf(stderr, "       %s seq2sdbg       build sdbg from megahit contigs + edges\n", program_name);
    fprintf(stderr, "       %s dumpversion       dump version\n", program_name);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        DisplayHelp(argv[0]);
        exit(1);
    }

    if (std::string(argv[1]) == "count")
        return main_kmer_count(argc - 1, argv + 1);

    if (std::string(argv[1]) == "read2sdbg")
        return main_read2sdbg(argc - 1, argv + 1);

    if (std::string(argv[1]) == "seq2sdbg")
        return main_seq2sdbg(argc - 1, argv + 1);
    else if (strcmp(argv[1], "dumpversion") == 0) {
        printf("%s\n", PACKAGE_VERSION);
        return 0;
    }
    else {
        DisplayHelp(argv[0]);
        exit(1);
    }

    return 0;
}