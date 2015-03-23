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

#include <assert.h>
#include <omp.h>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>

#include "options_description.h"
#include "lv2_gpu_functions.h"
#include "helper_functions-inl.h"
#include "sdbg_builder_util.h"
#include "timer.h"

struct Options {
    int kmer_k;
    int min_edge_freq;
    double host_mem;
    double gpu_mem;
    int max_read_length;
    int num_cpu_threads;
    int num_output_threads;
    std::string input_file;
    std::string output_prefix;
    int mem_flag;
    bool need_mercy;

    Options() {
        kmer_k = 21;
        min_edge_freq = 2;
        host_mem = 0;
        gpu_mem = 0;
        max_read_length = 120;
        num_cpu_threads = 0;
        num_output_threads = 0;
        input_file = "";
        output_prefix = "out";
        mem_flag = 1;
        need_mercy = false;
    }
} opt;

void ParsePhase1Option(int argc, char *argv[]) {
    OptionsDescription desc;

    desc.AddOption("kmer_k", "k", opt.kmer_k, "kmer size");
    desc.AddOption("min_kmer_frequency", "m", opt.min_edge_freq, "min frequency to output an edge");
    desc.AddOption("host_mem", "", opt.host_mem, "Max memory to be used. 90% of the free memory is recommended.");
    desc.AddOption("gpu_mem", "", opt.gpu_mem, "gpu memory to be used. 0 for auto detect.");
    desc.AddOption("max_read_length", "", opt.max_read_length, "max read length");
    desc.AddOption("num_cpu_threads", "", opt.num_cpu_threads, "number of CPU threads. At least 2.");
    desc.AddOption("num_output_threads", "", opt.num_output_threads, "number of threads for output. Must be less than num_cpu_threads");
    desc.AddOption("input_file", "", opt.input_file, "input fastx file, can be gzip'ed. \"-\" for stdin.");
    desc.AddOption("output_prefix", "", opt.output_prefix, "output prefix");
    desc.AddOption("mem_flag", "", opt.mem_flag, "memory options. 0: minimize memory usage; 1: automatically use moderate memory; other: use all available mem specified by '--host_mem'");
    desc.AddOption("need_mercy", "", opt.need_mercy, "to add mercy edges.");

    try {
        desc.Parse(argc, argv);
        if (opt.input_file == "") {
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
#ifndef DISABLE_GPU
            size_t free_gpu_mem, total_gpu_mem;
            get_cuda_memory(free_gpu_mem, total_gpu_mem);
            opt.gpu_mem = free_gpu_mem;
#else 
            // we "simulate" the GTX680 here
            opt.gpu_mem = 4243689472ULL;
#endif
        }

        if (opt.num_cpu_threads == 1) {
            throw std::logic_error("Number of CPU threads is at least 2!");
        }
        if (opt.num_output_threads >= opt.num_cpu_threads) {
            throw std::logic_error("Number of output threads must be less than number of CPU threads!");
        }
    } catch (std::exception &e) {
        std::cerr << e.what() << std::endl;
        std::cerr << "Usage: builder count --input_file fastx_file -o out" << std::endl;
        std::cerr << "Options:" << std::endl;
        std::cerr << desc << std::endl;
        exit(1);
    }
}

static AutoMaxRssRecorder recorder;

int main(int argc, char** argv) {
    struct global_data_t globals;

    ParsePhase1Option(argc, argv);

    globals.kmer_k = opt.kmer_k;
    globals.kmer_freq_threshold = opt.min_edge_freq;
    globals.max_read_length = opt.max_read_length;
    globals.host_mem = opt.host_mem;
    globals.gpu_mem = opt.gpu_mem;
    globals.num_cpu_threads = opt.num_cpu_threads;
    globals.phase1_num_output_threads = opt.num_output_threads;
    globals.phase2_num_output_threads = opt.num_output_threads;
    globals.input_file = opt.input_file.c_str();
    globals.output_prefix = opt.output_prefix.c_str();
    globals.mem_flag = opt.mem_flag;
    globals.need_mercy = opt.need_mercy;

    log ("Host memory to be used: %ld\n", globals.host_mem);
    log ("Number CPU threads: %d\n", globals.num_cpu_threads);

    phase1::Phase1Entry(globals);
    phase2::Phase2Entry(globals);

    return 0;
}
