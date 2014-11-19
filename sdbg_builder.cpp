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
// #include <sys/sysinfo.h>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>

#include "options_description.h"
#include "lv2_gpu_functions.h"
#include "helper_functions-inl.h"
#include "sdbg_builder_util.h"

struct Phase1Options {
    int kmer_k;
    int min_edge_freq;
    double host_mem;
    double gpu_mem;
    int max_read_length;
    int num_cpu_threads;
    int num_output_threads;
    std::string input_file;
    std::string output_prefix;

    Phase1Options() {
        kmer_k = 21;
        min_edge_freq = 2;
        host_mem = 0;
        gpu_mem = 0;
        max_read_length = 120;
        num_cpu_threads = 0;
        num_output_threads = 0;
        input_file = "";
        output_prefix = "out";
    }
} phase1_options;

struct Phase2Options {
    bool need_mercy;
    double host_mem;
    double gpu_mem;
    int num_edge_files;
    int num_cpu_threads;
    int num_output_threads;
    int max_read_length;
    std::string input_prefix;
    std::string output_prefix;

    Phase2Options() {
        need_mercy = false;
        host_mem = 0;
        gpu_mem = 0;
        num_edge_files = 0;
        num_cpu_threads = 0;
        num_output_threads = 0;
        max_read_length = 120;

        input_prefix = "";
        output_prefix = "out";
    }
} phase2_options;

void ParsePhase1Option(int argc, char *argv[]) {
    OptionsDescription desc;

    desc.AddOption("kmer_k", "k", phase1_options.kmer_k, "kmer size");
    desc.AddOption("min_kmer_frequency", "m", phase1_options.min_edge_freq, "min frequency to output an edge");
    desc.AddOption("host_mem", "", phase1_options.host_mem, "memory to be used. No more than 95% of the free memory is recommended. 0 for auto detect.");
    desc.AddOption("gpu_mem", "", phase1_options.gpu_mem, "gpu memory to be used. 0 for auto detect.");
    desc.AddOption("max_read_length", "", phase1_options.max_read_length, "max read length");
    desc.AddOption("num_cpu_threads", "", phase1_options.num_cpu_threads, "number of CPU threads. At least 2.");
    desc.AddOption("num_output_threads", "", phase1_options.num_output_threads, "number of threads for output. Must be less than num_cpu_threads");
    desc.AddOption("input_file", "", phase1_options.input_file, "input fastx file, can be gzip'ed. \"-\" for stdin.");
    desc.AddOption("output_prefix", "", phase1_options.output_prefix, "output prefix");

    try {
        desc.Parse(argc, argv);
        if (phase1_options.input_file == "") {
            throw std::logic_error("No input file!");
        }

        if (phase1_options.num_cpu_threads == 0) {
            phase1_options.num_cpu_threads = omp_get_max_threads();
        }

        if (phase1_options.num_output_threads == 0) {
            phase1_options.num_output_threads = std::max(1, phase1_options.num_cpu_threads / 3);
        }

        if (phase1_options.host_mem == 0) {
            throw std::logic_error("Please specify the host memory!");
            // struct sysinfo s_info;
            // sysinfo(&s_info);
            // phase1_options.host_mem = (s_info.freeram + s_info.bufferram) * 0.95;
        }

        if (phase1_options.gpu_mem == 0) {
#ifndef DISABLE_GPU
            size_t free_gpu_mem, total_gpu_mem;
            get_cuda_memory(free_gpu_mem, total_gpu_mem);
            phase1_options.gpu_mem = free_gpu_mem;
#else 
            // we "simulate" the GTX680 here
            phase1_options.gpu_mem = 4243689472ULL;
#endif
        }

        if (phase1_options.num_cpu_threads == 1) {
            throw std::logic_error("Number of CPU threads is at least 2!");
        }
        if (phase1_options.num_output_threads >= phase1_options.num_cpu_threads) {
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

void ParsePhase2Option(int argc, char *argv[]) {
    OptionsDescription desc;

    desc.AddOption("host_mem", "", phase2_options.host_mem, "memory to be used. No more than 95% of the free memory is recommended. 0 for auto detect.");
    desc.AddOption("gpu_mem", "", phase2_options.gpu_mem, "gpu memory to be used. 0 for auto detect.");
    desc.AddOption("num_cpu_threads", "t", phase2_options.num_cpu_threads, "number of CPU threads. At least 2.");
    desc.AddOption("num_output_threads", "", phase2_options.num_output_threads, "number of threads for output. Must be less than num_cpu_threads");
    desc.AddOption("input_prefix", "", phase2_options.input_prefix, "files input_prefix.edges.* output by count module, can be gzip'ed.");
    desc.AddOption("num_edge_files", "", phase2_options.num_edge_files, "the number of files with name input_prefix.edges.*");
    desc.AddOption("output_prefix", "o", phase2_options.output_prefix, "output prefix");
    desc.AddOption("need_mercy", "", phase2_options.need_mercy, "to add mercy edges. The file input_prefix.cand output by count module should exist.");
    desc.AddOption("max_read_length", "", phase2_options.max_read_length, "max read length");

    try {
        desc.Parse(argc, argv);
        if (phase2_options.input_prefix == "") {
            throw std::logic_error("No input prefix!");
        }
        if (phase2_options.num_edge_files == 0) {
            throw std::logic_error("Number of edge files cannot be 0!");
        }

        if (phase2_options.num_cpu_threads == 0) {
            phase2_options.num_cpu_threads = omp_get_max_threads();
        }

        if (phase2_options.num_output_threads == 0) {
            phase2_options.num_output_threads = std::max(1, phase2_options.num_cpu_threads / 3);
        }

        if (phase2_options.host_mem == 0) {
            throw std::logic_error("Please specify the host memory!");
            // struct sysinfo s_info;
            // sysinfo(&s_info);
            // phase2_options.host_mem = (s_info.freeram + s_info.bufferram) * 0.95;
        }

        if (phase2_options.gpu_mem == 0) {
#ifndef DISABLE_GPU
            size_t free_gpu_mem, total_gpu_mem;
            get_cuda_memory(free_gpu_mem, total_gpu_mem);
            phase2_options.gpu_mem = free_gpu_mem;
#else 
            // we "simulate" the GTX680 here
            phase2_options.gpu_mem = 4243689472ULL;
#endif
        }

        if (phase2_options.num_cpu_threads == 1) {
            throw std::logic_error("Number of CPU threads is at least 2!");
        }
        if (phase2_options.num_output_threads >= phase2_options.num_cpu_threads) {
            throw std::logic_error("Number of output threads must be less than number of CPU threads!");
        }
    } catch (std::exception &e) {
        std::cerr << e.what() << std::endl;
        std::cerr << "Usage: builder build --input_prefix input --num_edge_files num -o out" << std::endl;
        std::cerr << "Options:" << std::endl;
        std::cerr << desc << std::endl;
        exit(1);
    }
}

void DisplayHelp(char *program_name) {
    fprintf(stderr, "Usage: \n");
    fprintf(stderr, "    1. Counting & output solid edges: \n");
    fprintf(stderr, "       type \"%s count\" for help.\n", program_name);
    fprintf(stderr, "    2. Build Succinct dBG from solid edges: \n");
    fprintf(stderr, "       type \"%s build\" for help.\n", program_name);
}

int main(int argc, char** argv) {
    struct global_data_t globals;

    if (argc < 2 || 
        (std::string(argv[1]) != "count" && std::string(argv[1]) != "build") ) {
        DisplayHelp(argv[0]);
        exit(1);
    }

    if (std::string(argv[1]) == "count") {
        ParsePhase1Option(argc - 1, argv + 1);

        globals.kmer_k = phase1_options.kmer_k;
        globals.kmer_freq_threshold = phase1_options.min_edge_freq;
        globals.max_read_length = phase1_options.max_read_length;
        globals.host_mem = phase1_options.host_mem;
        globals.gpu_mem = phase1_options.gpu_mem;
        globals.num_cpu_threads = phase1_options.num_cpu_threads;
        globals.phase1_num_output_threads = phase1_options.num_output_threads;
        globals.input_file = phase1_options.input_file.c_str();
        globals.output_prefix = phase1_options.output_prefix.c_str();

        log ("Host memory to be used: %ld\n", globals.host_mem);
        log ("Number CPU threads: %d\n", globals.num_cpu_threads);

#ifndef DISABLE_GPU
        log ("GPU memory to be used: %ld\n",  globals.gpu_mem);
        if (globals.gpu_mem < 2147483648LL) {
            err("Warning, maybe not enough GPU memory. At least 2G is recommended. Process will be continue though.\n");
        }
#endif

        phase1::Phase1Entry(globals);
    } else if (std::string(argv[1]) == "build") {
        ParsePhase2Option(argc - 1, argv + 1);

        globals.need_mercy = phase2_options.need_mercy;
        globals.host_mem = phase2_options.host_mem;
        globals.gpu_mem = phase2_options.gpu_mem;
        globals.num_cpu_threads = phase2_options.num_cpu_threads;
        globals.phase1_num_output_threads = phase2_options.num_edge_files;
        globals.phase2_num_output_threads = phase2_options.num_output_threads;
        globals.phase2_input_prefix = phase2_options.input_prefix.c_str();
        globals.output_prefix = phase2_options.output_prefix.c_str();
        globals.max_read_length = phase2_options.max_read_length;

        log ("Host memory to be used: %ld\n", globals.host_mem);
        log ("Number CPU threads: %d\n", globals.num_cpu_threads);

#ifndef DISABLE_GPU
        log ("GPU memory to be used: %ld\n",  globals.gpu_mem);
        if (globals.gpu_mem < 2147483648LL) {
            err("Warning, maybe not enough GPU memory. At least 2G is recommended. Process will be continue though.\n");
        }
#endif
        phase2::Phase2Entry(globals);
    }

    return 0;
}
