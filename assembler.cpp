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

#include <omp.h>
#include <assert.h>
#include <string>
#include <algorithm>
#include <iostream>
#include <stdexcept>

#include "succinct_dbg.h"
#include "assembly_algorithms.h"
#include "timer.h"
#include "options_description.h"
#include "mem_file_checker-inl.h"

using std::string;

struct AssemblerOptions {
    string sdbg_name;
    string output_prefix;
    string final_contig_file_name;
    int num_cpu_threads;

    int max_tip_len;
    int min_final_contig_len;
    bool is_final_round;
    bool no_bubble;
    double bubble_remove_ratio;
    bool remove_low_local;
    double low_local_ratio;
    int mercy_threshold;

    AssemblerOptions() {
        output_prefix = "out";
        num_cpu_threads = 0;
        max_tip_len = -1;
        min_final_contig_len = 200;
        no_bubble = false;
        bubble_remove_ratio = 1;
        remove_low_local = false;
        low_local_ratio = 0.2;
        is_final_round = false;
    }

    string contig_file() {
        return output_prefix + ".contigs.fa";
    }

    string multi_file() {
        return output_prefix + ".multi";
    }

    string final_contig_file() {
        return output_prefix + ".final.contigs.fa";
    }

    string addi_contig_file() {
        return output_prefix + ".addi.fa";
    }

    string addi_multi_file() {
        return output_prefix + ".addi.multi";
    }

} options;

void ParseOption(int argc, char *argv[]) {
    OptionsDescription desc;

    desc.AddOption("sdbg_name", "s", options.sdbg_name, "succinct de Bruijn graph name");
    desc.AddOption("output_prefix", "o", options.output_prefix, "output prefix");
    desc.AddOption("num_cpu_threads", "t", options.num_cpu_threads, "number of cpu threads");
    desc.AddOption("max_tip_len", "", options.max_tip_len, "max length for tips to be removed. -1 for 2k");
    desc.AddOption("min_final_contig_len", "", options.min_final_contig_len, "min length to output a final contig");
    desc.AddOption("no_bubble", "", options.no_bubble, "do not remove bubbles");
    desc.AddOption("bubble_remove_ratio", "", options.bubble_remove_ratio, "bubbles with multiplicities lower than this ratio times to highest of its group will be removed");
    desc.AddOption("remove_low_local", "", options.remove_low_local, "remove low local depth contigs progressively");
    desc.AddOption("low_local_ratio", "", options.low_local_ratio, "ratio to define low depth contigs");
    desc.AddOption("is_final_round", "", options.is_final_round, "this is the last iteration");

    try {
        desc.Parse(argc, argv);
        if (options.sdbg_name == "") {
            throw std::logic_error("no succinct de Bruijn graph name!");
        }
    } catch (std::exception &e) {
        std::cerr << e.what() << std::endl;
        std::cerr << "Usage: " << argv[0] << " -s sdbg_name -o output_prefix" << std::endl;
        std::cerr << "options:" << std::endl;
        std::cerr << desc << std::endl;
        exit(1);
    }
}

int main(int argc, char **argv) {
    // set stdout line buffered
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    ParseOption(argc, argv);

    SuccinctDBG dbg;  
    xtimer_t timer;

    { // graph loading
        timer.reset();
        timer.start();
        printf("Loading succinct de Bruijn graph: %s\n", options.sdbg_name.c_str());
        dbg.LoadFromFile(options.sdbg_name.c_str());
        timer.stop();
        printf("Done. Time elapsed: %lf\n", timer.elapsed());
        printf("Number of Edges: %lld\n", (long long)dbg.size);;
        printf("K value: %d\n", dbg.kmer_k);
    }

    { // set parameters
        if (options.num_cpu_threads == 0) {
            options.num_cpu_threads = omp_get_max_threads();
        }
        omp_set_num_threads(options.num_cpu_threads);
        printf("Number of CPU threads: %d\n", options.num_cpu_threads);

        if (options.max_tip_len == -1) {
            options.max_tip_len = dbg.kmer_k * 2;
        }
    }

    if (options.max_tip_len > 0) { // tips removal
        timer.reset();
        timer.start();
        assembly_algorithms::RemoveTips(dbg, options.max_tip_len, options.min_final_contig_len);
        timer.stop();
        printf("Tips removal done! Time elapsed(sec): %lf\n", timer.elapsed());
    }

    if (!options.no_bubble) { // merge bubbles
        timer.reset();
        timer.start();
        int64_t num_bubbles = assembly_algorithms::PopBubbles(dbg, dbg.kmer_k + 2, options.bubble_remove_ratio);
        timer.stop();
        printf("Number of bubbles: %lld. Time elapsed: %lf\n", (long long)num_bubbles, timer.elapsed());
    }


    FILE *out_contig_file = OpenFileAndCheck(options.contig_file().c_str(), "w");
    FILE *out_multi_file = OpenFileAndCheck(options.multi_file().c_str(), "wb");
    FILE *out_final_contig_file = OpenFileAndCheck(options.final_contig_file().c_str(), "w");
    assert(out_contig_file != NULL);
    assert(out_multi_file != NULL);
    assert(out_final_contig_file != NULL);

    if (options.remove_low_local) { // remove local low depth
        timer.reset();
        timer.start();

        printf("Removing low local coverage...\n");
        if (!options.is_final_round) {
            FILE *out_addi_contig_file = OpenFileAndCheck(options.addi_contig_file().c_str(), "w");
            FILE *out_addi_multi_file = OpenFileAndCheck(options.addi_multi_file().c_str(), "w");
            assert(out_addi_multi_file != NULL);
            assert(out_addi_contig_file != NULL);

            // FILE *out_final_contig_file = NULL; // uncomment to avoid output final contigs
            assembly_algorithms::RemoveLowLocalAndOutputChanged(
                dbg, out_contig_file, 
                out_multi_file, 
                out_final_contig_file,
                out_addi_contig_file,
                out_addi_multi_file, 
                2, 
                dbg.kmer_k * 2, 
                options.low_local_ratio,
                options.min_final_contig_len);
            
            fclose (out_addi_contig_file);
            fclose(out_addi_multi_file);
        } else {
            assembly_algorithms::RemoveLowLocalAndOutputFinal(
                dbg, 
                out_final_contig_file, 
                2, 
                dbg.kmer_k * 2, 
                options.low_local_ratio,
                options.min_final_contig_len);
        }
        timer.stop();
        printf("Done! Time elapsed(sec.): %lf\n", timer.elapsed());
    } else {
        printf("Assembly after bubble merging...\n");
        timer.reset();
        timer.start();
        if (!options.is_final_round) {
            // FILE *out_final_contig_file = NULL; // uncomment to avoid output final contigs
            assembly_algorithms::AssembleFromUnitigGraph(
                dbg, 
                out_contig_file, 
                out_multi_file,
                out_final_contig_file,
                options.min_final_contig_len);
        } else {
            assembly_algorithms::AssembleFinalFromUnitigGraph(
                dbg, 
                out_final_contig_file, 
                options.min_final_contig_len);
        }

        timer.stop();
        printf("Done! Time elapsed(sec.): %lf\n", timer.elapsed());
    }

    fclose(out_contig_file);
    fclose(out_multi_file);
    fclose(out_final_contig_file);

    return 0;
}