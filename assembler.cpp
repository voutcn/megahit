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
#include <map>
#include <stdexcept>

#include "succinct_dbg.h"
#include "assembly_algorithms.h"
#include "utils.h"
#include "options_description.h"
#include "mem_file_checker-inl.h"
#include "unitig_graph.h"

using std::string;

struct AssemblerOptions {
    string sdbg_name;
    string output_prefix;
    string final_contig_file_name;
    int num_cpu_threads;

    int max_tip_len;
    int min_final_len;
    double min_depth;
    bool is_final_round;
    bool no_bubble;
    int merge_level;
    int prune_level;
    bool excessive_prune;
    double local_low_ratio;
    bool need_fastg;

    AssemblerOptions() {
        output_prefix = "out";
        num_cpu_threads = 0;
        max_tip_len = -1;
        min_final_len = 200;
        no_bubble = false;
        merge_level = 20;
        prune_level = 2;
        local_low_ratio = 0.2;
        min_depth = 1.5;
        is_final_round = false;
    }

    string contig_file() {
        return output_prefix + ".contigs.fa";
    }

    string final_contig_file() {
        return output_prefix + ".final.contigs.fa";
    }

    string addi_contig_file() {
        return output_prefix + ".addi.fa";
    }

} opt;

void ParseOption(int argc, char *argv[]) {
    OptionsDescription desc;

    desc.AddOption("sdbg_name", "s", opt.sdbg_name, "succinct de Bruijn graph name");
    desc.AddOption("output_prefix", "o", opt.output_prefix, "output prefix");
    desc.AddOption("num_cpu_threads", "t", opt.num_cpu_threads, "number of cpu threads");
    desc.AddOption("max_tip_len", "", opt.max_tip_len, "max length for tips to be removed. -1 for 2k");
    desc.AddOption("min_final_len", "", opt.min_final_len, "min length of a final contig");
    desc.AddOption("no_bubble", "", opt.no_bubble, "do not remove bubbles");
    desc.AddOption("merge_level", "", opt.merge_level, "strength of merging");
    desc.AddOption("prune_level", "", opt.prune_level, "strength of low local depth contig pruning (0-2)");
    desc.AddOption("local_low_ratio", "", opt.local_low_ratio, "ratio to define low depth contigs");
    desc.AddOption("min_depth", "", opt.min_depth, "if prune_level is 2, permanently remove low local coverage unitigs under this threshold");
    desc.AddOption("is_final_round", "", opt.is_final_round, "this is the last iteration");

    try {
        desc.Parse(argc, argv);
        if (opt.sdbg_name == "") {
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

void PrintStat(std::map<int64_t, int> &hist) {
    // total length
    int64_t total_length = 0;
    int64_t total_contigs = 0;
    int64_t average_length = 0;
    for (auto it = hist.begin(); it != hist.end(); ++it) {
        total_length += it->first * it->second;
        total_contigs += it->second;
    }

    if (total_contigs > 0) {
        average_length = total_length / total_contigs;
    }

    // N50
    int64_t n50 = -1;
    int64_t acc_length = 0;
    for (auto it = hist.rbegin(); it != hist.rend(); ++it) {
        acc_length += it->first * it->second;
        if (n50 == -1 && acc_length * 2 >= total_length) {
            n50 = it->first;
            break;
        }
    }

    printf("Total length: %lld, N50: %lld, Mean: %lld, number of contigs: %lld\n", (long long)total_length, (long long)n50, (long long)average_length, (long long)total_contigs);
    printf("Maximum length: %llu\n", (unsigned long long)(hist.size() > 0 ? hist.rbegin()->first : 0));
}

static AutoMaxRssRecorder recorder;

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
        printf("Loading succinct de Bruijn graph: %s\n", opt.sdbg_name.c_str());
        dbg.LoadFromFile(opt.sdbg_name.c_str());
        timer.stop();
        printf("Done. Time elapsed: %lf\n", timer.elapsed());
        printf("Number of Edges: %lld\n", (long long)dbg.size);;
        printf("K value: %d\n", dbg.kmer_k);
    }

    { // set parameters
        if (opt.num_cpu_threads == 0) {
            opt.num_cpu_threads = omp_get_max_threads();
        }
        omp_set_num_threads(opt.num_cpu_threads);
        printf("Number of CPU threads: %d\n", opt.num_cpu_threads);

        if (opt.max_tip_len == -1) {
            opt.max_tip_len = dbg.kmer_k * 2;
        }
    }

    if (opt.max_tip_len > 0) { // tips removal
        timer.reset();
        timer.start();
        assembly_algorithms::RemoveTips(dbg, opt.max_tip_len, opt.min_final_len);
        timer.stop();
        printf("Tips removal done! Time elapsed(sec): %lf\n", timer.elapsed());
    }

    // construct unitig graph
    timer.reset();
    timer.start();
    UnitigGraph unitig_graph(&dbg);
    unitig_graph.InitFromSdBG();
    timer.stop();
    printf("unitig graph size: %u, time for building: %lf\n", unitig_graph.size(), timer.elapsed());

    // remove bubbles
    if (!opt.no_bubble) {
        timer.reset();
        timer.start();
        unitig_graph.MergeBubbles(true);
        if (opt.merge_level > 0) {
            unitig_graph.MergeComplexBubbles(0.98, opt.merge_level, true);
        }
        timer.stop();
        printf("Time elapsed(sec): %lf\n", timer.elapsed());
    }

    // excessive pruning
    static const int kLocalWidth = 1000;
    int64_t num_removed = 0;

    if (opt.prune_level >= 2) {
        timer.reset();
        timer.start();
        unitig_graph.RemoveLocalLowDepth(opt.min_depth, opt.max_tip_len, kLocalWidth, std::min(opt.local_low_ratio, 0.1), num_removed, true);
        timer.stop();
        printf("Unitigs removed in excessive pruning: %lld, time: %lf\n", (long long)num_removed, timer.elapsed());   
    }

    // output contigs
    std::map<int64_t, int> histogram;
    FILE *out_contig_file = OpenFileAndCheck(opt.contig_file().c_str(), "w");
    FILE *out_final_contig_file = OpenFileAndCheck(opt.final_contig_file().c_str(), "w");

    if (!(opt.is_final_round && opt.prune_level >= 1)) { // otherwise output after local low depth pruning
        timer.reset();
        timer.start();
        histogram.clear();
        out_final_contig_file = NULL; // uncomment to avoid output final contigs

        unitig_graph.OutputContigs(out_contig_file, out_final_contig_file, histogram, false, opt.min_final_len);
        PrintStat(histogram);

        timer.stop();
        printf("Time to output: %lf\n", timer.elapsed());
    }

    // remove local low depth & output additional contigs
    if (opt.prune_level >= 1) {
        FILE *out_addi_contig_file = OpenFileAndCheck(opt.addi_contig_file().c_str(), "w");

        timer.reset();
        timer.start();
        num_removed = 0;
        double min_depth = opt.min_depth;

        while (min_depth < kMaxMulti_t) {
            if (!unitig_graph.RemoveLocalLowDepth(min_depth, opt.max_tip_len, kLocalWidth, opt.local_low_ratio, num_removed, opt.is_final_round)) {
                break;
            }

            min_depth *= 1.1;
        }
        
        if (opt.merge_level > 0)
            unitig_graph.MergeComplexBubbles(0.98, opt.merge_level, opt.is_final_round);

        timer.stop();
        printf("Number of local low depth unitigs removed: %lld, time: %lf\n", (long long)num_removed, timer.elapsed());

        histogram.clear();

        if (!opt.is_final_round) {
            unitig_graph.OutputContigs(out_addi_contig_file, NULL, histogram, true, 0);
        } else {
            unitig_graph.OutputContigs(out_contig_file, out_final_contig_file, histogram, false, opt.min_final_len);
        }

        PrintStat(histogram);

        fclose(out_addi_contig_file);
    }

    fclose(out_contig_file);
    fclose(out_final_contig_file);

    return 0;
}