/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong
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

struct asm_opt_t {
    string sdbg_name;
    string output_prefix;
    int num_cpu_threads;

    int max_tip_len;
    int min_standalone;
    double min_depth;
    bool is_final_round;
    bool no_bubble;
    int merge_len;
    double merge_similar;
    int prune_level;
    double low_local_ratio;
    bool output_standalone;

    asm_opt_t() {
        output_prefix = "out";
        num_cpu_threads = 0;
        max_tip_len = -1;
        min_standalone = 200;
        no_bubble = false;
        merge_len = 20;
        merge_similar = 0.98;
        prune_level = 2;
        low_local_ratio = 0.2;
        min_depth = 1.5;
        is_final_round = false;
        output_standalone = false;
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

};

static asm_opt_t opt;

void ParseAsmOption(int argc, char *argv[]) {
    OptionsDescription desc;

    desc.AddOption("sdbg_name", "s", opt.sdbg_name, "succinct de Bruijn graph name");
    desc.AddOption("output_prefix", "o", opt.output_prefix, "output prefix");
    desc.AddOption("num_cpu_threads", "t", opt.num_cpu_threads, "number of cpu threads");
    desc.AddOption("max_tip_len", "", opt.max_tip_len, "max length for tips to be removed. -1 for 2k");
    desc.AddOption("min_standalone", "", opt.min_standalone, "min length of a standalone contig to output to final.contigs.fa");
    desc.AddOption("no_bubble", "", opt.no_bubble, "do not remove bubbles");
    desc.AddOption("merge_len", "", opt.merge_len, "merge complex bubbles of length <= merge_len * k");
    desc.AddOption("merge_similar", "", opt.merge_similar, "min similarity of complex bubble merging");
    desc.AddOption("prune_level", "", opt.prune_level, "strength of low local depth contig pruning (0-2)");
    desc.AddOption("low_local_ratio", "", opt.low_local_ratio, "ratio to define low depth contigs");
    desc.AddOption("min_depth", "", opt.min_depth, "if prune_level is 2, permanently remove low local coverage unitigs under this threshold");
    desc.AddOption("is_final_round", "", opt.is_final_round, "this is the last iteration");
    desc.AddOption("output_standalone", "", opt.output_standalone, "output standalone contigs to *.final.contigs.fa");

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

    xlog("Total length: %lld, N50: %lld, Mean: %lld, number of contigs: %lld\n", (long long)total_length, (long long)n50, (long long)average_length, (long long)total_contigs);
    xlog("Maximum length: %llu\n", (unsigned long long)(hist.size() > 0 ? hist.rbegin()->first : 0));
}

int main_assemble(int argc, char **argv) {
    ParseAsmOption(argc, argv);

    SuccinctDBG dbg;
    xtimer_t timer;

    {
        // graph loading
        timer.reset();
        timer.start();
        xlog("Loading succinct de Bruijn graph: %s ", opt.sdbg_name.c_str());
        dbg.LoadFromFile(opt.sdbg_name.c_str());
        timer.stop();
        xlog_ext("Done. Time elapsed: %lf\n", timer.elapsed());
        xlog("Number of Edges: %lld; K value: %d\n", (long long)dbg.size, dbg.kmer_k);
    }

    {
        // set parameters
        if (opt.num_cpu_threads == 0) {
            opt.num_cpu_threads = omp_get_max_threads();
        }
        omp_set_num_threads(opt.num_cpu_threads);
        xlog("Number of CPU threads: %d\n", opt.num_cpu_threads);

        if (opt.max_tip_len == -1) {
            opt.max_tip_len = dbg.kmer_k * 2;
        }
    }

    if (opt.max_tip_len > 0) { // tips removal
        timer.reset();
        timer.start();
        assembly_algorithms::RemoveTips(dbg, opt.max_tip_len, opt.min_standalone);
        timer.stop();
        xlog("Tips removal done! Time elapsed(sec): %lf\n", timer.elapsed());
    }

    // construct unitig graph
    timer.reset();
    timer.start();
    UnitigGraph unitig_graph(&dbg);
    unitig_graph.InitFromSdBG();
    timer.stop();
    xlog("unitig graph size: %u, time for building: %lf\n", unitig_graph.size(), timer.elapsed());

    // remove bubbles
    if (!opt.no_bubble) {
        timer.reset();
        timer.start();
        uint32_t num_bubbles = unitig_graph.MergeBubbles(true);
        uint32_t num_complex_bubbles = 0;
        if (opt.merge_len > 0) {
            num_complex_bubbles += unitig_graph.MergeComplexBubbles(opt.merge_similar, opt.merge_len, true);
        }
        timer.stop();
        xlog("Number of bubbles/complex bubbles removed: %u/%u, Time elapsed(sec): %lf\n",
             num_bubbles, num_complex_bubbles, timer.elapsed());
    }

    // excessive pruning
    static const int kLocalWidth = 1000;
    int64_t num_removed = 0;

    if (opt.prune_level >= 2) {
        timer.reset();
        timer.start();
        unitig_graph.RemoveLocalLowDepth(opt.min_depth, opt.max_tip_len, kLocalWidth, std::min(opt.low_local_ratio, 0.1), num_removed, true);
        timer.stop();
        xlog("Unitigs removed in excessive pruning: %lld, time: %lf\n", (long long)num_removed, timer.elapsed());
    }

    // output contigs
    std::map<int64_t, int> histogram;
    FILE *out_contig_file = OpenFileAndCheck(opt.contig_file().c_str(), "w");
    FILE *out_final_contig_file = OpenFileAndCheck(opt.final_contig_file().c_str(), "w");

    if (!(opt.is_final_round && opt.prune_level >= 1)) { // otherwise output after local low depth pruning
        timer.reset();
        timer.start();
        histogram.clear();

        // unitig_graph.OutputContigs(out_contig_file, out_final_contig_file, histogram, false, opt.min_standalone);
        unitig_graph.OutputContigs(out_contig_file, opt.output_standalone ? out_final_contig_file : NULL,
                                   histogram, false, opt.min_standalone);

        PrintStat(histogram);

        timer.stop();
        xlog("Time to output: %lf\n", timer.elapsed());
    }

    // remove local low depth & output additional contigs
    if (opt.prune_level >= 1) {
        FILE *out_addi_contig_file = OpenFileAndCheck(opt.addi_contig_file().c_str(), "w");

        timer.reset();
        timer.start();
        num_removed = 0;
        double min_depth = opt.min_depth;

        while (min_depth < kMaxMulti_t) {
            if (!unitig_graph.RemoveLocalLowDepth(min_depth, opt.max_tip_len, kLocalWidth, opt.low_local_ratio, num_removed, opt.is_final_round)) {
                break;
            }

            min_depth *= 1.1;
        }

        uint32_t num_complex_bubbles = 0;
        if (opt.merge_len > 0)
            num_complex_bubbles = unitig_graph.MergeComplexBubbles(opt.merge_similar, opt.merge_len, opt.is_final_round);

        timer.stop();
        xlog("Number of local low depth unitigs removed: %lld, complex bubbles removed: %u, time: %lf\n",
             (long long)num_removed, num_complex_bubbles, timer.elapsed());

        histogram.clear();

        if (!opt.is_final_round) {
            unitig_graph.OutputContigs(out_addi_contig_file, NULL, histogram, true, 0);
        } else {
            // unitig_graph.OutputContigs(out_contig_file, out_final_contig_file, histogram, false, opt.min_standalone);
            unitig_graph.OutputContigs(out_contig_file, opt.output_standalone ? out_final_contig_file : NULL, histogram, false, opt.min_standalone);
        }

        PrintStat(histogram);
        fclose(out_addi_contig_file);
    }

    fclose(out_contig_file);
    fclose(out_final_contig_file);

    return 0;
}