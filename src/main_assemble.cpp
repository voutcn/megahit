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
#include <algorithm>
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <string>

#include "assembly/all_algo.h"
#include "assembly/contig_stat.h"
#include "assembly/contig_writer.h"
#include "utils/histgram.h"
#include "utils/options_description.h"
#include "utils/safe_open.h"
#include "utils/utils.h"

using std::string;

namespace {

struct Option {
  string sdbg_name;
  string output_prefix{"out"};
  int num_cpu_threads{0};

  int local_width{1000};
  int max_tip_len{-1};
  int min_standalone{200};
  double min_depth{-1};
  bool is_final_round{false};
  int bubble_level{2};
  int merge_len{20};
  double merge_similar{0.98};
  int prune_level{2};
  double low_local_ratio{0.2};
  bool output_standalone{false};
  bool careful_bubble{false};

  string contig_file() { return output_prefix + ".contigs.fa"; }
  string standalone_file() { return output_prefix + ".final.contigs.fa"; }
  string addi_contig_file() { return output_prefix + ".addi.fa"; }
  string bubble_file() { return output_prefix + ".bubble_seq.fa"; }
} opt;

void ParseAsmOption(int argc, char *argv[]) {
  OptionsDescription desc;

  desc.AddOption("sdbg_name", "s", opt.sdbg_name, "succinct de Bruijn graph name");
  desc.AddOption("output_prefix", "o", opt.output_prefix, "output prefix");
  desc.AddOption("num_cpu_threads", "t", opt.num_cpu_threads, "number of cpu threads");
  desc.AddOption("max_tip_len", "", opt.max_tip_len, "max length for tips to be removed. -1 for 2k");
  desc.AddOption("min_standalone", "", opt.min_standalone,
                 "min length of a standalone contig to output to final.contigs.fa");
  desc.AddOption("bubble_level", "", opt.bubble_level, "bubbles level 0-3");
  desc.AddOption("merge_len", "", opt.merge_len, "merge complex bubbles of length <= merge_len * k");
  desc.AddOption("merge_similar", "", opt.merge_similar, "min similarity of complex bubble merging");
  desc.AddOption("prune_level", "", opt.prune_level, "strength of low local depth contig pruning (0-3)");
  desc.AddOption("low_local_ratio", "", opt.low_local_ratio, "ratio to define low depth contigs");
  desc.AddOption("min_depth", "", opt.min_depth,
                 "if prune_level >= 2, permanently remove low local coverage unitigs under this threshold");
  desc.AddOption("is_final_round", "", opt.is_final_round, "this is the last iteration");
  desc.AddOption("output_standalone", "", opt.output_standalone, "output standalone contigs to *.final.contigs.fa");
  desc.AddOption("careful_bubble", "", opt.careful_bubble, "remove bubble carefully");

  try {
    desc.Parse(argc, argv);
    if (opt.sdbg_name.empty()) {
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

}  // namespace

int main_assemble(int argc, char **argv) {
  AutoMaxRssRecorder recorder;
  ParseAsmOption(argc, argv);

  SDBG dbg;
  SimpleTimer timer;

  // graph loading
  timer.reset();
  timer.start();
  xinfo("Loading succinct de Bruijn graph: %s ", opt.sdbg_name.c_str());
  dbg.LoadFromFile(opt.sdbg_name.c_str());
  timer.stop();
  xinfoc("Done. Time elapsed: %lf\n", timer.elapsed());
  xinfo("Number of Edges: %lld; K value: %d\n", (long long)dbg.size(), dbg.k());

  // set cpu threads
  if (opt.num_cpu_threads == 0) {
    opt.num_cpu_threads = omp_get_max_threads();
  }
  omp_set_num_threads(opt.num_cpu_threads);
  xinfo("Number of CPU threads: %d\n", opt.num_cpu_threads);

  // set tip len
  if (opt.max_tip_len == -1) {
    opt.max_tip_len = dbg.k() * 2;
  }
  // set min depth
  if (opt.min_depth <= 0) {
    opt.min_depth = sdbg_pruning::InferMinDepth(dbg);
    xinfo("min depth set to %.3lf\n", opt.min_depth);
  }

  // tips removal before building unitig graph
  if (opt.max_tip_len > 0) {
    timer.reset();
    timer.start();
    sdbg_pruning::RemoveTips(dbg, opt.max_tip_len);
    timer.stop();
    xinfo("Tips removal done! Time elapsed(sec): %lf\n", timer.elapsed());
  }

  // construct unitig graph
  timer.reset();
  timer.start();
  UnitigGraph graph(&dbg);
  timer.stop();
  xinfo("unitig graph size: %u, time for building: %lf\n", graph.size(), timer.elapsed());
  CalcAndPrintStat(graph);

  // set up bubble
  FILE *bubble_file = xfopen(opt.bubble_file().c_str(), "w");
  NaiveBubbleRemover naiver_bubble_remover;
  ComplexBubbleRemover complex_bubble_remover;
  complex_bubble_remover.set_merge_similarity(opt.merge_similar).set_merge_level(opt.merge_len);
  Histgram<int64_t> bubble_hist;
  if (opt.careful_bubble) {
    naiver_bubble_remover.set_careful_threshold(0.2).set_bubble_file(bubble_file).set_hist(bubble_hist);
    complex_bubble_remover.set_careful_threshold(0.2).set_bubble_file(bubble_file).set_hist(bubble_hist);
  }

  // graph cleaning for 5 rounds
  for (int round = 1; round <= 5; ++round) {
    bool changed = false;
    if (round > 1) {
      timer.reset();
      timer.start();
      uint32_t num_tips = RemoveTips(graph, opt.max_tip_len);
      changed |= num_tips > 0;
      timer.stop();
      xinfo("Tips removed: %u, time: %lf\n", num_tips, timer.elapsed());
    }
    // remove bubbles
    if (opt.bubble_level >= 1) {
      timer.reset();
      timer.start();
      uint32_t num_bubbles = naiver_bubble_remover.PopBubbles(graph, true);
      timer.stop();
      xinfo("Number of bubbles removed: %u, Time elapsed(sec): %lf\n", num_bubbles, timer.elapsed());
      changed |= num_bubbles > 0;
    }
    // remove complex bubbles
    if (opt.bubble_level >= 2) {
      timer.reset();
      timer.start();
      uint32_t num_bubbles = complex_bubble_remover.PopBubbles(graph, true);
      timer.stop();
      xinfo("Number of complex bubbles removed: %u, Time elapsed(sec): %lf\n", num_bubbles, timer.elapsed());
      changed |= num_bubbles > 0;
    }

    // disconnect
    timer.reset();
    timer.start();
    uint32_t num_disconnected = DisconnectWeakLinks(graph, 0.1);
    timer.stop();
    xinfo("Number unitigs disconnected: %u, time: %lf\n", num_disconnected, timer.elapsed());
    changed |= num_disconnected > 0;

    // excessive pruning
    uint32_t num_excessive_pruned = 0;
    if (opt.prune_level >= 3) {
      timer.reset();
      timer.start();
      num_excessive_pruned = RemoveLowDepth(graph, opt.min_depth);
      num_excessive_pruned += naiver_bubble_remover.PopBubbles(graph, true);
      if (opt.bubble_level >= 2 && opt.merge_len > 0) {
        num_excessive_pruned += complex_bubble_remover.PopBubbles(graph, true);
      }
      timer.stop();
      xinfo("Unitigs removed in (more-)excessive pruning: %llu, time: %lf\n", num_excessive_pruned, timer.elapsed());
    } else if (opt.prune_level >= 2) {
      timer.reset();
      timer.start();
      RemoveLocalLowDepth(graph, opt.min_depth, opt.max_tip_len, opt.local_width, std::min(opt.low_local_ratio, 0.1),
                          true, &num_excessive_pruned);
      timer.stop();
      xinfo("Unitigs removed in excessive pruning: %llu, time: %lf\n", num_excessive_pruned, timer.elapsed());
    }
    if (!changed) break;
  }

  ContigStat stat = CalcAndPrintStat(graph);

  // output contigs
  FILE *contig_file = xfopen(opt.contig_file().c_str(), "w");
  FILE *standalone_file = xfopen(opt.standalone_file().c_str(), "w");

  if (!(opt.is_final_round && opt.prune_level >= 1)) {  // otherwise output after local low depth pruning
    timer.reset();
    timer.start();
    FILE *out_contig_info = xfopen((opt.contig_file() + ".info").c_str(), "w");
    fprintf(out_contig_info, "%lu %lu\n", graph.size() + bubble_hist.size(), stat["total size"] + bubble_hist.sum());
    fclose(out_contig_info);

    OutputContigs(graph, contig_file, opt.output_standalone ? standalone_file : nullptr, false, opt.min_standalone);
    timer.stop();
    xinfo("Time to output: %lf\n", timer.elapsed());
  }

  // remove local low depth & output as contigs
  if (opt.prune_level >= 1) {
    FILE *add_contig_file = xfopen(opt.addi_contig_file().c_str(), "w");
    FILE *add_contig_info = xfopen((opt.addi_contig_file() + ".info").c_str(), "w");

    timer.reset();
    timer.start();
    uint32_t num_removed = IterateLocalLowDepth(graph, opt.min_depth, opt.max_tip_len, opt.local_width,
                                                opt.low_local_ratio, opt.is_final_round);

    uint32_t n_bubbles = 0;
    if (opt.bubble_level >= 2 && opt.merge_len > 0) {
      complex_bubble_remover.set_bubble_file(nullptr);
      n_bubbles = complex_bubble_remover.PopBubbles(graph, false);
      timer.stop();
    }
    xinfo("Number of local low depth unitigs removed: %lu, complex bubbles removed: %u, time: %lf\n", num_removed,
          n_bubbles, timer.elapsed());
    CalcAndPrintStat(graph);

    if (!opt.is_final_round) {
      OutputContigs(graph, add_contig_file, nullptr, true, 0);
    } else {
      OutputContigs(graph, contig_file, opt.output_standalone ? standalone_file : nullptr, false, opt.min_standalone);
    }

    auto stat_changed = CalcAndPrintStat(graph, false, true);
    fprintf(add_contig_info, "%lu %lu\n", stat_changed["number contigs"], stat_changed["total size"]);

    fclose(add_contig_file);
    fclose(add_contig_info);
  }

  fclose(contig_file);
  fclose(standalone_file);
  fclose(bubble_file);

  return 0;
}