//
// Created by vout on 11/21/18.
//

#ifndef MEGAHIT_CONTIG_STAT_H
#define MEGAHIT_CONTIG_STAT_H

#include <map>
#include "utils.h"
#include "histgram.h"
#include "unitig_graph.h"

using ContigStat = std::map<std::string, uint64_t>;

inline ContigStat CalcAndPrintStat(UnitigGraph &graph, bool print = true, bool changed_only = false) {
  uint32_t n_isolated = 0, n_looped = 0;
  Histgram<uint64_t> hist;
#pragma omp parallel for reduction(+: n_looped, n_isolated)
  for (UnitigGraph::size_type i = 0; i < graph.size(); ++i) {
    auto adapter = graph.MakeVertexAdapter(i);
    if (changed_only && !adapter.is_changed()) {
      continue;
    }
    hist.insert(adapter.length() + graph.k());
    n_looped += adapter.is_loop();
    n_isolated += adapter.forsure_standalone() || (graph.InDegree(adapter) == 0 && graph.OutDegree(adapter) == 0);
  }
  uint64_t total_size = hist.sum();
  ContigStat stat = {
      {"Max", hist.maximum()},
      {"Min", hist.minimum()},
      {"N50", hist.Nx(0.5 * total_size)},
      {"total size", total_size},
      {"number contigs", hist.size()},
      {"number looped", n_looped},
      {"number isolated", n_isolated}
  };

  if (print) {
    xinfo("")
    for (auto &kv: stat) {
      xinfoc("%s: %lu, ", kv.first.c_str(), kv.second);
    }
    xinfoc("\n");
  }
  return stat;
}

#endif //MEGAHIT_CONTIG_STAT_H
