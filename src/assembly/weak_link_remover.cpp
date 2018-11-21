//
// Created by vout on 11/21/18.
//

#include "weak_link_remover.h"
#include "unitig_graph.h"

uint32_t DisconnectWeakLinks(UnitigGraph &graph, double local_ratio = 0.1) {
  uint32_t num_disconnected = 0;
#pragma omp parallel for reduction(+: num_disconnected)
  for (UnitigGraph::size_type i = 0; i < graph.size(); ++i) {
    auto adapter = graph.MakeVertexAdapter(i);
    if (adapter.is_loop() || adapter.is_palindrome()) {
      continue;
    }
    for (int strand = 0; strand < 2; ++strand, adapter.ReverseComplement()) {
      UnitigGraph::VertexAdapter nexts[4];
      double depths[4];
      double total_depth = 0;
      int degree = graph.GetNextAdapters(adapter, nexts);
      if (degree <= 1) {
        continue;
      }
      for (int j = 0; j < degree; ++j) {
        depths[j] = nexts[j].avg_depth();
        total_depth += depths[j];
      }
      for (int j = 0; j < degree; ++j) {
        if (depths[j] <= local_ratio * total_depth) {
          num_disconnected += !nexts[j].is_to_disconnect();
          nexts[j].MarkToDisconnect();
        }
      }
    }
  }
  graph.Refresh(false);
  return num_disconnected;
}