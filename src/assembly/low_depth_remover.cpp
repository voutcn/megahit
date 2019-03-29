//
// Created by vout on 11/21/18.
//

#include "low_depth_remover.h"
#include "unitig_graph.h"

namespace {

double LocalDepth(UnitigGraph &graph, UnitigGraph::VertexAdapter &adapter, uint32_t local_width) {
  double total_depth = 0;
  uint64_t num_added_edges = 0;

  for (int strand = 0; strand < 2; ++strand, adapter.ReverseComplement()) {
    UnitigGraph::VertexAdapter outs[4];
    int degree = graph.GetNextAdapters(adapter, outs);

    for (int i = 0; i < degree; ++i) {
      if (outs[i].length() <= local_width) {
        num_added_edges += outs[i].length();
        total_depth += outs[i].total_depth();
      } else {
        num_added_edges += local_width;
        total_depth += outs[i].avg_depth() * local_width;
      }
    }
  }

  if (num_added_edges == 0) {
    return 0;
  } else {
    return total_depth / num_added_edges;
  }
}

}

bool RemoveLocalLowDepth(UnitigGraph &graph, double min_depth, uint32_t max_len,
                         uint32_t local_width, double local_ratio,
                         bool permanent_rm, uint32_t *num_removed) {
  bool is_changed = false;
  bool need_refresh = false;
  uint32_t removed = 0;

#pragma omp parallel for reduction(+: removed)
  for (UnitigGraph::size_type i = 0; i < graph.size(); ++i) {
    auto adapter = graph.MakeVertexAdapter(i);
    if (adapter.forsure_standalone() || adapter.length() > max_len) {
      continue;
    }
    int indegree = graph.InDegree(adapter);
    int outdegree = graph.OutDegree(adapter);
    if (indegree + outdegree == 0) {
      continue;
    }

    if ((indegree <= 1 && outdegree <= 1) || indegree == 0 || outdegree == 0) {
      double depth = adapter.avg_depth();
      if (is_changed && depth > min_depth)
        continue;
      double mean = LocalDepth(graph, adapter, local_width);
      double threshold = min_depth;

      if (min_depth < mean * local_ratio)
        is_changed = true;
      else
        threshold = mean * local_ratio;

      if (depth < threshold) {
        is_changed = true;
        need_refresh = true;
        adapter.set_to_delete();
        ++removed;
      }
    }
  }

  if (need_refresh) {
    bool set_changed = !permanent_rm;
    graph.Refresh(set_changed);
  }
  *num_removed = removed;
  return is_changed;
}

uint32_t IterateLocalLowDepth(UnitigGraph &graph, double min_depth, uint32_t min_len,
                              uint32_t local_width, double local_ratio, bool permanent_rm) {
  uint32_t num_removed = 0;
  while (min_depth < kMaxMul) {
    if (!RemoveLocalLowDepth(graph, min_depth, min_len, local_width,
                             local_ratio, permanent_rm, &num_removed)) {
      break;
    }
    min_depth *= 1.1;
  }
  return num_removed;
}

uint32_t RemoveLowDepth(UnitigGraph &graph, double min_depth) {
  uint32_t num_removed = 0;
#pragma omp parallel for reduction(+: num_removed)
  for (UnitigGraph::size_type i = 0; i < graph.size(); ++i) {
    auto adapter = graph.MakeVertexAdapter(i);
    if (adapter.avg_depth() < min_depth) {
      adapter.set_to_delete();
      ++num_removed;
    }
  }
  graph.Refresh(false);
  return num_removed;
}