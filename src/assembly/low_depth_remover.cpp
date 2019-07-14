//
// Created by vout on 11/21/18.
//

#include "low_depth_remover.h"
#include "unitig_graph.h"

namespace {

double LocalDepth(UnitigGraph &graph, UnitigGraph::VertexAdapter &adapter,
                  uint32_t local_width) {
  double total_depth = 0;
  uint64_t num_added_edges = 0;

  for (int strand = 0; strand < 2; ++strand, adapter.ReverseComplement()) {
    UnitigGraph::VertexAdapter outs[4];
    int degree = graph.GetNextAdapters(adapter, outs);

    for (int i = 0; i < degree; ++i) {
      if (outs[i].GetLength() <= local_width) {
        num_added_edges += outs[i].GetLength();
        total_depth += outs[i].GetTotalDepth();
      } else {
        num_added_edges += local_width;
        total_depth += outs[i].GetAvgDepth() * local_width;
      }
    }
  }

  if (num_added_edges == 0) {
    return 0;
  } else {
    return total_depth / num_added_edges;
  }
}

}  // namespace

bool RemoveLocalLowDepth(UnitigGraph &graph, double min_depth, uint32_t max_len,
                         uint32_t local_width, double local_ratio,
                         bool permanent_rm, uint32_t *num_removed) {
  bool need_refresh = false;
  uint32_t removed = 0;
  std::atomic_bool is_changed{false};

#pragma omp parallel for reduction(+ : removed) reduction(|| : need_refresh)
  for (UnitigGraph::size_type i = 0; i < graph.size(); ++i) {
    auto adapter = graph.MakeVertexAdapter(i);
    if (adapter.IsStandalone() || adapter.GetLength() > max_len) {
      continue;
    }
    int indegree = graph.InDegree(adapter);
    int outdegree = graph.OutDegree(adapter);
    if (indegree + outdegree == 0) {
      continue;
    }

    if ((indegree <= 1 && outdegree <= 1) || indegree == 0 || outdegree == 0) {
      double depth = adapter.GetAvgDepth();
      if (is_changed.load(std::memory_order_relaxed) && depth > min_depth)
        continue;
      double mean = LocalDepth(graph, adapter, local_width);
      double threshold = min_depth;

      if (min_depth < mean * local_ratio)
        is_changed.store(true, std::memory_order_relaxed);
      else
        threshold = mean * local_ratio;

      if (depth < threshold) {
        is_changed.store(true, std::memory_order_relaxed);
        need_refresh = true;
        bool success = adapter.SetToDelete();
        assert(success);
        removed += success;
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

uint32_t IterateLocalLowDepth(UnitigGraph &graph, double min_depth,
                              uint32_t min_len, uint32_t local_width,
                              double local_ratio, bool permanent_rm) {
  uint32_t total_removed = 0;
  while (min_depth < kMaxMul) {
    uint32_t num_removed = 0;
    if (!RemoveLocalLowDepth(graph, min_depth, min_len, local_width,
                             local_ratio, permanent_rm, &num_removed)) {
      break;
    }
    total_removed += num_removed;
    min_depth *= 1.1;
  }
  return total_removed;
}

uint32_t RemoveLowDepth(UnitigGraph &graph, double min_depth) {
  uint32_t num_removed = 0;
#pragma omp parallel for reduction(+ : num_removed)
  for (UnitigGraph::size_type i = 0; i < graph.size(); ++i) {
    auto adapter = graph.MakeVertexAdapter(i);
    if (adapter.GetAvgDepth() < min_depth) {
      bool success = adapter.SetToDelete();
      assert(success);
      num_removed += success;
    }
  }
  graph.Refresh(false);
  return num_removed;
}