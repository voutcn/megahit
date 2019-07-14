//
// Created by vout on 11/21/18.
//

#include "bubble_remover.h"
#include <utils/utils.h>

namespace {  // helper function

double GetSimilarity(const std::string &a, const std::string &b,
                     double min_similarity) {
  int n = a.length();
  int m = b.length();
  int max_indel = std::max(n, m) * (1 - min_similarity);

  if (abs(n - m) > max_indel) {
    return 0;
  }
  if (max_indel < 1) {
    return 0;
  }
  std::vector<int> dp[2];
  for (int i = 0; i < 2; ++i) {
    dp[i].resize(max_indel * 2 + 1, 0);
  }

#define IDX(j, i) ((j) - (i) + max_indel)

  for (int j = 0; j <= max_indel; ++j) {
    dp[0][IDX(j, 0)] = j;
  }
  for (int i = 1; i <= n; ++i) {
    std::fill(dp[i & 1].begin(), dp[i & 1].end(), 0x3f3f3f3f);
    if (i - max_indel <= 0) {
      dp[i & 1][IDX(0, i)] = i;
    }
    for (int j = std::max(i - max_indel, 1); j <= m && j <= i + max_indel;
         ++j) {
      dp[i & 1][IDX(j, i)] =
          std::min(dp[i & 1][IDX(j, i)],
                   dp[(i ^ 1) & 1][IDX(j - 1, i - 1)] + (a[i - 1] != b[j - 1]));
      if (j > i - max_indel) {
        dp[i & 1][IDX(j, i)] =
            std::min(dp[i & 1][IDX(j, i)], dp[i & 1][IDX(j - 1, i)] + 1);
      }
      if (j < i + max_indel) {
        dp[i & 1][IDX(j, i)] =
            std::min(dp[i & 1][IDX(j, i)], dp[(i ^ 1) & 1][IDX(j, i - 1)] + 1);
      }
    }
  }
  return 1 - dp[n & 1][IDX(m, n)] * 1.0 / std::max(n, m);
#undef IDX
}

}  // namespace

int BaseBubbleRemover::SearchAndPopBubble(UnitigGraph &graph,
                                          UnitigGraph::VertexAdapter &adapter,
                                          uint32_t max_len,
                                          const checker_type &checker) {
  UnitigGraph::VertexAdapter right;
  UnitigGraph::VertexAdapter middle[4];
  UnitigGraph::VertexAdapter possible_right[4];

  int degree = graph.GetNextAdapters(adapter, middle);
  if (degree <= 1) {
    return 0;
  }

  for (int j = 0; j < degree; ++j) {
    if (middle[j].GetLength() > max_len) {
      return 0;
    }
  }

  for (int j = 0; j < degree; ++j) {
    if (graph.InDegree(middle[j]) != 1 ||
        graph.GetNextAdapters(middle[j], possible_right) != 1) {
      return 0;
    }
    if (j == 0) {
      right = possible_right[0];
      if (right.canonical_id() < adapter.canonical_id() ||
          graph.InDegree(right) != degree) {
        return 0;
      }
    } else {
      if (right.b() != possible_right[0].b()) {
        return 0;
      }
    }
  }

  std::sort(middle, middle + degree,
            [](const UnitigGraph::VertexAdapter &a,
               const UnitigGraph::VertexAdapter &b) {
              if (a.GetAvgDepth() != b.GetAvgDepth())
                return a.GetAvgDepth() > b.GetAvgDepth();
              return a.canonical_id() < b.canonical_id();
            });

  for (int j = 1; j < degree; ++j) {
    if (!checker(middle[0], middle[j])) {
      return 0;
    }
  }

  bool careful_merged = false;
  int num_removed = 0;
  for (int j = 1; j < degree; ++j) {
    bool success = middle[j].SetToDelete();
    assert(success || adapter.canonical_id() == right.canonical_id() ||
           adapter.IsPalindrome());
    num_removed += success;
    if (bubble_file_ && middle[j].GetAvgDepth() >=
                            middle[0].GetAvgDepth() * careful_threshold_) {
      std::string label = graph.VertexToDNAString(middle[j]);
      bubble_file_->WriteContig(label, graph.k(), 0, 0,
                                middle[j].GetAvgDepth());
      careful_merged = true;
    }
  }

  if (careful_merged) {
    std::string left_label = graph.VertexToDNAString(adapter);
    std::string right_label = graph.VertexToDNAString(right);
    bubble_file_->WriteContig(left_label, graph.k(), 0, 0,
                              adapter.GetAvgDepth());
    bubble_file_->WriteContig(right_label, graph.k(), 0, 0,
                              right.GetAvgDepth());
  }
  return num_removed;
}

size_t BaseBubbleRemover::PopBubbles(UnitigGraph &graph, bool permanent_rm,
                                     uint32_t max_len,
                                     const checker_type &checker) {
  uint32_t num_removed = 0;
#pragma omp parallel for reduction(+ : num_removed)
  for (UnitigGraph::size_type i = 0; i < graph.size(); ++i) {
    UnitigGraph::VertexAdapter adapter = graph.MakeVertexAdapter(i);
    if (adapter.IsStandalone()) {
      continue;
    }
    for (int strand = 0; strand < 2; ++strand, adapter.ReverseComplement()) {
      num_removed += SearchAndPopBubble(graph, adapter, max_len, checker);
    }
  }
  graph.Refresh(!permanent_rm);
  return num_removed;
}

size_t ComplexBubbleRemover::PopBubbles(UnitigGraph &graph, bool permanent_rm) {
  uint32_t k = graph.k();
  double sim = similarity_;
  uint32_t max_len = lround(merge_level_ * k / sim);
  if (max_len * (1 - similarity_) < 1) {
    return 0;
  }

  auto checker = [&graph, k, sim](const UnitigGraph::VertexAdapter &a,
                                  const UnitigGraph::VertexAdapter &b) -> bool {
    return (b.GetLength() + k - 1) * sim <= (a.GetLength() + k - 1) &&
           (a.GetLength() + k - 1) * sim <= (b.GetLength() + k - 1) &&
           GetSimilarity(graph.VertexToDNAString(a), graph.VertexToDNAString(b),
                         sim) >= sim;
  };
  return BaseBubbleRemover::PopBubbles(graph, permanent_rm, max_len, checker);
}