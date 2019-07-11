//
// Created by vout on 11/21/18.
//

#ifndef MEGAHIT_BUBBLE_REMOVER_H
#define MEGAHIT_BUBBLE_REMOVER_H

#include <cstdint>
#include "contig_output.h"
#include "unitig_graph.h"
#include "utils/histgram.h"

class BaseBubbleRemover {
 public:
  using checker_type = std::function<bool(const UnitigGraph::VertexAdapter &,
                                          const UnitigGraph::VertexAdapter &)>;

 public:
  BaseBubbleRemover &SetWriter(ContigWriter *bubble_file) {
    bubble_file_ = bubble_file;
    return *this;
  }
  BaseBubbleRemover &SetCarefulThreshold(double threshold) {
    careful_threshold_ = threshold;
    return *this;
  }

 private:
  ContigWriter *bubble_file_{};
  double careful_threshold_{1 + 1e-3};

 protected:
  size_t PopBubbles(UnitigGraph &graph, bool permanent_rm, uint32_t max_len,
                    const checker_type &checker);
  int SearchAndPopBubble(UnitigGraph &graph,
                         UnitigGraph::VertexAdapter &adapter, uint32_t max_len,
                         const checker_type &checker);
};

class NaiveBubbleRemover : public BaseBubbleRemover {
 public:
  size_t PopBubbles(UnitigGraph &graph, bool permanent_rm) {
    return BaseBubbleRemover::PopBubbles(graph, permanent_rm, graph.k() + 2,
                                         Check);
  }

 private:
  static constexpr bool Check(const UnitigGraph::VertexAdapter &a,
                              const UnitigGraph::VertexAdapter &b) {
    return true;
  }
};

class ComplexBubbleRemover : public BaseBubbleRemover {
 private:
  int merge_level_{20};
  double similarity_{0.95};

 public:
  ComplexBubbleRemover &SetMergeLevel(int merge_level) {
    merge_level_ = merge_level;
    return *this;
  }
  ComplexBubbleRemover &SetMergeSimilarity(double similarity) {
    similarity_ = similarity;
    return *this;
  }
  size_t PopBubbles(UnitigGraph &graph, bool permanent_rm);
};

#endif  // MEGAHIT_BUBBLE_REMOVER_H
