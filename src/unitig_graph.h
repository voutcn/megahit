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

#ifndef UNITIG_GRAPH_H_
#define UNITIG_GRAPH_H_

#include <deque>
#include <limits>
#include <functional>
#include <assert.h>
#include <omp.h>

#include "sparsepp/sparsepp/spp.h"
#include "histgram.h"
#include "sdbg/sdbg.h"
#include "assembly/unitig_graph.h"

class UnitigGraph : public assembly::UnitigGraph {
 public:
  typedef uint32_t size_type;
  typedef assembly::UnitigGraph::Vertex Vertex;
  typedef assembly::UnitigGraph::VertexAdapter VertexAdapter;
 public:
  UnitigGraph(SuccinctDBG *sdbg): assembly::UnitigGraph(sdbg) {}
  int64_t RemoveLowDepth(double min_depth);
  bool RemoveLocalLowDepth(double min_depth,
                           uint32_t min_len,
                           uint32_t local_width,
                           double local_ratio,
                           int64_t &num_removed,
                           bool permanent_rm = false);
  uint32_t MergeBubbles(bool permanent_rm, bool careful, FILE *bubble_file, Histgram<int64_t> &hist);
  uint32_t MergeComplexBubbles(double similarity,
                               int merge_level,
                               bool permanent_rm,
                               bool careful,
                               FILE *bubble_file,
                               Histgram<int64_t> &hist);
  uint32_t DisconnectWeakLinks(double local_ratio);
  uint32_t RemoveTips(uint32_t max_tip_len);
  // output
  void OutputContigs(FILE *contig_file,
                     FILE *final_file,
                     Histgram<int64_t> &hist,
                     bool change_only,
                     uint32_t min_final_standalone);

 private:
  uint32_t MergeSimpleBubbles(
      bool permanent_rm,
      bool careful,
      FILE *bubble_file,
      Histgram<int64_t> &hist,
      uint32_t max_bubble_len,
      double careful_threshold,
      const std::function<bool(const VertexAdapter &, const VertexAdapter &)> &check_mergable);
  double LocalDepth(size_type id, uint32_t local_width);
};

#endif // UNITIG_GRAPH_H_