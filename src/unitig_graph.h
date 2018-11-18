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


struct UnitigVertex {
  UnitigVertex(int64_t start_node, int64_t end_node,
               int64_t rev_start_node, int64_t rev_end_node, int64_t depth, uint32_t length) :
      start_node(start_node), end_node(end_node), rev_start_node(rev_start_node),
      rev_end_node(rev_end_node), depth(depth), length(length) {
    is_deleted = false;
    is_changed = false;
    is_marked = false;
    is_dead = false;
    is_loop = false;
    is_palindrome = false;
    Normalize();
  }

  uint64_t start_node, end_node;
  uint64_t rev_start_node, rev_end_node;
  int64_t depth: 60; // if is_loop, depth is equal to average depth
  bool is_deleted: 1;
  bool is_changed: 1;
  bool is_marked: 1;
  bool is_dead: 1;
  bool is_loop: 1;
  uint32_t length: 30;
  bool is_palindrome: 1;

  UnitigVertex ReverseComplement() const {
    UnitigVertex ret = *this;
    std::swap(ret.start_node, ret.rev_start_node);
    std::swap(ret.end_node, ret.rev_end_node);
    return ret;
  }

  void Normalize() {
    if (start_node < rev_start_node) {
      std::swap(start_node, rev_start_node);
      std::swap(end_node, rev_end_node);
    }
  }

  int64_t Representation() const {
    return std::max(start_node, rev_start_node);
  }

  bool operator<(const UnitigVertex &rhs) const {
    if (is_deleted != rhs.is_deleted) return is_deleted < rhs.is_deleted;
    return Representation() < rhs.Representation();
  }

  double avg_depth() const {
    return depth * 1. / length;
  }
};

class UnitigGraph {
 public:
  typedef uint32_t size_type;
 public:
  struct VertexAdapter {
    VertexAdapter() = default;
    VertexAdapter(UnitigVertex *vertex, size_type id, bool is_rc)
    : vertex(vertex), id(id), is_rc(is_rc) {}
    UnitigVertex *vertex{};
    size_type id;
    bool is_rc{false};
    bool valid() {
      return vertex != nullptr;
    }
    uint64_t start() const {
      return is_rc ? vertex->rev_start_node : vertex->start_node;
    }
    uint64_t end() const {
      return is_rc ? vertex->rev_end_node : vertex->end_node;
    }
    uint64_t rc_start() const {
      return is_rc ? vertex->start_node : vertex->rev_start_node;
    }
    uint64_t rc_end() const {
      return is_rc ? vertex->end_node : vertex->rev_end_node;
    }
  };
 public:
  UnitigGraph(SuccinctDBG *sdbg) : sdbg_(sdbg) {}
  ~UnitigGraph() = default;

  void BuildFromSDBG();
  uint32_t size() {
    return vertices_.size();
  }
  int64_t RemoveLowDepth(double min_depth);
  bool RemoveLocalLowDepth(double min_depth,
                           int min_len,
                           int local_width,
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
  uint32_t RemoveTips(int max_tip_len);
  // output
  void OutputContigs(FILE *contig_file,
                     FILE *final_file,
                     Histgram<int64_t> &hist,
                     bool change_only,
                     int min_final_standalone);

 private:
  uint32_t MergeSimpleBubbles(
      bool permanent_rm,
      bool careful,
      FILE *bubble_file,
      Histgram<int64_t> &hist,
      int max_bubble_len,
      double careful_threshold,
      const std::function<bool(const UnitigVertex &, const UnitigVertex &)> &check_mergable);
  double LocalDepth(size_type id, int local_width);
  void Refresh(bool set_changed = true);

 private:

  int GetNextAdapters(const VertexAdapter &v, VertexAdapter *out = nullptr) {
    uint64_t out_sdbg_id[4];
    int degree = sdbg_->OutgoingEdges(v.end(), out_sdbg_id);
    if (out) {
      for (int i = 0; i < degree; ++i) {
        auto vid = start_node_map_.at(out_sdbg_id[i]);
        auto &vertex = vertices_[vid];
        out[i] = VertexAdapter(&vertex, vid, vertex.start_node != out_sdbg_id[i]);
      }
    }
    return degree;
  }

  int GetPrevAdapters(const VertexAdapter &v, VertexAdapter *out = nullptr) {
    uint64_t out_sdbg_id[4];
    int degree = sdbg_->OutgoingEdges(v.rc_end(), out_sdbg_id);
    if (out) {
      for (int i = 0; i < degree; ++i) {
        size_type vid = start_node_map_.at(out_sdbg_id[i]);
        auto &vertex = vertices_[vid];
        out[i] = VertexAdapter(&vertex, vid, vertex.rev_start_node != out_sdbg_id[i]);
      }
    }
    return degree;
  }

  VertexAdapter UniqueNextAdapter(const VertexAdapter &v) {
    int64_t next_sdbg_id = sdbg_->UniqueNextEdge(v.end());
    if (next_sdbg_id == -1) {
      return VertexAdapter();
    } else {
      size_type vid = start_node_map_.at(next_sdbg_id);
      auto &vertex = vertices_[vid];
      return VertexAdapter(&vertex, vid, vertex.start_node != next_sdbg_id);
    }
  }

  VertexAdapter UniquePrevAdapter(const VertexAdapter &v) {
    int64_t next_sdbg_id = sdbg_->UniqueNextEdge(v.rc_end());
    if (next_sdbg_id == -1) {
      return VertexAdapter();
    } else {
      size_type vid = start_node_map_.at(next_sdbg_id);
      auto &vertex = vertices_[vid];
      return VertexAdapter(&vertex, vid, vertex.rev_start_node != next_sdbg_id);
    }
  }

  VertexAdapter NextSimplePathAdapter(const VertexAdapter &v) {
    int64_t next_sdbg_id = sdbg_->NextSimplePathEdge(v.end());
    if (next_sdbg_id != -1) {
      size_type vid = start_node_map_.at(next_sdbg_id);
      auto &vertex = vertices_[vid];
      return VertexAdapter(&vertex, vid, vertex.start_node != next_sdbg_id);
    } else {
      return VertexAdapter();
    }
  }

  VertexAdapter PrevSimplePathAdapter(const VertexAdapter &v) {
    int64_t next_sdbg_id = sdbg_->NextSimplePathEdge(v.rc_end());
    if (next_sdbg_id != -1) {
      size_type vid = start_node_map_.at(next_sdbg_id);
      auto &vertex = vertices_[vid];
      return VertexAdapter(&vertex, vid, vertex.rev_start_node != next_sdbg_id);
    } else {
      return VertexAdapter();
    }
  }

 private:
  static const size_t kMaxNumVertices = std::numeric_limits<size_type>::max();
  SuccinctDBG *sdbg_;
  spp::sparse_hash_map<int64_t, size_type> start_node_map_;
  std::deque<UnitigVertex> vertices_;
  kmlib::AtomicBitVector<> locks_;
};

#endif // UNITIG_GRAPH_H_