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

#include "unitig_graph.h"

#include <omp.h>
#include <queue>
#include <set>
#include <vector>
#include <algorithm>
#include <string>
#include <utility>
#include <map>

#include "utils.h"
#include "definitions.h"
#include "sdbg/sdbg.h"
#include "assembly_algorithms.h"
#include "kmlib/kmbitvector.h"

// -- helper functions --
static inline char Complement(char c) {
  if (c >= 0 && c < 4) {
    return 3 - c;
  }

  switch (c) {
    case 'A': {
      return 'T';
    }

    case 'C': {
      return 'G';
    }

    case 'G': {
      return 'C';
    }

    case 'T': {
      return 'A';
    }

    default: {
      assert(false);
    }
  }
}

static inline void ReverseComplement(std::string &s) {
  int i, j;

  for (i = 0, j = s.length() - 1; i < j; ++i, --j) {
    std::swap(s[i], s[j]);
    s[i] = Complement(s[i]);
    s[j] = Complement(s[j]);
  }

  if (i == j) {
    s[i] = Complement(s[i]);
  }
}

std::string VertexToDNAString(const SuccinctDBG *sdbg_, const UnitigVertex &v) {
  std::string label;
  int64_t cur_edge = v.end_node;

  for (int i = 1; i < v.length; ++i) {
    int8_t cur_char = sdbg_->GetW(cur_edge);
    assert(1 <= cur_char && cur_char <= 8);
    label.append(1, "ACGT"[cur_char > 4 ? (cur_char - 5) : (cur_char - 1)]);

    cur_edge = sdbg_->PrevSimplePathEdge(cur_edge);
    // assert(cur_edge != -1);
    if (cur_edge == -1) {
      xfatal("%lld, %lld, %lld, %lld, (%lld, %lld), %d, %d\n",
             v.start_node,
             v.end_node,
             v.rev_start_node,
             v.rev_end_node,
             sdbg_->EdgeReverseComplement(v.end_node),
             sdbg_->EdgeReverseComplement(v.start_node),
             v.length,
             i);
    }
  }

  int8_t cur_char = sdbg_->GetW(cur_edge);
  label.append(1, "ACGT"[cur_char > 4 ? (cur_char - 5) : (cur_char - 1)]);

  if (cur_edge != v.start_node) {
    xfatal("fwd: %lld, %lld, rev: %lld, %lld, (%lld, %lld) length: %d\n",
           v.start_node,
           v.end_node,
           v.rev_start_node,
           v.rev_end_node,
           sdbg_->EdgeReverseComplement(v.end_node),
           sdbg_->EdgeReverseComplement(v.start_node),
           v.length);
  }

  uint8_t seq[kMaxK];
  sdbg_->Label(v.start_node, seq);

  for (int i = sdbg_->k() - 1; i >= 0; --i) {
    assert(seq[i] >= 1 && seq[i] <= 4);
    label.append(1, "ACGT"[seq[i] - 1]);
  }

  std::reverse(label.begin(), label.end());
  return label;
}

void FoldPalindrome(std::string &s, int kmer_k, bool is_loop) {
  if (is_loop) {
    for (unsigned i = 1; i + kmer_k <= s.length(); ++i) {
      std::string rc = s.substr(i, kmer_k);
      ReverseComplement(rc);

      if (rc == s.substr(i - 1, kmer_k)) {
        assert(i <= s.length() / 2);
        s = s.substr(i, s.length() / 2);
        break;
      }
    }
  } else {
    int num_edges = s.length() - kmer_k;
    assert(num_edges % 2 == 1);
    s.resize((num_edges - 1) / 2 + kmer_k + 1);
  }
}

void WriteContig(const std::string &label,
                 int k_size,
                 long long id,
                 int flag,
                 double multiplicity,
                 FILE *file) {
  std::string rev_label(label);
  ReverseComplement(rev_label);
  fprintf(file, ">k%d_%lld flag=%d multi=%.4lf len=%d\n%s\n",
          k_size,
          id,
          flag,
          multiplicity,
          (int) label.length(),
          label < rev_label ? label.c_str() : rev_label.c_str());
}

double GetSimilarity(const std::string &a, const std::string &b, double min_similar) {
  int n = a.length();
  int m = b.length();
  int max_indel = std::max(n, m) * (1 - min_similar);

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
    std::fill(dp[i & 1].begin(), dp[i & 1].end(), 99999999);

    if (i - max_indel <= 0) {
      dp[i & 1][IDX(0, i)] = i;
    }

    for (int j = std::max(i - max_indel, 1); j <= m && j <= i + max_indel; ++j) {
      // assert(IDX(j,i) >= 0 && IDX(j,i) < max_indel * 2 + 1);
      dp[i & 1][IDX(j, i)] =
          std::min(dp[i & 1][IDX(j, i)], dp[(i ^ 1) & 1][IDX(j - 1, i - 1)] + (a[i - 1] != b[j - 1]));

      if (j > i - max_indel) {
        // assert(IDX(j-1,i) >= 0 && IDX(j-1,i) < max_indel * 2 + 1);
        dp[i & 1][IDX(j, i)] = std::min(dp[i & 1][IDX(j, i)], dp[i & 1][IDX(j - 1, i)] + 1);
      }

      if (j < i + max_indel) {
        // assert(IDX(j,i-1) >= 0 && IDX(j,i-1) < max_indel * 2 + 1);
        dp[i & 1][IDX(j, i)] = std::min(dp[i & 1][IDX(j, i)], dp[(i ^ 1) & 1][IDX(j, i - 1)] + 1);
      }
    }
  }

  return 1 - dp[n & 1][IDX(m, n)] * 1.0 / std::max(n, m);
#undef IDX
}

// -- end of helper functions --

void UnitigGraph::BuildFromSDBG() {
  start_node_map_.clear();
  vertices_.clear();

  omp_lock_t path_lock;
  omp_init_lock(&path_lock);
  AtomicBitVector marked(sdbg_->size());
  size_t count_palindrome = 0;

  // assemble simple paths
#pragma omp parallel for reduction(+: count_palindrome)
  for (size_t edge_idx = 0; edge_idx < sdbg_->size(); ++edge_idx) {
    if (sdbg_->IsValidEdge(edge_idx) && sdbg_->NextSimplePathEdge(edge_idx) == -1 && marked.try_lock(edge_idx)) {
      bool will_be_added = true;
      int64_t cur_edge = edge_idx, prev_edge;
      int64_t depth = sdbg_->EdgeMultiplicity(edge_idx);
      uint32_t length = 1;

      while ((prev_edge = sdbg_->PrevSimplePathEdge(cur_edge)) != -1) {
        cur_edge = prev_edge;

        if (!marked.try_lock(cur_edge)) {
          will_be_added = false;
          break;
        }

        depth += sdbg_->EdgeMultiplicity(cur_edge);
        ++length;
      }

      if (!will_be_added) {
        continue;
      }

      int64_t rc_start = sdbg_->EdgeReverseComplement(edge_idx);
      int64_t rc_end = -1;
      assert(rc_start != -1);

      if (!marked.try_lock(rc_start)) {
        rc_end = sdbg_->EdgeReverseComplement(cur_edge);

        if (std::max(static_cast<int64_t>(edge_idx), cur_edge) < std::max(rc_start, rc_end)) {
          will_be_added = false;
        }
      } else {
        // lock through the rc path
        int64_t rc_cur_edge = rc_start;
        rc_end = rc_cur_edge;
        bool extend_full = true;

        while ((rc_cur_edge = sdbg_->NextSimplePathEdge(rc_cur_edge)) != -1) {
          rc_end = rc_cur_edge;

          if (!marked.try_lock(rc_cur_edge)) {
            extend_full = false;
            break;
          }
        }

        if (!extend_full) {
          rc_end = sdbg_->EdgeReverseComplement(cur_edge);
        }
      }

      if (!will_be_added) {
        continue;
      }

      omp_set_lock(&path_lock);
      vertices_.emplace_back(cur_edge, edge_idx, rc_start, rc_end, depth, length);
      vertices_.back().is_palindrome = cur_edge == rc_start;
      if (cur_edge == rc_start) {
        assert(edge_idx == rc_end);
      }
      count_palindrome += vertices_.back().is_palindrome;
      omp_unset_lock(&path_lock);
    }
  }
  xlog("Graph size without loops: %lu, palindrome: %lu\n", vertices_.size(), count_palindrome);

  // assemble looped paths
#pragma omp parallel for
  for (size_t edge_idx = 0; edge_idx < sdbg_->size(); ++edge_idx) {
    if (!marked.at(edge_idx) && sdbg_->IsValidEdge(edge_idx)) {
      omp_set_lock(&path_lock);

      if (!marked.at(edge_idx)) {
        int64_t cur_edge = edge_idx;
        int64_t depth = sdbg_->EdgeMultiplicity(edge_idx);
        uint32_t length = 0;

        bool rc_marked =
            marked.at(sdbg_->EdgeReverseComplement(edge_idx)); // whether it is marked before entering the loop

        while (!marked.at(cur_edge)) {
          marked.set(cur_edge);
          depth += sdbg_->EdgeMultiplicity(cur_edge);
          ++length;

          cur_edge = sdbg_->PrevSimplePathEdge(cur_edge);
          assert(cur_edge != -1);
        }

        assert(cur_edge == edge_idx);

        if (!rc_marked) {
          int64_t start = sdbg_->NextSimplePathEdge(edge_idx);
          int64_t end = edge_idx;
          vertices_.emplace_back(start,
                                 end,
                                 sdbg_->EdgeReverseComplement(end),
                                 sdbg_->EdgeReverseComplement(start),
                                 depth,
                                 length);
          vertices_.back().is_loop = true;
          vertices_.back().is_deleted =
              true; // loop path will not process to further steps, but still should be there for output

          if (marked.at(sdbg_->EdgeReverseComplement(edge_idx))) {
            // this loop is palindrome
            vertices_.back().is_palindrome = true;
          }
        }
      }

      omp_unset_lock(&path_lock);
    }
  }

  if (vertices_.size() >= kMaxNumVertices) {
    xfatal("[ERROR] Too many vertices in the unitig graph (%llu >= %llu)\n",
           (unsigned long long) vertices_.size(), (unsigned long long) kMaxNumVertices);
  }

  // free memory for hash table construction
  sdbg_->FreeMultiplicity();
  {
    AtomicBitVector empty_abv;
    marked.swap(empty_abv);
  }

  start_node_map_.reserve(vertices_.size() * 2);
  for (size_type i = 0; i < vertices_.size(); ++i) {
    if (!vertices_[i].is_deleted) {
      start_node_map_[vertices_[i].start_node] = i;
      start_node_map_[vertices_[i].rev_start_node] = i;
    }
  }
  locks_.reset(vertices_.size());
  omp_destroy_lock(&path_lock);
}

uint32_t UnitigGraph::RemoveTips(int max_tip_len) {
  uint32_t num_removed = 0;
  for (int thre = 2; thre < max_tip_len; thre = std::min(thre * 2, max_tip_len)) {
#pragma omp parallel for schedule(static) reduction(+: num_removed)
    for (size_type i = 0; i < vertices_.size(); ++i) {
      if (vertices_[i].is_deleted || vertices_[i].length >= thre) {
        continue;
      }
      VertexAdapter adapter(&vertices_[i], i, false);
      VertexAdapter nexts[4], prevs[4];
      int outd = GetNextAdapters(adapter);
      int ind = GetNextAdapters(adapter);

      if (ind + outd == 0) {
        vertices_[i].is_dead = true;
      } else if (outd == 1 && ind == 0) {
        if (!UniquePrevAdapter(nexts[0]).valid()) {
          vertices_[i].is_dead = true;
        }
      } else if (outd == 0 && ind == 1) {
        if (!UniqueNextAdapter(prevs[0]).valid()) {
          vertices_[i].is_dead = true;
        }
      }
      num_removed += vertices_[i].is_dead;
    }

    Refresh(false);
    if (thre >= max_tip_len) { break; }
  }
  return num_removed;
}

uint32_t UnitigGraph::MergeBubbles(bool permanent_rm, bool careful, FILE *bubble_file, Histgram<int64_t> &hist) {
  int max_bubble_len = sdbg_->k() + 2; // allow 1 indel
  return MergeSimpleBubbles(
      permanent_rm, careful, bubble_file, hist, max_bubble_len, 0.2,
      [=](const UnitigVertex &a, const UnitigVertex &b) {
        return a.length <= max_bubble_len && b.length <= max_bubble_len;
      }
  );
}

uint32_t UnitigGraph::MergeComplexBubbles(double similarity,
                                          int merge_level,
                                          bool permanent_rm,
                                          bool careful,
                                          FILE *bubble_file,
                                          Histgram<int64_t> &hist) {
  int max_bubble_len = sdbg_->k() * merge_level / similarity + 0.5;
  if (max_bubble_len * (1 - similarity) < 1) {
    return 0;
  }
  return MergeSimpleBubbles(
      permanent_rm, careful, bubble_file, hist, max_bubble_len, 0,
      [&, this](const UnitigVertex &a, const UnitigVertex &b) -> bool {
        return a.length <= max_bubble_len && b.length <= max_bubble_len &&
            (b.length + sdbg_->k() - 1) * similarity <= (a.length + sdbg_->k() - 1) &&
            (a.length + sdbg_->k() - 1) * similarity <= (b.length + sdbg_->k() - 1) &&
            GetSimilarity(VertexToDNAString(sdbg_, a), VertexToDNAString(sdbg_, b), similarity) >= similarity;
      }
  );
}

int64_t UnitigGraph::RemoveLowDepth(double min_depth) {
  int64_t num_removed = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+: num_removed)

  for (size_type i = 0; i < vertices_.size(); ++i) {
    if (!vertices_[i].is_deleted && vertices_[i].depth < min_depth) {
      vertices_[i].is_dead = true;
      ++num_removed;
    }
  }

  Refresh(false);
  return num_removed;
}

bool UnitigGraph::RemoveLocalLowDepth(double min_depth,
                                      int min_len,
                                      int local_width,
                                      double local_ratio,
                                      int64_t &num_removed,
                                      bool permanent_rm) {
  bool is_changed = false;
  bool need_refresh = false;
  int64_t num_removed_ = 0;

#pragma omp parallel for schedule(dynamic, 1) reduction(+: num_removed_)

  for (size_type i = 0; i < vertices_.size(); ++i) {
    if (vertices_[i].is_deleted || vertices_[i].length >= min_len) {
      continue;
    }
    assert(vertices_[i].length > 0);
    VertexAdapter adapter(&vertices_[i], i, false);
    int indegree = GetNextAdapters(adapter);
    int outdegree = GetPrevAdapters(adapter);
    if (indegree + outdegree == 0) {
      continue;
    }

    if ((indegree <= 1 && outdegree <= 1) || indegree == 0 || outdegree == 0) {
      double depth = (double) vertices_[i].depth / vertices_[i].length;
      if (is_changed && depth > min_depth)
        continue;
      double mean = LocalDepth(i, local_width);
      double threshold = min_depth;

      if (min_depth < mean * local_ratio)
        is_changed = true;
      else
        threshold = mean * local_ratio;

      if (depth < threshold) {
        is_changed = true;
        need_refresh = true;
        vertices_[i].is_dead = true;
        ++num_removed_;
      }
    }
  }

  if (need_refresh) {
    bool set_changed = !permanent_rm;
    Refresh(set_changed);
  }
  num_removed += num_removed_;
  return is_changed;
}

uint32_t UnitigGraph::DisconnectWeakLinks(double local_ratio = 0.1) {
  // see metaspades paper
  vector<VertexAdapter> to_remove;
  std::mutex mutex;
#pragma omp parallel for schedule(dynamic, 1)
  for (size_type i = 0; i < vertices_.size(); ++i) {
    if (vertices_[i].is_deleted) {
      continue;
    }

    for (int strand = 0; strand < 2; ++strand) {
      VertexAdapter adapter(&vertices_[i], i, strand != 0);
      VertexAdapter next[4];
      double depths[4];
      double total_depth = 0;
      int degree = GetNextAdapters(adapter, next);

      if (degree <= 1) {
        continue;
      }

      for (int j = 0; j < degree; ++j) {
        depths[j] = next[j].vertex->avg_depth();
        total_depth += depths[j];
      }

      for (int j = 0; j < degree; ++j) {
        if (depths[j] <= total_depth * local_ratio) {
          std::lock_guard<std::mutex> lk(mutex);
          to_remove.push_back(next[j]);
        }
      }
    }
  }
  std::sort(to_remove.begin(), to_remove.end(),
            [](const VertexAdapter &a, const VertexAdapter &b) {
              return a.id < b.id;
            }
  );

  size_type num_disconnected = 0;
#pragma omp parallel for reduction(+: num_disconnected)
  for (size_type i = 0; i < to_remove.size(); ++i) {
    if ((i > 0 && to_remove[i].id == to_remove[i - 1].id) || to_remove[i].vertex->is_palindrome) {
      continue; // already handled by i - 1
    }

    bool dual_disconnected = i + 1 < to_remove.size() && to_remove[i + 1].id == to_remove[i].id;
    if (to_remove[i].vertex->length <= 1 + dual_disconnected) {
      to_remove[i].vertex->is_dead = true;
      continue;
    }

    auto &adapter = to_remove[i];
    auto old_start = adapter.start();
    auto old_end = adapter.end();
    auto old_rc_start = adapter.rc_start();
    auto old_rc_end = adapter.rc_end();
    size_type id = adapter.id;
    UnitigVertex &vertex = *adapter.vertex;
    int64_t new_start, new_end, new_rc_start, new_rc_end;

    if (dual_disconnected) {
      new_start = sdbg_->NextSimplePathEdge(old_start);
      new_rc_end = sdbg_->PrevSimplePathEdge(old_rc_end);
      new_rc_start = sdbg_->NextSimplePathEdge(old_rc_start);
      new_end = sdbg_->PrevSimplePathEdge(old_end);
    } else {
      new_start = sdbg_->NextSimplePathEdge(old_start);
      new_rc_end = sdbg_->PrevSimplePathEdge(old_rc_end);
      new_rc_start = old_rc_start;
      new_end = old_end;
    }

    assert(new_start != -1);
    assert(new_end != -1);
    assert(new_rc_start != -1);
    assert(new_rc_end != -1);

    sdbg_->SetInvalidEdge(old_start);
    sdbg_->SetInvalidEdge(old_rc_end);
    if (dual_disconnected) {
      sdbg_->SetInvalidEdge(old_rc_start);
      sdbg_->SetInvalidEdge(old_end);
    }

    vertex.start_node = new_start;
    vertex.end_node = new_end;
    vertex.rev_start_node = new_rc_start;
    vertex.rev_end_node = new_rc_end;
    vertex.depth *= (vertex.length - 1. - dual_disconnected) / vertex.length;
    vertex.length -= 1 + dual_disconnected;
    num_disconnected++;

    {
      std::lock_guard<std::mutex> lk(mutex);
      start_node_map_.erase(old_start);
      start_node_map_[new_start] = id;
      if (dual_disconnected) {
        start_node_map_.erase(old_rc_start);
        start_node_map_[new_rc_start] = id;
      }
    }
  }

  Refresh(true);
}

double UnitigGraph::LocalDepth(size_type id, int local_width) {
  double total_depth = 0;
  double num_added_edges = 0;

  for (int strand = 0; strand < 2; ++strand) {
    VertexAdapter adapter(&vertices_[id], id, strand != 0);
    VertexAdapter outs[4];
    int degree = GetNextAdapters(adapter, outs);

    for (int i = 0; i < degree; ++i) {
      auto &next_vertex = *outs[i].vertex;
      if (next_vertex.length <= local_width) {
        num_added_edges += next_vertex.length;
        total_depth += next_vertex.depth;
      } else {
        num_added_edges += local_width;
        total_depth += (double) next_vertex.depth * local_width / next_vertex.length;
      }
    }
  }

  if (num_added_edges == 0) {
    return 0;
  } else {
    return total_depth / num_added_edges;
  }
}

void UnitigGraph::Refresh(bool set_changed) {
#pragma omp parallel for
  for (size_type i = 0; i < vertices_.size(); ++i) {
    if (vertices_[i].is_dead && !vertices_[i].is_deleted) {
      int64_t cur_edge = vertices_[i].end_node, prev_edge;

      while (cur_edge != vertices_[i].start_node) {
        prev_edge = sdbg_->UniquePrevEdge(cur_edge);
        sdbg_->SetInvalidEdge(cur_edge);
        cur_edge = prev_edge;
        assert(cur_edge != -1);
      }

      sdbg_->SetInvalidEdge(cur_edge);

      if (vertices_[i].rev_end_node != vertices_[i].end_node) {
        cur_edge = vertices_[i].rev_end_node;

        while (cur_edge != vertices_[i].rev_start_node) {
          prev_edge = sdbg_->UniquePrevEdge(cur_edge);
          sdbg_->SetInvalidEdge(cur_edge);
          cur_edge = prev_edge;
          assert(cur_edge != -1);
        }
        sdbg_->SetInvalidEdge(cur_edge);
      }
      vertices_[i].is_deleted = true;
    }
  }

#pragma omp parallel for
  for (size_type i = 0; i < vertices_.size(); ++i) {
    if (vertices_[i].is_deleted) {
      continue;
    }
    for (int strand = 0, done = false; !done && strand < 2; ++strand) {
      VertexAdapter adapter(&vertices_[i], i, strand != 0);
      if (PrevSimplePathAdapter(adapter).valid()) {
        continue;
      } else {
        done = true;
      }
      if (!locks_.try_lock(i)) {
        break;
      }
      std::vector<VertexAdapter> linear_path;
      for (auto cur = NextSimplePathAdapter(adapter); cur.valid(); cur = NextSimplePathAdapter(cur)) {
        linear_path.emplace_back(cur);
      }

      if (linear_path.empty()) {
        adapter.vertex->is_marked = true;
        break;
      }

      if (linear_path.back().id != adapter.id &&
          !locks_.try_lock(linear_path.back().id)) {
        if (linear_path.back().id > i) {
          locks_.unlock(i);
          break;
        } else {
          locks_.lock(linear_path.back().id);
        }
      }
      adapter.vertex->is_marked = true;
      auto length = adapter.vertex->length;
      auto depth = adapter.vertex->depth;
      for (auto &v: linear_path) {
        length += v.vertex->length;
        depth += v.vertex->depth;
        v.vertex->is_deleted = true;
      }

      auto new_start = adapter.start();
      auto new_rc_end = adapter.rc_end();
      auto new_rc_start = linear_path.back().rc_start();
      auto new_end = linear_path.back().end();

      vertices_[i].start_node = new_start;
      vertices_[i].end_node = new_end;
      vertices_[i].rev_start_node = new_rc_start;
      vertices_[i].rev_end_node = new_rc_end;
      vertices_[i].is_palindrome = new_start == new_rc_start;
      vertices_[i].is_changed |= set_changed;
      vertices_[i].depth = depth;
      vertices_[i].length = length;
      vertices_[i].is_deleted = false;
    }
  }

  // looped path
  std::mutex mutex;
#pragma omp parallel for
  for (size_type i = 0; i < vertices_.size(); ++i) {
    if (!vertices_[i].is_deleted && !vertices_[i].is_marked) {
      std::lock_guard<std::mutex> lk(mutex);
      if (vertices_[i].is_deleted) {
        continue;
      }
      VertexAdapter adapter(&vertices_[i], i, false);
      uint32_t length = adapter.vertex->length;
      int64_t depth = adapter.vertex->depth;

      adapter.vertex->is_changed |= set_changed;
      adapter.vertex->is_loop = true;
      adapter.vertex->is_deleted = true;
      bool is_palindrome = false;

      for (auto next = NextSimplePathAdapter(adapter); next.start() != adapter.start();
           next = NextSimplePathAdapter(next)) {
        assert(next.valid());
        if (next.vertex->is_deleted) {
          is_palindrome = true;
        }
        length += next.vertex->length;
        depth += next.vertex->depth;
        next.vertex->is_deleted = true;
      }

      vertices_[i].depth = depth;
      vertices_[i].length = length;
      vertices_[i].is_palindrome = is_palindrome;
      vertices_[i].end_node = sdbg_->PrevSimplePathEdge(vertices_[i].start_node);
      vertices_[i].rev_start_node = sdbg_->EdgeReverseComplement(vertices_[i].end_node);
      vertices_[i].rev_end_node = sdbg_->EdgeReverseComplement(vertices_[i].start_node);
    }
  }

#pragma omp parallel for
  for (size_type i = 0; i < vertices_.size(); ++i) {
    if (!vertices_[i].is_deleted) {
      vertices_[i].is_marked = false;
      vertices_[i].Normalize();
      start_node_map_.at(vertices_[i].start_node) = i;
      start_node_map_.at(vertices_[i].rev_start_node) = i;
    } else {
      assert(!vertices_[i].is_marked);
    }
  }
  locks_.reset();
}

void UnitigGraph::OutputContigs(FILE *contig_file, FILE *final_file, Histgram<int64_t> &histo,
                                bool change_only, int min_final_standalone) {
  histo.clear();
  assert(!(change_only && final_file != NULL)); // if output changed contigs, must not output final contigs

#pragma omp parallel for
  for (size_type i = 0; i < vertices_.size(); ++i) {
    if (vertices_[i].is_deleted && !vertices_[i].is_loop) {
      continue;
    }

    double multi = std::min((double) kMaxMul, (double) vertices_[i].depth / vertices_[i].length);

    if (change_only) {
      multi = 1;
    }

    std::string label = VertexToDNAString(sdbg_, vertices_[i]);

    if (vertices_[i].is_palindrome && vertices_[i].is_loop) {
      FoldPalindrome(label, sdbg_->k(), vertices_[i].is_loop);
    }

    histo.insert(label.length());

    if (change_only && !vertices_[i].is_changed) {
      continue;
    }

    if (vertices_[i].is_loop) {
      int flag = contig_flag::kLoop | contig_flag::kIsolated;
      FILE *out_file = contig_file;

      if (vertices_[i].is_palindrome) {
        flag = contig_flag::kIsolated;
      }

      if (final_file != NULL) {
        if (label.length() < (unsigned) min_final_standalone) {
          continue;
        } else {
          out_file = final_file;
        }
      }

      WriteContig(label, sdbg_->k(), i, flag, multi, out_file);

    } else {
      FILE *out_file = contig_file;
      int flag = 0;

      int indegree = sdbg_->EdgeIndegree(vertices_[i].start_node);
      int outdegree = sdbg_->EdgeOutdegree(vertices_[i].end_node);

      if (indegree == 0 && outdegree == 0) {
        vertices_[i].is_deleted = true;

        if (vertices_[i].is_palindrome) {
          FoldPalindrome(label, sdbg_->k(), vertices_[i].is_loop);
        }

        flag = contig_flag::kIsolated;

        if (final_file != NULL) {
          if (label.length() < (unsigned) min_final_standalone) {
            continue;
          } else {
            out_file = final_file;
          }
        }
      }

      WriteContig(label, sdbg_->k(), i, flag, multi, out_file);
    }
  }
}

uint32_t UnitigGraph::MergeSimpleBubbles(
    bool permanent_rm,
    bool careful,
    FILE *bubble_file,
    Histgram<int64_t> &hist,
    int max_bubble_len,
    double careful_threshold,
    const std::function<bool(const UnitigVertex &, const UnitigVertex &)> &check_mergable) {
  uint32_t num_removed = 0;
  std::vector<VertexAdapter> branches;

#pragma omp parallel for private(branches) reduction(+: num_removed)
  for (size_type i = 0; i < vertices_.size(); ++i) {
    if (vertices_[i].is_deleted) {
      continue;
    }
    bool handled_by_reverse = false;
    for (int strand = 0; !handled_by_reverse && strand < 2; ++strand) {
      VertexAdapter left = {&vertices_[i], i, strand != 0};
      VertexAdapter right;
      VertexAdapter middle[4];
      int degree = GetNextAdapters(left, middle);
      if (degree <= 1) {
        continue;
      }

      branches.clear();
      bool bubble_valid = true;

      for (int j = 0; j < degree; ++j) {
        auto rhs = UniqueNextAdapter(middle[j]);
        if (!rhs.valid()) {
          bubble_valid = false;
          break;
        }
        auto lhs = UniquePrevAdapter(middle[j]);
        if (!lhs.valid()) {
          bubble_valid = false;
          break;
        }

        if (j == 0) {
          right = rhs;
          if (rhs.vertex->Representation() < vertices_[i].Representation()) {
            handled_by_reverse = true;
            bubble_valid = false;
            break;
          }
          if (GetPrevAdapters(rhs, nullptr) != degree) {
            bubble_valid = false;
            break;
          }
        } else {
          if (right.start() != rhs.start()) {
            bubble_valid = false;
            break;
          }
        }
        branches.emplace_back(middle[j]);
      }

      if (!bubble_valid) {
        continue;
      }

      std::sort(branches.begin(), branches.end(),
                [](const VertexAdapter &a, const VertexAdapter &b) {
                  return a.vertex->avg_depth() > b.vertex->avg_depth();
                });
      bool careful_merged = false;

      auto &v0 = *(branches[0].vertex);
      for (int j = 1; j < degree; ++j) {
        auto &vj = *(branches[j].vertex);
        if (!check_mergable(v0, vj)) {
          bubble_valid = false;
          break;
        }
      }
      if (!bubble_valid) {
        continue;
      }
      for (int j = 1; j < degree; ++j) {
        auto &vj = *(branches[j].vertex);
        num_removed += 1;
        vj.is_dead = true;

        if (careful && vj.avg_depth() >= v0.avg_depth() * careful_threshold) {
          careful_merged = true;
          string label = VertexToDNAString(sdbg_, vj);
          WriteContig(label, sdbg_->k(), 0, 0, v0.avg_depth(), bubble_file);
          hist.insert(label.length());
        }
      }

      if (careful_merged) {
        UnitigVertex &leftVertex = *left.vertex;
        UnitigVertex &rightVertex = *right.vertex;
        string left_label = VertexToDNAString(sdbg_, leftVertex);
        string right_label = VertexToDNAString(sdbg_, rightVertex);
        WriteContig(left_label,
                    sdbg_->k(),
                    0,
                    0,
                    leftVertex.depth * 1.0 / leftVertex.length,
                    bubble_file);
        WriteContig(right_label,
                    sdbg_->k(),
                    0,
                    0,
                    rightVertex.depth * 1.0 / rightVertex.length,
                    bubble_file);
        hist.insert(left_label.length());
        hist.insert(right_label.length());
      }

      if (left.vertex == right.vertex || left.vertex->is_palindrome) {
        break;
      }
    }
  }
  Refresh(!permanent_rm);
  return num_removed;
}
