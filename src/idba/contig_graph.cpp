/**
 * @file contig_graph.cpp
 * @brief
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-26
 */

#include "contig_graph.h"

#include <deque>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <sstream>

#include "contig_graph_branch_group.h"
#include "sequence.h"

using namespace std;

void ContigGraph::Initialize(const deque<Sequence> &contigs,
                             const deque<ContigInfo> &contig_infos) {
  vertices_.clear();
  vertices_.resize(contigs.size());

  for (int64_t i = 0; i < (int64_t)contigs.size(); ++i) {
    vertices_[i].clear();
    vertices_[i].set_contig(contigs[i]);
    vertices_[i].set_contig_info(contig_infos[i]);
    vertices_[i].set_id(i);
  }
  RefreshEdges();
}

void ContigGraph::Refresh() {
  RefreshVertices();
  RefreshEdges();
}

void ContigGraph::RefreshVertices() {
  uint64_t index = 0;
  for (unsigned i = 0; i < vertices_.size(); ++i) {
    if (!vertices_[i].status().IsDead()) {
      vertices_[index].swap(vertices_[i]);
      vertices_[index].set_id(index);
      ++index;
    }
  }
  vertices_.resize(index);
}

void ContigGraph::RefreshEdges() {
  BuildBeginIdbaKmerMap();

  uint64_t total_degree = 0;

  for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i) {
    for (int strand = 0; strand < 2; ++strand) {
      ContigGraphVertexAdaptor current(&vertices_[i], strand);

      for (int x = 0; x < 4; ++x) {
        if (current.out_edges()[x]) {
          IdbaKmer kmer = current.end_kmer(kmer_size_);
          kmer.ShiftAppend(x);
          if (FindVertexAdaptorByBeginIdbaKmer(kmer).is_null())
            current.out_edges().Remove(x);
        }
      }

      total_degree += current.out_edges().size();
    }

    if (vertices_[i].contig().size() == kmer_size_ &&
        vertices_[i].contig().IsPalindrome()) {
      vertices_[i].in_edges() =
          vertices_[i].out_edges() | vertices_[i].out_edges();
      vertices_[i].out_edges() = vertices_[i].in_edges();
    }
  }

  num_edges_ = total_degree / 2;
}

void ContigGraph::ClearStatus() {
  for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    vertices_[i].status().clear();
}

int64_t ContigGraph::Trim(int min_length) {
  uint64_t old_num_vertices = vertices_.size();

  for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i) {
    if (vertices_[i].contig().size() == kmer_size_ &&
        vertices_[i].contig().IsPalindrome())
      continue;

    if ((vertices_[i].in_edges().empty() || vertices_[i].out_edges().empty()) &&
        vertices_[i].contig().size() < min_length + kmer_size_ - 1 &&
        (vertices_[i].in_edges().size() + vertices_[i].out_edges().size() <=
         1)) {
      vertices_[i].status().SetDeadFlag();
    }
  }
  Refresh();
  MergeSimplePaths();

  return old_num_vertices - vertices_.size();
}

int64_t ContigGraph::RemoveDeadEnd(int min_length) {
  uint64_t num_deadend = 0;
  int l = 1;
  while (true) {
    l = min(2 * l, min_length);
    num_deadend += Trim(l);

    if (l == min_length) break;
  }
  num_deadend += Trim(min_length);
  return num_deadend;
}

int64_t ContigGraph::RemoveBubble() {
  deque<ContigGraphVertexAdaptor> candidates;

  for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i) {
    for (int strand = 0; strand < 2; ++strand) {
      ContigGraphVertexAdaptor current(&vertices_[i], strand);

      if (current.out_edges().size() > 1 &&
          current.contig_size() > kmer_size_) {
        ContigGraphBranchGroup branch_group(this, current, 4, kmer_size_ + 2);

        if (branch_group.Search()) {
          ContigGraphVertexAdaptor begin = branch_group.begin();
          ContigGraphVertexAdaptor end = branch_group.end();

          begin.ReverseComplement();
          end.ReverseComplement();
          std::swap(begin, end);
          ContigGraphBranchGroup rev_branch_group(this, begin, 4,
                                                  kmer_size_ + 2);

          if (rev_branch_group.Search() && rev_branch_group.end() == end) {
            candidates.push_back(current);
          }
        }
      }
    }
  }

  int64_t bubble = 0;
  for (unsigned i = 0; i < candidates.size(); ++i) {
    ContigGraphVertexAdaptor current = candidates[i];

    if (current.out_edges().size() > 1) {
      ContigGraphBranchGroup branch_group(this, current, 4, kmer_size_ + 2);

      if (branch_group.Search()) {
        ContigGraphVertexAdaptor begin = branch_group.begin();
        ContigGraphVertexAdaptor end = branch_group.end();

        begin.ReverseComplement();
        end.ReverseComplement();
        std::swap(begin, end);
        ContigGraphBranchGroup rev_branch_group(this, begin, 4, kmer_size_ + 2);

        if (rev_branch_group.Search() && rev_branch_group.end() == end) {
          branch_group.Merge();
          ++bubble;
        }
      }
    }
  }

  Refresh();
  MergeSimplePaths();

  return bubble;
}

double ContigGraph::IterateCoverage(int min_length, double min_cover,
                                    double max_cover, double factor) {
  min_cover = min(min_cover, max_cover);
  while (true) {
    RemoveLowCoverage(min_cover, min_length);
    min_cover *= factor;
    if (min_cover >= max_cover) break;
  }
  return min_cover;
}

bool ContigGraph::RemoveLowCoverage(double min_cover, int min_length) {
  bool is_changed = false;

  for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i) {
    ContigGraphVertexAdaptor current(&vertices_[i]);

    if (current.contig_size() < min_length + kmer_size_ - 1 &&
        ((current.in_edges().size() <= 1 && current.out_edges().size() <= 1) ||
         current.in_edges().size() == 0 || current.out_edges().size() == 0)) {
      if (current.coverage() < min_cover) {
        is_changed = true;
        current.status().SetDeadFlag();
      }
    }
  }

  Refresh();
  // Trim(min_length);
  MergeSimplePaths();

  return is_changed;
}

void ContigGraph::MergeSimplePaths() {
  deque<Sequence> contigs;
  deque<ContigInfo> contig_infos;
  Assemble(contigs, contig_infos);
  Initialize(contigs, contig_infos);
}

int64_t ContigGraph::Assemble(deque<Sequence> &contigs,
                              deque<ContigInfo> &contig_infos) {
  contigs.clear();
  contig_infos.clear();

  for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i) {
    if (vertices_[i].contig().size() == kmer_size_ &&
        vertices_[i].contig().IsPalindrome()) {
      vertices_[i].status().Lock(1);

      Sequence contig = vertices_[i].contig();
      // ContigInfo contig_info(vertices_[i].kmer_count(),
      // vertices_[i].in_edges(), vertices_[i].out_edges());
      ContigInfo contig_info;
      contig_info.set_kmer_count(vertices_[i].kmer_count());
      contig_info.in_edges() = vertices_[i].in_edges();
      contig_info.out_edges() = vertices_[i].out_edges();

      contigs.push_back(contig);
      contig_infos.push_back(contig_info);
    }
  }

  for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i) {
    if (!vertices_[i].status().Lock(0)) continue;

    ContigGraphPath path;
    path.Append(ContigGraphVertexAdaptor(&vertices_[i]), 0);

    Sequence contig;
    ContigInfo contig_info;
    for (int strand = 0; strand < 2; ++strand) {
      while (true) {
        ContigGraphVertexAdaptor current = path.back();
        ContigGraphVertexAdaptor next;

        if (!GetNextVertexAdaptor(current, next)) break;

        if (IsPalindromeLoop(path, next)) break;

        if (IsLoop(path, next)) goto FAIL;

        if (!next.status().LockPreempt(0)) goto FAIL;

        path.Append(next, -kmer_size_ + 1);
      }

      path.ReverseComplement();
    }

    path.Assemble(contig, contig_info);
    contigs.push_back(contig);
    contig_infos.push_back(contig_info);
  FAIL:;
  }

  ClearStatus();

  return contigs.size();
}

struct SearchNode {
  ContigGraphVertexAdaptor node;
  int distance;
  int label;
};

void ContigGraph::BuildBeginIdbaKmerMap() {
  begin_kmer_map_.clear();
  begin_kmer_map_.reserve(vertices_.size() * 2);

  for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i) {
    for (int strand = 0; strand < 2; ++strand) {
      ContigGraphVertexAdaptor current(&vertices_[i], strand);
      IdbaKmer kmer = current.begin_kmer(kmer_size_);

      IdbaKmer key = kmer.unique_format();
      begin_kmer_map_[key] = i;
    }
  }
}

bool ContigGraph::CycleDetect(ContigGraphVertexAdaptor current,
                              map<int, int> &status) {
  if (status[current.id()] == 0) {
    bool flag = false;
    status[current.id()] = 1;
    for (int x = 0; x < 4; ++x) {
      if (current.out_edges()[x]) {
        if (CycleDetect(GetNeighbor(current, x), status)) flag = true;
      }
    }
    status[current.id()] = 2;
    return flag;
  } else if (status[current.id()] == 1)
    return true;
  else
    return false;
}

void ContigGraph::TopSortDFS(deque<ContigGraphVertexAdaptor> &order,
                             ContigGraphVertexAdaptor current,
                             map<int, int> &status) {
  if (status[current.id()] == 0) {
    status[current.id()] = 1;
    for (int x = 0; x < 4; ++x) {
      if (current.out_edges()[x])
        TopSortDFS(order, GetNeighbor(current, x), status);
    }
    order.push_back(current);
  }
}

int ContigGraph::GetDepth(ContigGraphVertexAdaptor current, int depth,
                          int &maximum, int min_length) {
  if (depth > maximum) maximum = depth;

  if (maximum >= min_length) return min_length;

  deque<ContigGraphVertexAdaptor> neighbors;
  GetNeighbors(current, neighbors);
  for (unsigned i = 0; i < neighbors.size(); ++i) {
    if (neighbors[i].status().IsDead()) continue;

    GetDepth(neighbors[i], depth - kmer_size_ + 1 + neighbors[i].contig_size(),
             maximum, min_length);
  }

  return min(maximum, min_length);
}
