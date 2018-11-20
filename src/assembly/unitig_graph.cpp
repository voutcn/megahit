//
// Created by vout on 11/19/18.
//

#include "unitig_graph.h"
#include <omp.h>
#include <mutex>
#include "utils.h"

namespace assembly {

UnitigGraph::UnitigGraph(SuccinctDBG *sdbg) : sdbg_(sdbg) {
  start_node_map_.clear();
  vertices_.clear();
  std::mutex path_lock;
  locks_.reset(sdbg_->size());
  size_t count_palindrome = 0;

  // assemble simple paths
#pragma omp parallel for reduction(+: count_palindrome)
  for (size_t edge_idx = 0; edge_idx < sdbg_->size(); ++edge_idx) {
    if (sdbg_->IsValidEdge(edge_idx) && sdbg_->NextSimplePathEdge(edge_idx) == SDBG::kNullID
        && locks_.try_lock(edge_idx)) {
      bool will_be_added = true;
      uint64_t cur_edge = edge_idx;
      uint64_t prev_edge;
      int64_t depth = sdbg_->EdgeMultiplicity(edge_idx);
      uint32_t length = 1;

      while ((prev_edge = sdbg_->PrevSimplePathEdge(cur_edge)) != SDBG::kNullID) {
        cur_edge = prev_edge;
        if (!locks_.try_lock(cur_edge)) {
          will_be_added = false;
          break;
        }
        depth += sdbg_->EdgeMultiplicity(cur_edge);
        ++length;
      }

      if (!will_be_added) {
        continue;
      }

      uint64_t rc_start = sdbg_->EdgeReverseComplement(edge_idx);
      uint64_t rc_end;
      assert(rc_start != SDBG::kNullID);

      if (!locks_.try_lock(rc_start)) {
        rc_end = sdbg_->EdgeReverseComplement(cur_edge);
        if (std::max(edge_idx, cur_edge) < std::max(rc_start, rc_end)) {
          will_be_added = false;
        }
      } else {
        // lock through the rc path
        uint64_t rc_cur_edge = rc_start;
        rc_end = rc_cur_edge;
        bool extend_full = true;
        while ((rc_cur_edge = sdbg_->NextSimplePathEdge(rc_cur_edge)) != SDBG::kNullID) {
          rc_end = rc_cur_edge;
          if (!locks_.try_lock(rc_cur_edge)) {
            extend_full = false;
            break;
          }
        }
        if (!extend_full) {
          rc_end = sdbg_->EdgeReverseComplement(cur_edge);
        }
      }

      if (will_be_added) {
        std::lock_guard<std::mutex> lk(path_lock);
        vertices_.emplace_back(cur_edge, edge_idx, rc_start, rc_end, depth, length);
        count_palindrome += vertices_.back().is_palindrome;
      }
    }
  }
  xlog("Graph size without loops: %lu, palindrome: %lu\n", vertices_.size(), count_palindrome);

  // assemble looped paths
#pragma omp parallel for
  for (size_t edge_idx = 0; edge_idx < sdbg_->size(); ++edge_idx) {
    if (!locks_.at(edge_idx) && sdbg_->IsValidEdge(edge_idx)) {
      std::lock_guard<std::mutex> lk(path_lock);
      if (!locks_.at(edge_idx)) {
        uint64_t cur_edge = edge_idx;
        uint64_t rc_edge = sdbg_->EdgeReverseComplement(edge_idx);
        uint64_t depth = sdbg_->EdgeMultiplicity(edge_idx);
        uint32_t length = 0;
        // whether it is marked before entering the loop
        bool rc_marked = locks_.at(rc_edge);

        while (!locks_.at(cur_edge)) {
          locks_.set(cur_edge);
          depth += sdbg_->EdgeMultiplicity(cur_edge);
          ++length;
          cur_edge = sdbg_->PrevSimplePathEdge(cur_edge);
          assert(cur_edge != SDBG::kNullID);
        }
        assert(cur_edge == edge_idx);

        if (!rc_marked) {
          uint64_t start = sdbg_->NextSimplePathEdge(edge_idx);
          uint64_t end = edge_idx;
          vertices_.emplace_back(start, end, sdbg_->EdgeReverseComplement(end),
                                 sdbg_->EdgeReverseComplement(start), depth, length, true);
          if (locks_.at(rc_edge)) {
            // this loop is palindrome
            vertices_.back().is_palindrome = true;
          }
        }
      }
    }
  }

  if (vertices_.size() >= kMaxNumVertices) {
    xfatal("Too many vertices in the unitig graph (%llu >= %llu), "
           "you may increase the kmer size to remove tons of erroneous kmers.\n",
           static_cast<unsigned long long>(vertices_.size()),
           static_cast<unsigned long long>(kMaxNumVertices));
  }

  sdbg_->FreeMultiplicity();

  start_node_map_.reserve(vertices_.size() * 2);
  for (size_type i = 0; i < vertices_.size(); ++i) {
    start_node_map_[vertices_[i].strand_info[0].start] = i;
    start_node_map_[vertices_[i].strand_info[1].start] = i;
  }
  locks_.reset(vertices_.size());
}

void UnitigGraph::RefreshDisconnected() {
  std::mutex mutex;
#pragma omp parallel for
  for (size_type i = 0; i < vertices_.size(); ++i) {
    auto adapter = MakeSudoAdapter(i);
    if (adapter.is_to_delete() || adapter.is_palindrome() || adapter.is_loop()) {
      continue;
    }

    uint8_t to_disconnect = adapter.is_to_disconnected();
    adapter.ReverseComplement();
    uint8_t rc_to_disconnect = adapter.is_to_disconnected();
    adapter.ReverseComplement();

    if (!to_disconnect && !rc_to_disconnect) {
      continue;
    }

    if (adapter.length() <= to_disconnect + rc_to_disconnect) {
      adapter.MarkToDelete();
      continue;
    }

    auto old_start = adapter.start();
    auto old_end = adapter.end();
    auto old_rc_start = adapter.rc_start();
    auto old_rc_end = adapter.rc_end();
    uint64_t new_start, new_end, new_rc_start, new_rc_end;

    if (to_disconnect) {
      new_start = sdbg_->NextSimplePathEdge(old_start);
      new_rc_end = sdbg_->PrevSimplePathEdge(old_rc_end);
      assert(new_start != SDBG::kNullID && new_rc_end != SDBG::kNullID);
      sdbg_->SetInvalidEdge(old_start);
      sdbg_->SetInvalidEdge(old_rc_end);
    } else {
      new_start = old_start;
      new_rc_end = old_rc_end;
    }

    if (rc_to_disconnect) {
      new_rc_start = sdbg_->NextSimplePathEdge(old_rc_start);
      new_end = sdbg_->PrevSimplePathEdge(old_end);
      assert(new_rc_start != SDBG::kNullID && new_end != SDBG::kNullID);
      sdbg_->SetInvalidEdge(old_rc_start);
      sdbg_->SetInvalidEdge(old_end);
    } else {
      new_rc_start = old_rc_start;
      new_end = old_end;
    }

    uint32_t new_length = adapter.length() - to_disconnect - rc_to_disconnect;
    uint64_t new_total_depth = adapter.depth() * new_length + 0.5;
    vertices_[i] = Vertex(new_start, new_end, new_rc_start, new_rc_end,
        new_total_depth, new_length);

    std::lock_guard<std::mutex> lk(mutex);
    if (to_disconnect) {
      start_node_map_.erase(old_start);
      start_node_map_[new_start] = i;
    }
    if (rc_to_disconnect) {
      start_node_map_.erase(old_rc_start);
      start_node_map_[new_rc_start] = i;
    }
  }
}

void UnitigGraph::Refresh(bool set_changed) {
  RefreshDisconnected();
#pragma omp parallel for
  for (size_type i = 0; i < vertices_.size(); ++i) {
    if (vertices_[i].to_delete && !vertices_[i].is_deleted) {
      auto adapter = MakeSudoAdapter(i);
      for (int strand = 0; strand < 2; ++strand) {
        uint64_t cur_edge = adapter.end(), prev_edge;
        while (cur_edge != adapter.start()) {
          prev_edge = sdbg_->UniquePrevEdge(cur_edge);
          sdbg_->SetInvalidEdge(cur_edge);
          cur_edge = prev_edge;
          assert(cur_edge != SDBG::kNullID);
        }
        sdbg_->SetInvalidEdge(cur_edge);
        if (adapter.is_palindrome()) {
          break;
        } else {
          adapter.ReverseComplement();
        }
      }
      adapter.MarkDeleted();
    }
  }
#pragma omp parallel for
  for (size_type i = 0; i < vertices_.size(); ++i) {
    for (auto &strand : vertices_[i].strand_info) {
      strand.cached_out_degree = Vertex::kUnknownDegree;
    }
  }

  locks_.reset();
#pragma omp parallel for
  for (size_type i = 0; i < vertices_.size(); ++i) {
    auto adapter = MakeSudoAdapter(i);
    if (adapter.is_deleted()) {
      continue;
    }
    for (int strand = 0; strand < 2; ++strand) {
      if (adapter.PrevSimplePathAdapter().valid()) {
        adapter.ReverseComplement();
        continue;
      }
      if (!locks_.try_lock(i)) {
        break;
      }
      std::vector<SudoVertexAdapter> linear_path;
      for (auto cur = adapter.NextSimplePathAdapter(); cur.valid(); cur = cur.NextSimplePathAdapter()) {
        linear_path.emplace_back(cur);
      }

      if (linear_path.empty()) {
        adapter.SetFlag(1);
        break;
      }

      if (linear_path.back().id() != i && !locks_.try_lock(linear_path.back().id())) {
        if (linear_path.back().id() > i) {
          locks_.unlock(i);
          break;
        } else {
          locks_.lock(linear_path.back().id());
        }
      }

      auto new_length = adapter.length();
      auto new_total_depth = adapter.total_depth();
      for (auto &v: linear_path) {
        new_length += v.length();
        new_total_depth += v.total_depth();
        v.MarkDeleted();
      }

      auto new_start = adapter.start();
      auto new_rc_end = adapter.rc_end();
      auto new_rc_start = linear_path.back().rc_start();
      auto new_end = linear_path.back().end();

      vertices_[i] = Vertex(new_start, new_end, new_rc_start, new_rc_end, new_total_depth, new_length);
      vertices_[i].is_changed |= set_changed;
      vertices_[i].flag = 1;
      break;
    }
  }

  // looped path
  std::mutex mutex;
#pragma omp parallel for
  for (size_type i = 0; i < vertices_.size(); ++i) {
    if (!vertices_[i].is_deleted && !vertices_[i].flag) {
      std::lock_guard<std::mutex> lk(mutex);
      auto adapter = MakeSudoAdapter(i);
      if (adapter.is_deleted()) {
        continue;
      }
      uint32_t length = adapter.length();
      uint64_t total_depth = adapter.total_depth();

      adapter.MarkDeleted();
      bool is_palindrome = false;

      for (auto next = adapter.NextSimplePathAdapter(); next.start() != adapter.start();
           next = next.NextSimplePathAdapter()) {
        assert(next.valid());
        length += next.length();
        total_depth += next.total_depth();
        is_palindrome |= next.is_deleted();
        next.MarkDeleted();
      }

      auto new_start = adapter.start();
      auto new_end = sdbg_->PrevSimplePathEdge(new_start);
      auto new_rc_start = sdbg_->EdgeReverseComplement(new_end);
      auto new_rc_end = sdbg_->EdgeReverseComplement(new_start);

      vertices_[i] = Vertex(new_start, new_end, new_rc_start, new_rc_end, total_depth, length, true);
      vertices_[i].is_palindrome = is_palindrome;
      vertices_[i].is_deleted = true;
      vertices_[i].is_changed |= set_changed;
    }
  }

  vertices_.resize(
      std::remove_if(
          vertices_.begin(), vertices_.end(),
          [](const Vertex &a) { return a.is_deleted; }) - vertices_.begin());

#pragma omp parallel for
  for (size_type i = 0; i < vertices_.size(); ++i) {
    assert(vertices_[i].is_loop || vertices_[i].flag);
    vertices_[i].flag = 0;
    start_node_map_.at(vertices_[i].strand_info[0].start) = i;
    start_node_map_.at(vertices_[i].strand_info[1].start) = i;
  }
}

}