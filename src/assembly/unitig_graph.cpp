//
// Created by vout on 11/19/18.
//

#include "unitig_graph.h"
#include <omp.h>
#include <mutex>
#include "utils.h"

namespace assembly {

UnitigGraph::UnitigGraph(SuccinctDBG *sdbg)
    : sdbg_(sdbg), adapter_impl_(this), sudo_adapter_impl_(this) {
  id_map_.clear();
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
          assert(rc_end != SDBG::kNullID);
        }
      }

      if (will_be_added) {
        std::lock_guard<std::mutex> lk(path_lock);
        vertices_.emplace_back(cur_edge, edge_idx, rc_start, rc_end, depth, length);
        count_palindrome += cur_edge == rc_start;
        assert(sdbg_->IsValidEdge(cur_edge));
        assert(sdbg_->IsValidEdge(edge_idx));
        assert(sdbg_->IsValidEdge(rc_start));
        assert(sdbg_->IsValidEdge(rc_end));
        assert(VertexAdapter(vertices_.back()).begin() == cur_edge);
        assert(VertexAdapter(vertices_.back()).end() == edge_idx);
        assert(VertexAdapter(vertices_.back()).rbegin() == rc_start);
        assert(VertexAdapter(vertices_.back()).rend() == rc_end);
        assert(sdbg_->EdgeOutdegree(edge_idx) != 1);
        assert(sdbg_->EdgeOutdegree(rc_end) != 1);
      }
    }
  }
  xinfo("Graph size without loops: %lu, palindrome: %lu\n", vertices_.size(), count_palindrome);

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
  id_map_.reserve(vertices_.size() * 2);

  for (size_type i = 0; i < vertices_.size(); ++i) {
    VertexAdapter adapter(vertices_[i]);
    id_map_[adapter.begin()] = i;
    id_map_[adapter.rbegin()] = i;
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

    uint8_t to_disconnect = adapter.is_to_disconnect();
    adapter.ReverseComplement();
    uint8_t rc_to_disconnect = adapter.is_to_disconnect();
    adapter.ReverseComplement();

    if (!to_disconnect && !rc_to_disconnect) {
      continue;
    }

    if (adapter.length() <= to_disconnect + rc_to_disconnect) {
      adapter.MarkToDelete();
      continue;
    }

    auto old_start = adapter.begin();
    auto old_end = adapter.end();
    auto old_rc_start = adapter.rbegin();
    auto old_rc_end = adapter.rend();
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
    uint64_t new_total_depth = adapter.avg_depth() * new_length + 0.5;
    vertices_[i] = UnitigGraphVertex(new_start, new_end, new_rc_start, new_rc_end,
                                     new_total_depth, new_length);

    std::lock_guard<std::mutex> lk(mutex);
    if (to_disconnect) {
      id_map_.erase(old_start);
      id_map_[new_start] = i;
    }
    if (rc_to_disconnect) {
      id_map_.erase(old_rc_start);
      id_map_[new_rc_start] = i;
    }
  }
}

void UnitigGraph::Refresh(bool set_changed) {
  static const uint8_t kDeleted = 0x1;
  static const uint8_t kVisited = 0x2;
  RefreshDisconnected();
#pragma omp parallel for
  for (size_type i = 0; i < vertices_.size(); ++i) {
    auto adapter = MakeSudoAdapter(i);
    if (!adapter.is_to_delete()) {
      continue;
    }
    for (int strand = 0; strand < 2; ++strand, adapter.ReverseComplement()) {
      uint64_t cur_edge = adapter.end(), prev_edge;
      while (cur_edge != adapter.begin()) {
        prev_edge = sdbg_->UniquePrevEdge(cur_edge);
        sdbg_->SetInvalidEdge(cur_edge);
        cur_edge = prev_edge;
        assert(cur_edge != SDBG::kNullID);
      }
      sdbg_->SetInvalidEdge(cur_edge);
      if (adapter.is_palindrome()) {
        break;
      }
    }
    adapter.set_flag(kDeleted);
  }
#pragma omp parallel for
  for (size_type i = 0; i < vertices_.size(); ++i) {
    MakeSudoAdapter(i).clear_cache();
  }

  locks_.reset();
#pragma omp parallel for
  for (size_type i = 0; i < vertices_.size(); ++i) {
    auto adapter = MakeSudoAdapter(i);
    if (adapter.flag() & kDeleted) {
      continue;
    }
    for (int strand = 0; strand < 2; ++strand, adapter.ReverseComplement()) {
      if (PrevSimplePathAdapter(adapter).valid()) {
        continue;
      }
      if (!locks_.try_lock(i)) {
        break;
      }
      std::vector<SudoVertexAdapter> linear_path;
      for (auto cur = NextSimplePathAdapter(adapter); cur.valid(); cur = NextSimplePathAdapter(cur)) {
        linear_path.emplace_back(cur);
      }

      if (linear_path.empty()) {
        adapter.set_flag(kVisited);
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
      adapter.set_flag(kVisited);

      for (auto &v: linear_path) {
        new_length += v.length();
        new_total_depth += v.total_depth();
        v.set_flag(kDeleted);
      }

      auto new_start = adapter.begin();
      auto new_rc_end = adapter.rend();
      auto new_rc_start = linear_path.back().rbegin();
      auto new_end = linear_path.back().end();

      adapter.set_start_end(new_start, new_end, new_rc_start, new_rc_end);
      adapter.set_length(new_length);
      adapter.set_total_depth(new_total_depth);
      if (set_changed) adapter.set_changed();
      break;
    }
  }

  // looped path
  std::mutex mutex;
#pragma omp parallel for
  for (size_type i = 0; i < vertices_.size(); ++i) {
    auto adapter = MakeSudoAdapter(i);
    if (adapter.flag()) {
      std::lock_guard<std::mutex> lk(mutex);
      if (adapter.flag()) {
        continue;
      }

      SudoVertexAdapter next_adapter = adapter;
      uint32_t length = adapter.length();
      uint64_t total_depth = adapter.total_depth();
      while ((next_adapter = NextSimplePathAdapter(next_adapter)).begin() != adapter.begin()) {
        next_adapter.set_flag(kDeleted);
        length += next_adapter.length();
        total_depth += next_adapter.total_depth();
      }

      auto prev_adapter = PrevSimplePathAdapter(adapter);
      auto new_start = adapter.begin();
      auto new_end = prev_adapter.end();
      auto new_rc_start = prev_adapter.rbegin();
      auto new_rc_end = adapter.rend();

      adapter.set_start_end(new_start, new_end, new_rc_start, new_rc_end);
      adapter.set_length(length);
      adapter.set_total_depth(total_depth);
      adapter.set_loop();
      if (set_changed) adapter.set_changed();
    }
  }

  vertices_.resize(
      std::remove_if(
          vertices_.begin(), vertices_.end(),
          [](UnitigGraphVertex &a) { return SudoVertexAdapter(a).flag() & kDeleted; }
      ) - vertices_.begin());

#pragma omp parallel for
  for (size_type i = 0; i < vertices_.size(); ++i) {
    auto adapter = MakeSudoAdapter(i);
    adapter.set_flag(0);
    assert(id_map_.count(adapter.begin()));
    assert(id_map_.count(adapter.rbegin()));
    id_map_.at(adapter.begin()) = i;
    id_map_.at(adapter.rbegin()) = i;
  }
}

std::string UnitigGraph::VertexToDNAString(const VertexAdapter &v) {
  std::string label;
  label.reserve(k() + v.length());
  uint64_t cur_edge = v.end();

  for (unsigned i = 1; i < v.length(); ++i) {
    int8_t cur_char = sdbg_->GetW(cur_edge);
    label.push_back("ACGT"[cur_char > 4 ? (cur_char - 5) : (cur_char - 1)]);

    cur_edge = sdbg_->PrevSimplePathEdge(cur_edge);
    if (cur_edge == SDBG::kNullID) {
      xfatal("%lld, %lld, %lld, %lld, (%lld, %lld), %d, %d\n",
             v.begin(),
             v.end(),
             v.rbegin(),
             v.rend(),
             sdbg_->EdgeReverseComplement(v.end()),
             sdbg_->EdgeReverseComplement(v.begin()),
             v.length(),
             i);
    }
  }

  int8_t cur_char = sdbg_->GetW(cur_edge);
  label.push_back("ACGT"[cur_char > 4 ? (cur_char - 5) : (cur_char - 1)]);

  if (cur_edge != v.begin()) {
    xfatal("fwd: %lld, %lld, rev: %lld, %lld, (%lld, %lld) length: %d\n",
           v.begin(),
           v.end(),
           v.rbegin(),
           v.rend(),
           sdbg_->EdgeReverseComplement(v.end()),
           sdbg_->EdgeReverseComplement(v.begin()),
           v.length());
  }

  uint8_t seq[kMaxK];
  sdbg_->Label(v.begin(), seq);

  for (int i = sdbg_->k() - 1; i >= 0; --i) {
    assert(seq[i] >= 1 && seq[i] <= 4);
    label.append(1, "ACGT"[seq[i] - 1]);
  }

  std::reverse(label.begin(), label.end());
  return label;
}

}