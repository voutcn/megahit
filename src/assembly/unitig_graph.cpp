//
// Created by vout on 11/19/18.
//

#include "unitig_graph.h"
#include <omp.h>
#include <cmath>

#include "kmlib/kmbitvector.h"
#include "utils/mutex.h"
#include "utils/utils.h"

UnitigGraph::UnitigGraph(SDBG *sdbg)
    : sdbg_(sdbg), adapter_impl_(this), sudo_adapter_impl_(this) {
  id_map_.clear();
  vertices_.clear();
  SpinLock path_lock;
  AtomicBitVector locks(sdbg_->size());
  size_t count_palindrome = 0;
// assemble simple paths
#pragma omp parallel for reduction(+ : count_palindrome)
  for (uint64_t edge_idx = 0; edge_idx < sdbg_->size(); ++edge_idx) {
    if (sdbg_->IsValidEdge(edge_idx) &&
        sdbg_->NextSimplePathEdge(edge_idx) == SDBG::kNullID &&
        locks.try_lock(edge_idx)) {
      bool will_be_added = true;
      uint64_t cur_edge = edge_idx;
      uint64_t prev_edge;
      int64_t depth = sdbg_->EdgeMultiplicity(edge_idx);
      uint32_t length = 1;

      while ((prev_edge = sdbg_->PrevSimplePathEdge(cur_edge)) !=
             SDBG::kNullID) {
        cur_edge = prev_edge;
        if (!locks.try_lock(cur_edge)) {
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

      if (!locks.try_lock(rc_start)) {
        rc_end = sdbg_->EdgeReverseComplement(cur_edge);
        if (std::max(edge_idx, cur_edge) < std::max(rc_start, rc_end)) {
          will_be_added = false;
        }
      } else {
        // lock through the rc path
        uint64_t rc_cur_edge = rc_start;
        rc_end = rc_cur_edge;
        bool extend_full = true;
        while ((rc_cur_edge = sdbg_->NextSimplePathEdge(rc_cur_edge)) !=
               SDBG::kNullID) {
          rc_end = rc_cur_edge;
          if (!locks.try_lock(rc_cur_edge)) {
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
        std::lock_guard<SpinLock> lk(path_lock);
        vertices_.emplace_back(cur_edge, edge_idx, rc_start, rc_end, depth,
                               length);
        count_palindrome += cur_edge == rc_start;
      }
    }
  }
  xinfo("Graph size without loops: {}, palindrome: {}\n", vertices_.size(),
        count_palindrome);

  // assemble looped paths
  std::mutex loop_lock;
  size_t count_loop = 0;
#pragma omp parallel for
  for (size_t edge_idx = 0; edge_idx < sdbg_->size(); ++edge_idx) {
    if (!locks.at(edge_idx) && sdbg_->IsValidEdge(edge_idx)) {
      std::lock_guard<std::mutex> lk(loop_lock);
      if (!locks.at(edge_idx)) {
        uint64_t cur_edge = edge_idx;
        uint64_t rc_edge = sdbg_->EdgeReverseComplement(edge_idx);
        uint64_t depth = sdbg_->EdgeMultiplicity(edge_idx);
        uint32_t length = 0;
        // whether it is marked before entering the loop
        bool rc_marked = locks.at(rc_edge);

        while (!locks.at(cur_edge)) {
          locks.set(cur_edge);
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
                                 sdbg_->EdgeReverseComplement(start), depth,
                                 length, true);
          count_loop += 1;
        }
      }
    }
  }

  if (vertices_.size() >= kMaxNumVertices) {
    xfatal(
        "Too many vertices in the unitig graph ({} >= {}), "
        "you may increase the kmer size to remove tons of erroneous kmers.\n",
        vertices_.size(), kMaxNumVertices);
  }

  sdbg_->FreeMultiplicity();
  id_map_.reserve(vertices_.size() * 2 - count_palindrome);

  for (size_type i = 0; i < vertices_.size(); ++i) {
    VertexAdapter adapter(vertices_[i]);
    id_map_[adapter.b()] = i;
    id_map_[adapter.rb()] = i;
  }
  assert(vertices_.size() * 2 - count_palindrome >= id_map_.size());
}

void UnitigGraph::RefreshDisconnected() {
  SpinLock mutex;
#pragma omp parallel for
  for (size_type i = 0; i < vertices_.size(); ++i) {
    auto adapter = MakeSudoAdapter(i);
    if (adapter.IsToDelete() || adapter.IsPalindrome() || adapter.IsLoop()) {
      continue;
    }

    uint8_t to_disconnect = adapter.IsToDisconnect();
    adapter.ReverseComplement();
    uint8_t rc_to_disconnect = adapter.IsToDisconnect();
    adapter.ReverseComplement();

    if (!to_disconnect && !rc_to_disconnect) {
      continue;
    }

    if (adapter.GetLength() <= to_disconnect + rc_to_disconnect) {
      adapter.SetToDelete();
      continue;
    }

    auto old_start = adapter.b();
    auto old_end = adapter.e();
    auto old_rc_start = adapter.rb();
    auto old_rc_end = adapter.re();
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

    uint32_t new_length =
        adapter.GetLength() - to_disconnect - rc_to_disconnect;
    uint64_t new_total_depth = lround(adapter.GetAvgDepth() * new_length);
    adapter.SetBeginEnd(new_start, new_end, new_rc_start, new_rc_end);
    adapter.SetLength(new_length);
    adapter.SetTotalDepth(new_total_depth);

    std::lock_guard<SpinLock> lk(mutex);
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
    if (!adapter.IsToDelete()) {
      continue;
    }
    adapter.SetFlag(kDeleted);
    if (adapter.IsStandalone()) {
      continue;
    }
    for (int strand = 0; strand < 2; ++strand, adapter.ReverseComplement()) {
      uint64_t cur_edge = adapter.e();
      for (size_t j = 1; j < adapter.GetLength(); ++j) {
        auto prev = sdbg_->UniquePrevEdge(cur_edge);
        sdbg_->SetInvalidEdge(cur_edge);
        cur_edge = prev;
        assert(cur_edge != SDBG::kNullID);
      }
      assert(cur_edge == adapter.b());
      sdbg_->SetInvalidEdge(cur_edge);
      if (adapter.IsPalindrome()) {
        break;
      }
    }
  }

  AtomicBitVector locks(size());
#pragma omp parallel for
  for (size_type i = 0; i < vertices_.size(); ++i) {
    auto adapter = MakeSudoAdapter(i);
    if (adapter.IsStandalone() || (adapter.GetFlag() & kDeleted)) {
      continue;
    }
    for (int strand = 0; strand < 2; ++strand, adapter.ReverseComplement()) {
      if (PrevSimplePathAdapter(adapter).IsValid()) {
        continue;
      }
      if (!locks.try_lock(i)) {
        break;
      }
      std::vector<SudoVertexAdapter> linear_path;
      for (auto cur = NextSimplePathAdapter(adapter); cur.IsValid();
           cur = NextSimplePathAdapter(cur)) {
        linear_path.emplace_back(cur);
      }

      if (linear_path.empty()) {
        adapter.SetFlag(kVisited);
        break;
      }

      size_type back_id = linear_path.back().UnitigId();
      if (back_id != i && !locks.try_lock(back_id)) {
        if (back_id > i) {
          locks.unlock(i);
          break;
        } else {
          locks.lock(back_id);
        }
      }

      auto new_length = adapter.GetLength();
      auto new_total_depth = adapter.GetTotalDepth();
      adapter.SetFlag(kVisited);

      for (auto &v : linear_path) {
        new_length += v.GetLength();
        new_total_depth += v.GetTotalDepth();
        if (v.canonical_id() != adapter.canonical_id()) v.SetFlag(kDeleted);
      }

      auto new_start = adapter.b();
      auto new_rc_end = adapter.re();
      auto new_rc_start = linear_path.back().rb();
      auto new_end = linear_path.back().e();

      adapter.SetBeginEnd(new_start, new_end, new_rc_start, new_rc_end);
      adapter.SetLength(new_length);
      adapter.SetTotalDepth(new_total_depth);
      if (set_changed) adapter.SetChanged();
      break;
    }
  }

  // looped path
  std::mutex mutex;
#pragma omp parallel for
  for (size_type i = 0; i < vertices_.size(); ++i) {
    auto adapter = MakeSudoAdapter(i);
    if (!adapter.IsStandalone() && !adapter.GetFlag()) {
      std::lock_guard<std::mutex> lk(mutex);
      if (adapter.GetFlag()) {
        continue;
      }

      uint32_t length = adapter.GetLength();
      uint64_t total_depth = adapter.GetTotalDepth();
      SudoVertexAdapter next_adapter = adapter;
      while (true) {
        next_adapter = NextSimplePathAdapter(next_adapter);
        assert(next_adapter.IsValid());
        if (next_adapter.b() == adapter.b()) {
          break;
        }
        next_adapter.SetFlag(kDeleted);
        length += next_adapter.GetLength();
        total_depth += next_adapter.GetTotalDepth();
      }

      auto new_start = adapter.b();
      auto new_end = sdbg_->PrevSimplePathEdge(new_start);
      auto new_rc_end = adapter.re();
      auto new_rc_start = sdbg_->NextSimplePathEdge(new_rc_end);
      assert(new_start == sdbg_->EdgeReverseComplement(new_rc_end));
      assert(new_end == sdbg_->EdgeReverseComplement(new_rc_start));

      adapter.SetBeginEnd(new_start, new_end, new_rc_start, new_rc_end);
      adapter.SetLength(length);
      adapter.SetTotalDepth(total_depth);
      adapter.SetLooped();
      if (set_changed) adapter.SetChanged();
    }
  }

  vertices_.resize(std::remove_if(vertices_.begin(), vertices_.end(),
                                  [](UnitigGraphVertex &a) {
                                    return SudoVertexAdapter(a).GetFlag() &
                                           kDeleted;
                                  }) -
                   vertices_.begin());

  size_type num_changed = 0;
#pragma omp parallel for reduction(+ : num_changed)
  for (size_type i = 0; i < vertices_.size(); ++i) {
    auto adapter = MakeSudoAdapter(i);
    assert(adapter.IsStandalone() || adapter.GetFlag());
    adapter.SetFlag(0);
    id_map_.at(adapter.b()) = i;
    id_map_.at(adapter.rb()) = i;
    num_changed += adapter.IsChanged();
  }
}

std::string UnitigGraph::VertexToDNAString(VertexAdapter v) {
  v.ToUniqueFormat();
  std::string label;
  label.reserve(k() + v.GetLength());
  uint64_t cur_edge = v.e();

  for (unsigned i = 1; i < v.GetLength(); ++i) {
    int8_t cur_char = sdbg_->GetW(cur_edge);
    label.push_back("ACGT"[cur_char > 4 ? (cur_char - 5) : (cur_char - 1)]);

    cur_edge = sdbg_->PrevSimplePathEdge(cur_edge);
    if (cur_edge == SDBG::kNullID) {
      xfatal("{}, {}, {}, {}, ({}, {}), {}, {}\n", v.b(), v.e(), v.rb(), v.re(),
             sdbg_->EdgeReverseComplement(v.e()),
             sdbg_->EdgeReverseComplement(v.b()), v.GetLength(), i);
    }
  }

  int8_t cur_char = sdbg_->GetW(cur_edge);
  label.push_back("ACGT"[cur_char > 4 ? (cur_char - 5) : (cur_char - 1)]);

  if (cur_edge != v.b()) {
    xfatal("fwd: {}, {}, rev: {}, {}, ({}, {}) length: {}\n", v.b(), v.e(),
           v.rb(), v.re(), sdbg_->EdgeReverseComplement(v.e()),
           sdbg_->EdgeReverseComplement(v.b()), v.GetLength());
  }

  uint8_t seq[kMaxK];
  sdbg_->GetLabel(v.b(), seq);

  for (int i = sdbg_->k() - 1; i >= 0; --i) {
    assert(seq[i] >= 1 && seq[i] <= 4);
    label.append(1, "ACGT"[seq[i] - 1]);
  }

  std::reverse(label.begin(), label.end());
  return label;
}
