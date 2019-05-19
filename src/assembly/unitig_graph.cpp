//
// Created by vout on 11/19/18.
//

#include "unitig_graph.h"
#include <omp.h>
#include <mutex>
#include <cmath>

#include "kmlib/kmbitvector.h"
#include "utils/utils.h"

UnitigGraph::UnitigGraph(SDBG *sdbg)
    : sdbg_(sdbg), adapter_impl_(this), sudo_adapter_impl_(this) {
  id_map_.clear();
  vertices_.clear();
  std::mutex path_lock;
  AtomicBitVector locks(sdbg_->size());
  size_t count_palindrome = 0;
  // assemble simple paths
#pragma omp parallel for reduction(+: count_palindrome)
  for (uint64_t edge_idx = 0; edge_idx < sdbg_->size(); ++edge_idx) {
    if (sdbg_->IsValidEdge(edge_idx) && sdbg_->NextSimplePathEdge(edge_idx) == SDBG::kNullID
        && locks.try_lock(edge_idx)) {
      bool will_be_added = true;
      uint64_t cur_edge = edge_idx;
      uint64_t prev_edge;
      int64_t depth = sdbg_->EdgeMultiplicity(edge_idx);
      uint32_t length = 1;

      while ((prev_edge = sdbg_->PrevSimplePathEdge(cur_edge)) != SDBG::kNullID) {
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
        while ((rc_cur_edge = sdbg_->NextSimplePathEdge(rc_cur_edge)) != SDBG::kNullID) {
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
        std::lock_guard<std::mutex> lk(path_lock);
        vertices_.emplace_back(cur_edge, edge_idx, rc_start, rc_end, depth, length);
        count_palindrome += cur_edge == rc_start;
      }
    }
  }
  xinfo("Graph size without loops: %lu, palindrome: %lu\n", vertices_.size(), count_palindrome);

  // assemble looped paths
#pragma omp parallel for
  for (size_t edge_idx = 0; edge_idx < sdbg_->size(); ++edge_idx) {
    if (!locks.at(edge_idx) && sdbg_->IsValidEdge(edge_idx)) {
      std::lock_guard<std::mutex> lk(path_lock);
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
    id_map_[adapter.Begin()] = i;
    id_map_[adapter.RevBegin()] = i;
  }
}

void UnitigGraph::RefreshDisconnected() {
  std::mutex mutex;
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

    if (adapter.Length() <= to_disconnect + rc_to_disconnect) {
      adapter.SetToDelete();
      continue;
    }

    auto old_start = adapter.Begin();
    auto old_end = adapter.End();
    auto old_rc_start = adapter.RevBegin();
    auto old_rc_end = adapter.RevEnd();
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

    uint32_t new_length = adapter.Length() - to_disconnect - rc_to_disconnect;
    uint64_t new_total_depth = lround(adapter.AvgDepth() * new_length);
    adapter.SetBeginEnd(new_start, new_end, new_rc_start, new_rc_end);
    adapter.SetLength(new_length);
    adapter.SetTotalDepth(new_total_depth);

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
    if (!adapter.IsToDelete()) {
      continue;
    }
    adapter.SetFlag(kDeleted);
    if (adapter.ForSureStandalone()) {
      continue;
    }
    for (int strand = 0; strand < 2; ++strand, adapter.ReverseComplement()) {
      uint64_t cur_edge = adapter.End();
      for (size_t j = 1; j < adapter.Length(); ++j) {
        auto prev = sdbg_->UniquePrevEdge(cur_edge);
        sdbg_->SetInvalidEdge(cur_edge);
        cur_edge = prev;
        assert(cur_edge != SDBG::kNullID);
      }
      assert(cur_edge == adapter.Begin());
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
    if (adapter.ForSureStandalone() || (adapter.GetFlag() & kDeleted)) {
      continue;
    }
    for (int strand = 0; strand < 2; ++strand, adapter.ReverseComplement()) {
      if (PrevSimplePathAdapter(adapter).Valid()) {
        continue;
      }
      if (!locks.try_lock(i)) {
        break;
      }
      std::vector<SudoVertexAdapter> linear_path;
      for (auto cur = NextSimplePathAdapter(adapter); cur.Valid(); cur = NextSimplePathAdapter(cur)) {
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

      auto new_length = adapter.Length();
      auto new_total_depth = adapter.TotalDepth();
      adapter.SetFlag(kVisited);

      for (auto &v: linear_path) {
        new_length += v.Length();
        new_total_depth += v.TotalDepth();
        if (v.SdbgId() != adapter.SdbgId())
          v.SetFlag(kDeleted);
      }

      auto new_start = adapter.Begin();
      auto new_rc_end = adapter.RevEnd();
      auto new_rc_start = linear_path.back().RevBegin();
      auto new_end = linear_path.back().End();

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
    if (!adapter.ForSureStandalone() && !adapter.GetFlag()) {
      std::lock_guard<std::mutex> lk(mutex);
      if (adapter.GetFlag()) {
        continue;
      }

      uint32_t length = adapter.Length();
      uint64_t total_depth = adapter.TotalDepth();
      SudoVertexAdapter next_adapter = adapter;
      while (true) {
        next_adapter = NextSimplePathAdapter(next_adapter);
        assert(next_adapter.Valid());
        assert(!(next_adapter.GetFlag() & kDeleted));
        if (next_adapter.Begin() == adapter.Begin()) {
          break;
        }
        next_adapter.SetFlag(kDeleted);
        length += next_adapter.Length();
        total_depth += next_adapter.TotalDepth();
      }

      auto new_start = adapter.Begin();
      auto new_end = sdbg_->PrevSimplePathEdge(new_start);
      auto new_rc_end = adapter.RevEnd();
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

  vertices_.resize(
      std::remove_if(
          vertices_.begin(), vertices_.end(),
          [](UnitigGraphVertex &a) { return SudoVertexAdapter(a).GetFlag() & kDeleted; }
      ) - vertices_.begin());

  size_type num_changed = 0;
#pragma omp parallel for reduction(+: num_changed)
  for (size_type i = 0; i < vertices_.size(); ++i) {
    auto adapter = MakeSudoAdapter(i);
    assert(adapter.ForSureStandalone() || adapter.GetFlag());
    adapter.SetFlag(0);
    id_map_.at(adapter.Begin()) = i;
    id_map_.at(adapter.RevBegin()) = i;
    num_changed += adapter.IsChanged();
  }
}

std::string UnitigGraph::VertexToDNAString(VertexAdapter v) {
  v.ToUniqueFormat();
  std::string label;
  label.reserve(k() + v.Length());
  uint64_t cur_edge = v.End();

  for (unsigned i = 1; i < v.Length(); ++i) {
    int8_t cur_char = sdbg_->GetW(cur_edge);
    label.push_back("ACGT"[cur_char > 4 ? (cur_char - 5) : (cur_char - 1)]);

    cur_edge = sdbg_->PrevSimplePathEdge(cur_edge);
    if (cur_edge == SDBG::kNullID) {
      xfatal("%lld, %lld, %lld, %lld, (%lld, %lld), %d, %d\n",
             v.Begin(),
             v.End(),
             v.RevBegin(),
             v.RevEnd(),
             sdbg_->EdgeReverseComplement(v.End()),
             sdbg_->EdgeReverseComplement(v.Begin()),
             v.Length(),
             i);
    }
  }

  int8_t cur_char = sdbg_->GetW(cur_edge);
  label.push_back("ACGT"[cur_char > 4 ? (cur_char - 5) : (cur_char - 1)]);

  if (cur_edge != v.Begin()) {
    xfatal("fwd: %lld, %lld, rev: %lld, %lld, (%lld, %lld) length: %d\n",
           v.Begin(),
           v.End(),
           v.RevBegin(),
           v.RevEnd(),
           sdbg_->EdgeReverseComplement(v.End()),
           sdbg_->EdgeReverseComplement(v.Begin()),
           v.Length());
  }

  uint8_t seq[kMaxK];
  sdbg_->Label(v.Begin(), seq);

  for (int i = sdbg_->k() - 1; i >= 0; --i) {
    assert(seq[i] >= 1 && seq[i] <= 4);
    label.append(1, "ACGT"[seq[i] - 1]);
  }

  std::reverse(label.begin(), label.end());
  return label;
}
