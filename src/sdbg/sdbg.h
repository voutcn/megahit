//
// Created by vout on 11/5/18.
//

#ifndef MEGAHIT_SDBG_H
#define MEGAHIT_SDBG_H

#include "sdbg_def.h"
#include "sdbg_raw_content.h"

#include <cassert>
#include <vector>
#include <iostream>
#include "kmlib/kmrns.h"
#include "kmlib/bitvector.h"

using std::vector;

/**
 * Succicent De Bruijn graph
 */
class SDBG {
 public:
  SDBG() = default;
  ~SDBG() = default;

  void LoadFromFile(const char *dbg_name) {
    ReadSdbgFromFile(dbg_name, &content_);
    rs_is_tip_.Build(content_.tip.data(), content_.meta.size());
    rs_w_.Build(content_.w.data(), content_.meta.size());
    rs_last_.Build(content_.last.data(), content_.meta.size());
    invalid_.FromPtr(content_.meta.size(), content_.tip.data());

    // build f and prefix_look_up_
    prefix_look_up_.resize(content_.meta.bucket_size());
    std::fill(f_, f_ + kAlphabetSize + 2, 0);
    f_[0] = -1;
    for (auto it = content_.meta.sorted_begin(); it != content_.meta.sorted_end(); ++it) {
      f_[it->bucket_id / (content_.meta.bucket_size() / kAlphabetSize) + 2] += it->num_items;
      prefix_look_up_[it->bucket_id].first = it->accumulate_item_count;
      prefix_look_up_[it->bucket_id].second = it->accumulate_item_count + it->num_items - 1;
    }
    for (unsigned i = 2; i < kAlphabetSize + 2; ++i) {
      f_[i] += f_[i - 1];
    }

    for (unsigned i = 1; i < kAlphabetSize + 2; ++i) {
      rank_f_[i] = rs_last_.Rank(f_[i] - 1);
    }

    for (uint64_t i = 0; i < content_.meta.size(); ++i) {
      if (GetW(i) == 0) {
        SetInvalidEdge(i);
      }
    }
  }

  uint64_t size() const {
    return content_.meta.size();
  }

  size_t k() const {
    return content_.meta.k();
  }

  uint8_t GetW(uint64_t x) const {
    return content_.w[x];
  }

  bool IsLast(uint64_t x) const {
    return content_.last[x];
  }

  bool IsLastOrTip(uint64_t x) const {
    return ((content_.last.data()[x / 64] | content_.tip.data()[x / 64]) >> (x % 64)) & 1;
  }

  int64_t GetLastIndex(uint64_t x) const {
    return rs_last_.Succ(x);
  }

  uint8_t GetNodeLastChar(uint64_t x) const {
    for (uint8_t i = 1; i < kAlphabetSize + 2; ++i) {
      if (f_[i] > int64_t(x)) {
        return i - 1;
      }
    }
    return kAlphabetSize + 2;
  }

  bool IsValidEdge(uint64_t edge_id) const {
    return !invalid_.get(edge_id);
  }

  bool IsTip(uint64_t edge_id) const {
    return content_.tip[edge_id];
  }

  void SetValidEdge(uint64_t edge_id) {
    invalid_.unset(edge_id);
  }

  void SetInvalidEdge(uint64_t edge_id) {
    invalid_.set(edge_id);
  }

  int EdgeMultiplicity(uint64_t edge_id) const {
    if (content_.full_mul.size()) {
      return content_.full_mul[edge_id];
    }

    if (content_.small_mul[edge_id] != kSmallMulSentinel) {
      return content_.small_mul[edge_id];
    } else {
      return content_.large_mul.at(edge_id);
    }
  }

  uint64_t Forward(uint64_t edge_id) const { // the last edge edge_id points to
    uint8_t a = GetW(edge_id);
    if (a > 4) {
      a -= 4;
    }
    int64_t count_a = rs_w_.Rank(a, edge_id);
    return rs_last_.Select(rank_f_[a] + count_a - 1);
  }

  uint64_t Backward(uint64_t edge_id) const { // the first edge points to edge_id
    uint8_t a = GetNodeLastChar(edge_id);
    int64_t count_a = rs_last_.Rank(edge_id - 1) - rank_f_[a];
    return rs_w_.Select(a, count_a);
  }

  /**
   * Find the index given a sequence
   * @param seq the sequence encoded in [0-3]
   * @return the index of the sequence in the graph, -1 if not exists
   */
  int64_t IndexBinarySearch(const uint8_t *seq) const {
    uint64_t prefix = 0;
    for (uint64_t i = 0; (1U << (i * 2)) < prefix_look_up_.size(); ++i) {
      prefix = prefix * kAlphabetSize + seq[k() - 1 - i] - 1;
    }
    auto l = prefix_look_up_[prefix].first;
    auto r = prefix_look_up_[prefix].second;

    while (l <= r) {
      int cmp = 0;
      uint64_t mid = (l + r) / 2;
      uint64_t y = mid;

      for (int i = k() - 1; i >= 0; --i) {
        if (IsTip(y)) {
          uint32_t const *tip_node_seq = content_.tip_lables.data() +
              content_.meta.words_per_tip_label() * (rs_is_tip_.Rank(y) - 1);
          for (unsigned j = 0; j < static_cast<unsigned>(i); ++j) {
            auto c = ((tip_node_seq[j / kCharsPerLabelWord] >>
                                                            (kCharsPerLabelWord - 1 - j % kCharsPerLabelWord)
                                                                * kBitsPerChar) & 3) + 1;
            if (c < seq[i - j]) {
              cmp = -1;
              break;
            } else if (c > seq[i - j]) {
              cmp = 1;
              break;
            }
          }

          if (cmp == 0) {
            if (IsTip(mid)) {
              cmp = -1;
            } else {
              auto c = ((tip_node_seq[i / kCharsPerLabelWord] >>
                  (kCharsPerLabelWord - 1 - i % kCharsPerLabelWord) * kBitsPerChar) & 3) + 1;
              if (c < seq[0]) {
                cmp = -1;
                break;
              } else if (c > seq[0]) {
                cmp = 1;
                break;
              }
            }
          }
          break;
        }

        y = Backward(y);
        uint8_t c = GetW(y);

        if (c < seq[i]) {
          cmp = -1;
          break;
        } else if (c > seq[i]) {
          cmp = 1;
          break;
        }
      }

      if (cmp == 0) {
        return GetLastIndex(mid);
      } else if (cmp > 0) {
        r = mid - 1;
      } else {
        l = mid + 1;
      }
    }
    return -1;
  }
  /**
   * Fetch the label of an edge from the graph
   * @param id the index in the graph
   * @param seq the label will be written to this address
   * @return the length of label (always k)
   */
  uint32_t Label(uint64_t id, uint8_t *seq) const {
    int64_t x = id;
    for (int i = k() - 1; i >= 0; --i) {
      if (IsTip(x)) {
        uint32_t const *tip_node_seq = content_.tip_lables.data() +
            content_.meta.words_per_tip_label() * (rs_is_tip_.Rank(x) - 1);
        for (int j = 0; j <= i; ++j) {
          seq[i - j] = static_cast<uint8_t>(
              ((tip_node_seq[j / kCharsPerLabelWord] >>
              (kCharsPerLabelWord - 1 - j % kCharsPerLabelWord) * kBitsPerChar) & 3) + 1);
        }
        break;
      }

      x = Backward(x);
      seq[i] = GetW(x);
      assert(seq[i] > 0);

      if (seq[i] > 4) {
        seq[i] -= 4;
      }
    }
    return k();
  }

 private:
  static const uint8_t kFlagWriteOut = 0x1;
  static const uint8_t kFlagMustEq0 = 0x2;
  static const uint8_t kFlagMustEq1 = 0x4;
  /**
   * An internal function to collect incoming edges & in-degrees
   * @tparam flag
   * @param edge_id
   * @param incomings the incoming edges will be written in the address if kFlagWriteOut is set
   * @return in degree of the edge; -1 if edge or flag invalid
   */
  template<uint8_t flag = 0>
  int ComputeIncomings(uint64_t edge_id, uint64_t *incomings) const {
    if (!IsValidEdge(edge_id)) {
      return -1;
    }

    uint64_t first_income = Backward(edge_id);
    uint8_t c = GetW(first_income);
    int count_ones = IsLastOrTip(first_income);
    int indegree = IsValidEdge(first_income);

    if (flag & kFlagMustEq0) {
      if (indegree) return -1;
    }

    if (flag & kFlagWriteOut) {
      if (indegree > 0) {
        incomings[0] = first_income;
      }
    }

    for (uint64_t y = first_income + 1; count_ones < 5 && y < size(); ++y) {
      count_ones += IsLastOrTip(y);
      uint8_t cur_char = GetW(y);

      if (cur_char == c) {
        break;
      } else if (cur_char == c + 4 && IsValidEdge(y)) {
        if (flag & kFlagMustEq0) {
          return -1;
        } else if (flag & kFlagMustEq1) {
          if (indegree == 1) return -1;
        }
        if (flag & kFlagWriteOut) {
          incomings[indegree] = y;
        }
        ++indegree;
      }
    }
    return indegree;
  }
  /**
   * An internal function to collect outgoing edges & in-degrees
   * @tparam flag
   * @param edge_id
   * @param outgoings the outgoing edges will be written in the address if kFlagWriteOut is set
   * @return out degree of the edge; -1 if edge or flag invalid
   */
  template<uint8_t flag = 0>
  int ComputeOutgoings(uint64_t edge_id, uint64_t *outgoings) const {
    if (!IsValidEdge(edge_id)) {
      return -1;
    }
    uint64_t outdegree = 0;
    int64_t next_edge = Forward(edge_id);
    do {
      if (IsValidEdge(next_edge)) {
        if (flag & kFlagMustEq0) {
          return -1;
        } else if (flag & kFlagMustEq1) {
          if (outdegree == 1) return -1;
        }
        if (flag & kFlagWriteOut) {
          outgoings[outdegree] = next_edge;
        }
        ++outdegree;
      }
      --next_edge;
    } while (next_edge >= 0 && !IsLastOrTip(next_edge));

    return outdegree;
  }

 public:
  /**
   * the in-degree of a node/edge
   * @param edge_id
   * @return the in-degree. -1 if id invalid.
   */
  int EdgeIndegree(uint64_t edge_id) const {
    return ComputeIncomings(edge_id, nullptr);
  }
  /**
   * the out-degree of an edge
   * @param edge_id
   * @return the out-degree. -1 if id invalid
   */
  int EdgeOutdegree(uint64_t edge_id) const {
    return ComputeOutgoings(edge_id, nullptr);
  }
  /**
   * get all incoming edges of an edge
   * @param edge_id
   * @param incomings all incoming edges' id will be written here
   * @return in-degree
   */
  int IncomingEdges(uint64_t edge_id, uint64_t *incomings) const {
    return ComputeIncomings<kFlagWriteOut>(edge_id, incomings);
  }
  /**
   * get all outgoing edges of an edge
   * @param edge_id
   * @param outgoings all outgoing edges' id will be written here
   * @return out-degree
   */
  int OutgoingEdges(uint64_t edge_id, uint64_t *outgoings) const {
    return ComputeOutgoings<kFlagWriteOut>(edge_id, outgoings);
  }
  /**
   * a more efficient way to judge whether an edge's in-degree is 0
   * @param edge_id
   * @return true if the edge's in-degree is 0
   */
  bool EdgeIndegreeZero(uint64_t edge_id) const {
    return ComputeIncomings<kFlagMustEq0>(edge_id, nullptr) == 0;
  }
  /**
   * A more efficient way to judge whether an edge's out-degree is 0
   * @param edge_id
   * @return true if the edge's out-degree is 0
   */
  bool EdgeOutdegreeZero(uint64_t edge_id) const {
    return ComputeOutgoings<kFlagMustEq0>(edge_id, nullptr) == 0;
  }
  /**
   * @param edge_id
   * @return if the edge has only one outgoing edge, return that one; otherwise -1
   */
  int64_t UniqueNextEdge(uint64_t edge_id) const {
    uint64_t ret = 0;
    if (ComputeOutgoings<kFlagWriteOut | kFlagMustEq1>(edge_id, &ret) == 1) {
      return ret;
    } else {
      return -1;
    }
  }
  /**
   * @param edge_id
   * @return if the edge has only one incoming edge, return that one; otherwise -1
   */
  int64_t UniquePrevEdge(uint64_t edge_id) const {
    uint64_t ret = 0;
    if (ComputeIncomings<kFlagWriteOut | kFlagMustEq1>(edge_id, &ret) == 1) {
      return ret;
    } else {
      return -1;
    }
  }
  /**
   * @param edge_id
   * @return if the edge has only one incoming edge which is on a simple path, return that one; otherwise -1
   */
  int64_t PrevSimplePathEdge(uint64_t edge_id) const {
    int64_t prev_edge = UniquePrevEdge(edge_id);
    if (prev_edge != -1 && UniqueNextEdge(prev_edge) != -1) {
      return prev_edge;
    } else {
      return -1;
    }
  }
  /**
   * @param edge_id
   * @return if the edge has only one outgoing edge which is on a simple path, return that one; otherwise -1
   */
  int64_t NextSimplePathEdge(uint64_t edge_id) const {
    int64_t next_edge = UniqueNextEdge(edge_id);
    if (next_edge != -1 && UniquePrevEdge(next_edge) != -1) {
      return next_edge;
    } else {
      return -1;
    }
  }
  /**
   * @param edge_id
   * @return the index of the reverse-complement edge of the input edge
   */
  int64_t EdgeReverseComplement(uint64_t edge_id) const {
    if (!IsValidEdge(edge_id)) {
      return -1;
    }

    uint8_t seq[kMaxK + 1];
    assert(k() == Label(edge_id, seq));
    seq[k()] = GetW(edge_id);

    if (seq[k()] > 4) {
      seq[k()] -= 4;
    }

    int i, j;
    for (i = 0, j = k(); i < j; ++i, --j) {
      std::swap(seq[i], seq[j]);
      seq[i] = 5 - seq[i];
      seq[j] = 5 - seq[j];
    }
    if (i == j) {
      seq[i] = 5 - seq[i];
    }

    int64_t rev_node = IndexBinarySearch(seq);
    if (rev_node == -1) return -1;
    do {
      uint8_t edge_label = GetW(rev_node);
      if (edge_label == seq[k()] || edge_label - 4 == seq[k()]) {
        return rev_node;
      }
      --rev_node;
    } while (rev_node >= 0 && !IsLastOrTip(rev_node));

    return -1;
  }

  /**
   * free mulitiplicity of all edges to reduce memory
   * WARNING: use this with cautions
   * After that EdgeMultiplicty() are invalid
   */
  void FreeMultiplicity() {
    content_.large_mul = std::move(spp::sparse_hash_map<uint64_t, mul_t>());
    content_.small_mul = std::move(std::vector<small_mul_t>());
    content_.full_mul = std::move(std::vector<mul_t>());
  }

 private:
  SdbgRawContent content_;
  kmlib::AtomicBitVector<uint64_t> invalid_;
  std::vector<std::pair<int64_t, int64_t>> prefix_look_up_;
  int64_t f_[kAlphabetSize + 2]{};
  int64_t rank_f_[kAlphabetSize + 2]{}; // = rs_last_.Rank(f_[i] - 1)
  kmlib::RankAndSelect<kAlphabetSize, kWAlphabetSize> rs_w_;
  RankAndSelect1Bit rs_last_;
  Rank1Bit rs_is_tip_;
};

using SuccinctDBG = SDBG;

#endif //MEGAHIT_SDBG_H
