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

#ifndef SUCCINCT_DBG_H_
#define SUCCINCT_DBG_H_
#include <assert.h>
#include <vector>
#include "definitions.h"
#include "kmlib/kmrns.h"
#include "sdbg_multi_io.h"
#include "kmlib/bitvector.h"

#ifdef USE_KHASH
#include "khash.h"
#else
#include "sparsepp/sparsepp/spp.h"
#endif

using std::vector;

/**
 * Succicent De Bruijn graph
 * @tparam mul_type the type of multiplicity
 * @tparam small_mul_type the type of small multiplicity
 */
template<typename mul_type = uint16_t, typename small_mul_type = uint8_t>
class SDBG {
 public:
#ifdef USE_KHASH
  KHASH_MAP_INIT_INT64(k64v16, mul_type); // declare khash
#endif
  // constants
  static const unsigned kAlphabetSize = 4;
  static const unsigned kWAlphabetSize = 9;
  static const unsigned kWBitsPerChar = 4;
  static const unsigned kWCharsPerWord = sizeof(unsigned long long) * 8 / kWBitsPerChar;
  static const unsigned kWCharMask = 0xF;
  static const unsigned kMaxKmerK = kMaxK + 1;
  static const unsigned kCharsPerUint32 = 16;
  static const unsigned kBitsPerChar = 2;
  static const unsigned kBitsPerULL = 64;

 public:
  SDBG() = default;
  ~SDBG() = default;

  void LoadFromFile(const char *dbg_name) {
    SdbgReader sdbg_reader;
    sdbg_reader.set_file_prefix(std::string(dbg_name));
    sdbg_reader.read_info();
    sdbg_reader.init_files();

    for (int i = 0; i < 6; ++i) {
      f_[i] = sdbg_reader.f()[i];
    }

    k_ = sdbg_reader.kmer_size();
    size_ = sdbg_reader.num_items();
    num_tip_nodes_ = sdbg_reader.num_tips();
    uint32_per_tip_nodes_ = sdbg_reader.words_per_tip_label();
    prefix_lk_len_ = sdbg_reader.prefix_lkt_len();

    size_t word_needed_w = (size_ + kWCharsPerWord - 1) / kWCharsPerWord;
    size_t word_needed_last = (size_ + kBitsPerULL - 1) / kBitsPerULL;

    w_.resize(word_needed_w);
    last_.resize(word_needed_last);
    is_tip_.resize(word_needed_last);
    tip_node_seq_.resize(size_t(num_tip_nodes_) * sdbg_reader.words_per_tip_label());
    prefix_lkt_.resize(sdbg_reader.prefix_lkt_size() * 2);

    if (sdbg_reader.num_large_mul() > (1 << 30) ||
        sdbg_reader.num_large_mul() > size_ * 0.08) {
      large_mul_.resize(size_);
    } else {
      small_mul_.resize(size_);
#ifdef USE_KHASH
      kh_ = kh_init(k64v16);
#else
      large_mul_lookup_.reserve(sdbg_reader.num_large_mul());
#endif
    }

    for (int i = 0; i < sdbg_reader.prefix_lkt_size() * 2; ++i) {
      prefix_lkt_[i] = sdbg_reader.prefix_lkt(i);
    }

    ull_t packed_w = 0;
    ull_t packed_last = 0;
    ull_t packed_tip = 0;

    int64_t w_word_idx = 0;
    int64_t last_word_idx = 0;

    int w_word_offset = 0;
    int last_word_offset = 0;

    uint64_t tip_label_offset = 0;
    uint16_t item = 0;

    for (uint64_t i = 0; LIKELY(i < size_); ++i) {
      assert(sdbg_reader.NextItem(item));

      packed_w |= (unsigned long long) (item & 0xF) << w_word_offset;
      w_word_offset += kWBitsPerChar;

      if (w_word_offset == kBitsPerULL) {
        w_[w_word_idx++] = packed_w;
        w_word_offset = 0;
        packed_w = 0;
      }

      packed_last |= (unsigned long long) ((item >> 4) & 1) << last_word_offset;
      packed_tip |= (unsigned long long) ((item >> 5) & 1) << last_word_offset;

      last_word_offset++;

      if (last_word_offset == kBitsPerULL) {
        last_[last_word_idx] = packed_last;
        is_tip_[last_word_idx++] = packed_tip;
        last_word_offset = 0;
        packed_tip = packed_last = 0;
      }

      if (!small_mul_.empty())
        small_mul_[i] = item >> 8;
      else
        large_mul_[i] = item >> 8;

      if (UNLIKELY((item >> 8) == kMulti2Sp)) {
        multi_t mul = sdbg_reader.NextLargeMul();
        assert(mul >= kMulti2Sp);
        if (!small_mul_.empty()) {
#ifdef USE_KHASH
          int ret;
          khint_t k = kh_put(k64v16, kh_, i, &ret);
          kh_value(kh_, k) = mul;
#else
          large_mul_lookup_[i] = mul;
#endif
        } else {
          large_mul_[i] = mul;
        }
      }

      if (UNLIKELY((item >> 5) & 1)) {
        sdbg_reader.NextTipLabel(&tip_node_seq_[tip_label_offset]);
        tip_label_offset += sdbg_reader.words_per_tip_label();
      }
    }

    if (w_word_offset != 0) {
      w_[w_word_idx++] = packed_w;
    }

    if (last_word_offset != 0) {
      last_[last_word_idx] = packed_last;
      is_tip_[last_word_idx++] = packed_tip;
    }

    assert(!sdbg_reader.NextItem(item));
    assert(tip_label_offset == num_tip_nodes_ * sdbg_reader.words_per_tip_label());

    invalid_ = is_tip_;
    rs_is_tip_.Build(&is_tip_[0], size_);
    rs_w_.Build(&w_[0], size_);
    rs_last_.Build(&last_[0], size_);

    for (unsigned i = 1; i < kAlphabetSize + 2; ++i) {
      rank_f_[i] = rs_last_.Rank(f_[i] - 1);
    }

    for (uint64_t i = 0; i < size_; ++i) {
      if (GetW(i) == 0) {
        SetInvalidEdge(i);
      }
    }

    in_or_out_zero.reset(size_);
  }

  uint64_t size() const {
    return size_;
  }

  uint32_t k() const {
    return k_;
  }

  uint8_t GetW(uint64_t x) const {
    return (w_[x / kWCharsPerWord] >> (x % kWCharsPerWord * kWBitsPerChar)) & kWCharMask;
  }

  bool IsLast(uint64_t x) const {
    return (last_[x / 64] >> (x % 64)) & 1;
  }

  bool IsLastOrTip(uint64_t x) const {
    return ((last_[x / 64] | is_tip_[x / 64]) >> (x % 64)) & 1;
  }

  int64_t GetLastIndex(uint64_t x) const {
    return rs_last_.Succ(x);
  }

  uint8_t GetNodeLastChar(uint64_t x) const {
    for (uint8_t i = 1; i < 6; ++i) {
      if (f_[i] > x) {
        return i - 1;
      }
    }
    assert(false);
    return 6;
  }

  bool IsValidEdge(uint64_t edge_id) const {
    return !((invalid_[edge_id / 64] >> (edge_id % 64)) & 1);
  }

  bool IsTip(uint64_t edge_id) const {
    return (is_tip_[edge_id / 64] >> (edge_id % 64)) & 1;
  }

  void SetValidEdge(uint64_t edge_id) {
    __sync_fetch_and_and(&invalid_[edge_id / 64], ~(1ULL << (edge_id % 64)));
  }

  void SetInvalidEdge(uint64_t edge_id) {
    __sync_fetch_and_or(&invalid_[edge_id / 64], 1ULL << (edge_id % 64));
  }

  int EdgeMultiplicity(uint64_t edge_id) const {
    if (large_mul_.size()) {
      return large_mul_[edge_id];
    }

    if (LIKELY(small_mul_[edge_id] != kMulti2Sp)) {
      return small_mul_[edge_id];
    } else {
#ifdef USE_KHASH
      return kh_value(kh_, kh_get(k64v16, kh_, edge_id));
#else
      return large_mul_lookup_.at(edge_id);
#endif
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
    // only work if k > 8
    int pre = 0;
    for (uint32_t i = 0; i < prefix_lk_len_; ++i) {
      pre = pre * 4 + seq[k_ - 1 - i] - 1;
    }
    uint64_t l = prefix_lkt_[pre * 2];
    uint64_t r = prefix_lkt_[pre * 2 + 1];

    while (l <= r) {
      int cmp = 0;
      uint64_t mid = (l + r) / 2;
      uint64_t y = mid;

      for (int i = k_ - 1; i >= 0; --i) {
        if (IsTip(y)) {
          uint32_t const *tip_node_seq = &tip_node_seq_[uint32_per_tip_nodes_ * (rs_is_tip_.Rank(y) - 1)];

          for (int j = 0; j < i; ++j) {
            uint8_t
                c =
                (tip_node_seq[j / kCharsPerUint32] >> (kCharsPerUint32 - 1 - j % kCharsPerUint32) * kBitsPerChar) & 3;
            c++;

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
              uint8_t c =
                  (tip_node_seq[i / kCharsPerUint32] >> (kCharsPerUint32 - 1 - i % kCharsPerUint32) * kBitsPerChar) & 3;
              c++;

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

    for (int i = k_ - 1; i >= 0; --i) {
      if (IsTip(x)) {
        uint32_t const *tip_node_seq = &tip_node_seq_[uint32_per_tip_nodes_ * (rs_is_tip_.Rank(x) - 1)];

        for (int j = 0; j <= i; ++j) {
          seq[i - j] =
              (tip_node_seq[j / kCharsPerUint32] >> (kCharsPerUint32 - 1 - j % kCharsPerUint32) * kBitsPerChar) & 3;
          seq[i - j]++;
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

    return k_;
  }

 private:
  static const uint8_t kFlagWriteOut = 0x1;
  static const uint8_t kFlagMustEq0 = 0x2;
  static const uint8_t kFlagMustEq1 = 0x4;
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

    for (uint64_t y = first_income + 1; count_ones < 5 && y < size_; ++y) {
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
    uint64_t ret = -1;
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
    uint64_t ret = -1;
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

    uint8_t seq[kMaxKmerK + 1];
    assert(k_ == Label(edge_id, seq));
    seq[k_] = GetW(edge_id);

    if (seq[k_] > 4) {
      seq[k_] -= 4;
    }

    int i, j;

    for (i = 0, j = k_; i < j; ++i, --j) {
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

      if (edge_label == seq[k_] || edge_label - 4 == seq[k_]) {
        return rev_node;
      }

      --rev_node;
    } while (rev_node >= 0 && !IsLastOrTip(rev_node));

    return -1;
  }

  // WARNING: use this with cautions
  // After that EdgeMultiplicty() are invalid
  void FreeMultiplicity() {
    small_mul_ = std::move(std::vector<small_mul_type>());
    large_mul_ = std::move(std::vector<mul_type>());
#ifdef USE_KHASH
    kh_destroy(k64v16, kh_);
#else
    large_mul_lookup_ = std::move(spp::sparse_hash_map<uint64_t, multi_t>());
#endif
  }

 private:
  uint64_t size_{};
  uint32_t k_{};

  // main memory
  using ull_t = RankAndSelect4Bits::ull_t;
  std::vector<ull_t> w_;
  std::vector<ull_t> last_;
  std::vector<ull_t> is_tip_;
  std::vector<ull_t> invalid_;
  std::vector<uint32_t> tip_node_seq_;
  std::vector<uint64_t> prefix_lkt_;
  std::vector<small_mul_type> small_mul_;
  std::vector<mul_type> large_mul_;
  AtomicBitVector in_or_out_zero;
#ifdef USE_KASH
  khash_t(k64v16) *kh_;
#else
  spp::sparse_hash_map<uint64_t, multi_t> large_mul_lookup_;
#endif

  ull_t f_[kAlphabetSize + 2]{};
  ull_t rank_f_[kAlphabetSize + 2]{}; // = rs_last_.Rank(f_[i] - 1)

  uint64_t num_tip_nodes_{};
  uint32_t uint32_per_tip_nodes_{};
  uint32_t prefix_lk_len_{};

  // auxiliary memory
  RankAndSelect4Bits rs_w_;
  RankAndSelect1Bit rs_last_;
  Rank1Bit rs_is_tip_;
};

using SuccinctDBG = SDBG<>;

#endif // SUCCINCT_DBG_H_