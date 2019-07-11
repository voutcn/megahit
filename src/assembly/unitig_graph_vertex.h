//
// Created by vout on 11/20/18.
//

#ifndef MEGAHIT_UNITIG_GRAPH_VERTEX_H
#define MEGAHIT_UNITIG_GRAPH_VERTEX_H

#include <utils/atomic_wrapper.h>
#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstdint>

/**
 * store the metadata of a unitig vertex; the vertex is associate with an SDBG
 */
class UnitigGraphVertex {
 public:
  UnitigGraphVertex() = default;
  UnitigGraphVertex(uint64_t begin, uint64_t end, uint64_t rbegin,
                    uint64_t rend, uint64_t total_depth, uint32_t length,
                    bool is_loop = false)
      : strand_info{{begin, end}, {rbegin, rend}},
        length(length),
        total_depth(total_depth),
        is_looped(is_loop),
        is_palindrome(begin == rbegin),
        is_changed(false),
        flag(0) {}

 private:
  struct StrandInfo {
    StrandInfo(uint64_t begin = 0, uint64_t end = 0) : begin(begin), end(end) {}
    uint64_t begin : 48;
    uint64_t end : 48;
  } __attribute__((packed));
  StrandInfo strand_info[2];
  uint32_t length;
  uint64_t total_depth : 52;
  bool is_looped : 1;
  bool is_palindrome : 1;
  bool is_changed : 1;
  // status that can be modified by adapter during traversal and must be atomic
  AtomicWrapper<uint8_t> flag;  // bit 0-4: any flag; bit 5: marked as to
                                // delete; bit 6 & 7: marked as to disconnect
  static const unsigned kToDeleteBit = 5;
  static const unsigned kToDisconnectBit = 6;
  static const uint8_t kFlagMask = (1u << kToDeleteBit) - 1;

 public:
  /**
   * An adapter on unitig vertex to handle graph traversal, modification
   * and reverse complement issues
   */
  class Adapter {
   public:
    using size_type = uint32_t;
    Adapter() = default;
    Adapter(UnitigGraphVertex &vertex, int strand = 0,
            size_type id = static_cast<size_type>(-1))
        : vertex_(&vertex), strand_(strand), id_(id) {}
    void ReverseComplement() { strand_ ^= 1u; }
    size_type UnitigId() const { return id_; }
    bool IsValid() const { return vertex_ != nullptr; }
    uint32_t GetLength() const { return vertex_->length; }
    uint64_t GetTotalDepth() const { return vertex_->total_depth; }
    double GetAvgDepth() const {
      return static_cast<double>(vertex_->total_depth) / vertex_->length;
    }
    bool IsLoop() const { return vertex_->is_looped; }
    bool IsPalindrome() const { return vertex_->is_palindrome; }
    bool IsChanged() const { return vertex_->is_changed; }
    void ToUniqueFormat() {
      if (canonical_id() != b()) {
        ReverseComplement();
      }
    }
    bool IsStandalone() const { return IsLoop(); }
    bool SetToDelete() {
      uint8_t mask = 1u << kToDeleteBit;
      auto old_val = vertex_->flag.v.fetch_or(mask, std::memory_order_relaxed);
      return !(old_val & mask);
    }
    bool SetToDisconnect() {
      uint8_t mask = 1u << kToDisconnectBit << strand_;
      auto old_val = vertex_->flag.v.fetch_or(mask, std::memory_order_relaxed);
      return !(old_val & mask);
    }

    uint64_t canonical_id() const { return std::min(b(), rb()); };
    uint64_t b() const { return StrandInfo().begin; }
    uint64_t e() const { return StrandInfo().end; }
    uint64_t rb() const { return StrandInfo(1).begin; }
    uint64_t re() const { return StrandInfo(1).end; }

   protected:
    UnitigGraphVertex::StrandInfo &StrandInfo(uint8_t relative_strand = 0) {
      return vertex_->strand_info[strand_ ^ relative_strand];
    }
    const UnitigGraphVertex::StrandInfo &StrandInfo(
        uint8_t relative_strand = 0) const {
      return vertex_->strand_info[strand_ ^ relative_strand];
    }

   protected:
    UnitigGraphVertex *vertex_{nullptr};
    uint8_t strand_;
    uint32_t id_;  // record the id for quick access
  };

  /**
   * An adapter with more permission to change the data of the vertex
   */
  class SudoAdapter : public Adapter {
   public:
    SudoAdapter() = default;
    SudoAdapter(UnitigGraphVertex &vertex, int strand = 0,
                size_type id = static_cast<size_type>(-1))
        : Adapter(vertex, strand, id) {}

   public:
    bool IsToDelete() const {
      return vertex_->flag.v.load(std::memory_order_relaxed) &
             (1u << kToDeleteBit);
    }
    bool IsToDisconnect() const {
      return vertex_->flag.v.load(std::memory_order_relaxed) &
             (1u << kToDisconnectBit << strand_);
    }
    void SetBeginEnd(uint64_t start, uint64_t end, uint64_t rc_start,
                     uint64_t rc_end) {
      StrandInfo(0) = {start, end};
      StrandInfo(1) = {rc_start, rc_end};
      vertex_->is_palindrome = start == rc_start;
    };
    void SetLength(uint32_t length) { vertex_->length = length; }
    void SetTotalDepth(uint64_t total_depth) {
      vertex_->total_depth = total_depth;
    }
    void SetLooped() { vertex_->is_looped = true; }
    void SetChanged() { vertex_->is_changed = true; }
    uint8_t GetFlag() const {
      return vertex_->flag.v.load(std::memory_order_relaxed) & kFlagMask;
    }
    void SetFlag(uint8_t flag) {
      assert(flag <= kFlagMask);
      vertex_->flag.v.store(flag, std::memory_order_relaxed);
    }
  };
};

static_assert(sizeof(UnitigGraphVertex) <= 40, "");

#endif  // MEGAHIT_UNITIG_GRAPH_VERTEX_H
