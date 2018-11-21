//
// Created by vout on 11/20/18.
//

#ifndef MEGAHIT_UNITIG_GRAPH_VERTEX_H
#define MEGAHIT_UNITIG_GRAPH_VERTEX_H

#include <cstdint>
#include <atomic>
#include <algorithm>

/**
 * store the metadata of a unitig vertex; the vertex is associate with an SDBG
 */
class UnitigGraphVertex {
 public:
  static const uint8_t kUnknownDegree = 7;
  UnitigGraphVertex() = default;
  UnitigGraphVertex(uint64_t begin, uint64_t end, uint64_t rbegin, uint64_t rend,
                    uint64_t total_depth, uint32_t length, bool is_loop = false)
      : strand_info{{begin, end}, {rbegin, rend}},
        length(length), total_depth(total_depth), is_looped(is_loop), is_palindrome(begin == rbegin),
        is_changed(false), is_deleted(false), to_delete(false), flag(0) {}
 private:
  struct StrandInfo {
    StrandInfo(uint64_t begin = 0, uint64_t end = 0) :
        begin(begin), end(end), cached_out_degree(kUnknownDegree), to_disconnect(false) {}
    uint64_t begin, end: 60;
    uint8_t cached_out_degree: 3;
    bool to_disconnect: 1;
  };
  StrandInfo strand_info[2];
  uint32_t length{};
  uint64_t total_depth : 51;
  bool is_looped: 1;
  bool is_palindrome: 1;
  // status
  bool is_changed: 1;
  bool is_deleted: 1;
  bool to_delete: 1;
  uint8_t flag{0};

 public:
  /**
   * An adapter on unitig vertex to handle graph traversal, modification
   * and reverse complement issues
   */
  class Adapter {
   public:
    Adapter() = default;
    Adapter(UnitigGraphVertex &vertex, int strand = 0, uint32_t id = static_cast<uint32_t>(-1))
        : vertex_(&vertex), strand_(strand), id_(id) {}
    void ReverseComplement() { strand_ ^= 1; }
    uint32_t id() const { return id_; }
    int strand() const { return strand_; }
    bool valid() const { return vertex_ != nullptr; }
    uint32_t length() const { return vertex_->length; }
    uint64_t total_depth() const { return vertex_->total_depth; }
    double avg_depth() const { return static_cast<double>(vertex_->total_depth) / vertex_->length; }
    bool is_loop() const { return vertex_->is_looped; }
    bool is_palindrome() const { return vertex_->is_palindrome; }
    bool is_to_delete() const { return vertex_->to_delete; }
    bool is_changed() const { return vertex_->is_changed; }
    bool is_to_disconnect() const { return strand_info().to_disconnect; }
    uint8_t cached_out_degree() const { return strand_info().cached_out_degree; }
    uint64_t rep_id() const { return std::min(begin(), rbegin()); };

    void MarkToDelete() { vertex_->to_delete = true; }
    void MarkToDisconnect() { strand_info().to_disconnect = true; }

    uint64_t begin() const { return strand_info().begin; }
    uint64_t end() const { return strand_info().end; }
    uint64_t rbegin() const { return strand_info(1).begin; }
    uint64_t rend() const { return strand_info(1).end; }

   protected:
    UnitigGraphVertex::StrandInfo &strand_info(int relative_strand = 0) {
      return vertex_->strand_info[strand_ ^ relative_strand];
    }
    const UnitigGraphVertex::StrandInfo &strand_info(int relative_strand = 0) const {
      return vertex_->strand_info[strand_ ^ relative_strand];
    }

   protected:
    UnitigGraphVertex *vertex_{nullptr};
    uint8_t strand_: 1;
    uint32_t id_: 31;  // record the id for quick access
  };

  /**
   * An adapter with more permission to change the data of the vertex
   */
  class SudoAdapter : public Adapter {
   public:
    SudoAdapter() = default;
    SudoAdapter(UnitigGraphVertex &vertex, int strand = 0, uint32_t id = static_cast<uint32_t>(-1))
        : Adapter(vertex, strand, id) {}
   public:
    void set_start_end(uint64_t start, uint64_t end, uint64_t rc_start, uint64_t rc_end) {
      strand_info(0) = {start, end};
      strand_info(1) = {rc_start, rc_end};
      vertex_->is_palindrome = start == rc_start;
    };
    void set_length(uint32_t length) {
      vertex_->length = length;
    }
    void set_total_depth(uint64_t total_depth) {
      vertex_->total_depth = total_depth;
    }
    void set_cached_outdegree(uint8_t degree) {
      strand_info().cached_out_degree = degree;
    }
    void set_looped() {
      vertex_->is_looped = true;
    }
    void clear_cache() {
      strand_info(0).cached_out_degree = strand_info(1).cached_out_degree = kUnknownDegree;
    }
    uint8_t flag() const { return vertex_->flag; }
    void set_flag(uint8_t flag) { vertex_->flag = flag; }
    void set_changed() { vertex_->is_changed = true; }
  };
};

//static_assert(sizeof(UnitigGraphVertex) <= 44, "");
//static_assert(sizeof(UnitigGraphVertex::Adapter) <= sizeof(UnitigGraphVertex*) + sizeof(uint32_t), "");

#endif //MEGAHIT_UNITIG_GRAPH_VERTEX_H
