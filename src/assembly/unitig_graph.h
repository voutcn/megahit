//
// Created by vout on 11/10/18.
//

#ifndef MEGAHIT_UNITIG_GRAPH_H
#define MEGAHIT_UNITIG_GRAPH_H

#include <limits>
#include <deque>
#include <kmlib/kmbitvector.h>
#include "sdbg/sdbg.h"
#include "sparsepp/sparsepp/spp.h"

namespace assembly {

class UnitigGraph {
 public:
  using size_type = uint32_t;
  static const size_type kMaxNumVertices = std::numeric_limits<size_type>::max() - 1;
  static const size_type kNullVertexID = kMaxNumVertices + 1;

  /**
   * store the metadata of a unitig vertex; the vertex is associate with an SDBG
   */
  struct Vertex {
    Vertex() = default;
    Vertex(uint64_t start, uint64_t end, uint64_t rc_start, uint64_t rc_end,
           uint64_t total_depth, uint32_t length, bool is_loop = false)
        : strand_info{{start, end, kUnknownDegree, false}, {rc_start, rc_end, kUnknownDegree, false}},
          length(length), total_depth(total_depth), is_loop(is_loop), is_palindrome(start == rc_start),
          is_changed(false), is_deleted(false), to_delete(false), flag(0) {
      if (start < rc_start) {
        std::swap(strand_info[0], strand_info[1]);
      }
    }
    struct StrandInfo {
      uint64_t start, end: 60;
      uint8_t cached_out_degree: 3;
      bool to_disconnect: 1;
    };
    StrandInfo strand_info[2]{};
    uint32_t length{};
    uint64_t total_depth : 51;
    bool is_loop: 1;
    bool is_palindrome: 1;
    // status
    bool is_changed: 1;
    bool is_deleted: 1;
    bool to_delete: 1;
    uint8_t flag{0};
    static const uint8_t kUnknownDegree = 7;
  };

  /**
   * An adapter on unitig vertex to handle graph traversal, modification
   * and reverse complement issues
   */
  template<class AdapterType>
  class VertexAdapterBase {
   public:
    VertexAdapterBase() = default;
    VertexAdapterBase(UnitigGraph &graph, Vertex &vertex, size_type id, int strand)
        : graph_(&graph), vertex_(&vertex), id_(id), strand_(strand) {}
    void ReverseComplement() { strand_ ^= 1; }
    const Vertex &vertex() const { return *vertex_; }
    size_type id() const { return id_; }
    bool valid() const { return vertex_ != nullptr; }
    uint32_t length() const { return vertex_->length; }
    uint64_t total_depth() const { return vertex_->total_depth; }
    double depth() const { return static_cast<double>(vertex_->total_depth) / vertex_->length; }
    bool is_loop() const { return vertex_->is_loop; }
    bool is_palindrome() const { return vertex_->is_palindrome; }
    bool is_to_delete() const { return vertex_->to_delete; }
    bool is_deleted() const { return vertex_->is_deleted; }
    bool is_to_disconnected() const { return strand_info().to_disconnect; }
    uint64_t rep_id() const {
      return std::min(vertex_->strand_info[0].start, vertex_->strand_info[1].start);
    };

    void MarkToDelete() { vertex_->to_delete = true; }
    void MarkToDisconnect() { strand_info().to_disconnect = true; }

    uint64_t start() const { return vertex_->strand_info[strand_].start; }
    uint64_t end() const { return vertex_->strand_info[strand_].end; }
    uint64_t rc_start() const { return vertex_->strand_info[strand_ ^ 1].start; }
    uint64_t rc_end() const { return vertex_->strand_info[strand_ ^ 1].end; }

    int OutDegree() {
      if (strand_info().cached_out_degree != Vertex::kUnknownDegree) {
        return strand_info().cached_out_degree;
      } else {
        return GetNextAdapters(nullptr);
      }
    }
    int InDegree() {
      ReverseComplement();
      int degree = OutDegree();
      ReverseComplement();
      return degree;
    }
    int GetNextAdapters(AdapterType *out) {
      uint64_t next_starts[4];
      int degree = graph_->sdbg_->OutgoingEdges(end(), next_starts);
      if (out) {
        for (int i = 0; i < degree; ++i) {
          auto id = graph_->start_node_map_.at(next_starts[i]);
          auto &vertex = graph_->vertices_[id];
          out[i] = AdapterType(*graph_, vertex, id, vertex.strand_info[0].start != next_starts[i]);
        }
      }
      return strand_info().cached_out_degree = degree;
    }
    int GetPrevAdapters(AdapterType *out) {
      ReverseComplement();
      int degree = GetNextAdapters(out);
      if (out) {
        for (int i = 0; i < degree; ++i) {
          out[i].ReverseComplement();
        }
      }
      ReverseComplement();
      return degree;
    }
    AdapterType NextSimplePathAdapter() {
      uint64_t next_sdbg_id = graph_->sdbg_->NextSimplePathEdge(end());
      if (next_sdbg_id != SDBG::kNullID) {
        size_type id = graph_->start_node_map_.at(next_sdbg_id);
        Vertex &vertex = graph_->vertices_[id];
        return AdapterType{*graph_, vertex, id, vertex.strand_info[0].start != next_sdbg_id};
      } else {
        return AdapterType{};
      }
    }
    AdapterType PrevSimplePathAdapter() {
      ReverseComplement();
      auto adapter = NextSimplePathAdapter();
      ReverseComplement();
      adapter.ReverseComplement();
      return adapter;
    }

   protected:
    Vertex::StrandInfo &strand_info() {
      return vertex_->strand_info[strand_];
    }
    const Vertex::StrandInfo &strand_info() const {
      return vertex_->strand_info[strand_];
    }

   protected:
    UnitigGraph *graph_{nullptr};
    Vertex *vertex_{nullptr};
    size_type id_{};
    int strand_{};
  };

  /**
   * Basic adapter to be used outside the UnitigGraph class
   */
  class VertexAdapter : public VertexAdapterBase<VertexAdapter> {
   public:
    VertexAdapter() = default;
    VertexAdapter(UnitigGraph &graph, Vertex &vertex, size_type id, int strand)
        : VertexAdapterBase<VertexAdapter>(graph, vertex, id, strand) {}
  };

 private:
  /**
   * A VertexAdapter with more permission to change the data of the vertex
   */
  class SudoVertexAdapter : public VertexAdapterBase<SudoVertexAdapter> {
   public:
    SudoVertexAdapter() = default;
    SudoVertexAdapter(UnitigGraph &graph, Vertex &vertex, size_type id, int strand) :
        VertexAdapterBase<SudoVertexAdapter>(graph, vertex, id, strand) {}
    void MarkDeleted() { vertex_->is_deleted = true; }
    void SetFlag(uint8_t flag) { vertex_->flag = flag; }
  };

 public:
  explicit UnitigGraph(SuccinctDBG *sdbg);
  UnitigGraph(const UnitigGraph &) = delete;
  UnitigGraph(const UnitigGraph &&) = delete;
  ~UnitigGraph() = default;

  size_type size() const { return vertices_.size(); }
  VertexAdapter MakeVertexAdapter(size_type id, int strand = 0) {
    return {*this, vertices_[id], id, strand};
  }
  void Refresh(bool mark_changed = false);

 private:
  void RefreshDisconnected();
  SudoVertexAdapter MakeSudoAdapter(size_type id, int strand = 0) {
    return {*this, vertices_[id], id, strand};
  }
 protected:
  SuccinctDBG *sdbg_{};
  std::deque<Vertex> vertices_;
  spp::sparse_hash_map<uint64_t, size_type> start_node_map_;
  kmlib::AtomicBitVector<> locks_;
};
}

#endif //MEGAHIT_UNITIG_GRAPH_H
