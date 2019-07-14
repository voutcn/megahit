//
// Created by vout on 11/10/18.
//

#ifndef MEGAHIT_UNITIG_GRAPH_H
#define MEGAHIT_UNITIG_GRAPH_H

#include <deque>
#include <limits>
#include "parallel_hashmap/phmap.h"
#include "sdbg/sdbg.h"
#include "unitig_graph_vertex.h"

class UnitigGraph {
 public:
  using Vertex = UnitigGraphVertex;
  using VertexAdapter = UnitigGraphVertex::Adapter;
  using size_type = VertexAdapter::size_type;
  static const size_type kMaxNumVertices =
      std::numeric_limits<size_type>::max() - 1;
  static const size_type kNullVertexID = kMaxNumVertices + 1;

 public:
  explicit UnitigGraph(SDBG *sdbg);
  UnitigGraph(const UnitigGraph &) = delete;
  UnitigGraph(const UnitigGraph &&) = delete;
  ~UnitigGraph() = default;
  size_type size() const { return vertices_.size(); }
  size_t k() const { return sdbg_->k(); }

 public:
  void Refresh(bool mark_changed = false);
  std::string VertexToDNAString(VertexAdapter adapter);

 public:
  /*
   * Function for VertexAdapter obtaining & traversal
   */
  VertexAdapter MakeVertexAdapter(size_type id, int strand = 0) {
    return adapter_impl_.MakeVertexAdapter(id, strand);
  }
  int GetNextAdapters(VertexAdapter &adapter, VertexAdapter *out) {
    return adapter_impl_.GetNextAdapters(adapter, out);
  }
  int GetPrevAdapters(VertexAdapter &adapter, VertexAdapter *out) {
    return adapter_impl_.GetPrevAdapters(adapter, out);
  }
  int OutDegree(VertexAdapter &adapter) {
    return adapter_impl_.OutDegree(adapter);
  }
  int InDegree(VertexAdapter &adapter) {
    return adapter_impl_.InDegree(adapter);
  }

 private:
  /*
   * Function for SudoVertexAdapter obtaining & traversal
   */
  using SudoVertexAdapter = UnitigGraphVertex::SudoAdapter;
  SudoVertexAdapter MakeSudoAdapter(size_type id, int strand = 0) {
    return sudo_adapter_impl_.MakeVertexAdapter(id, strand);
  }
  int GetNextAdapters(SudoVertexAdapter &adapter, SudoVertexAdapter *out) {
    return sudo_adapter_impl_.GetNextAdapters(adapter, out);
  }
  int GetPrevAdapters(SudoVertexAdapter &adapter, SudoVertexAdapter *out) {
    return sudo_adapter_impl_.GetPrevAdapters(adapter, out);
  }
  int OutDegree(SudoVertexAdapter &adapter) {
    return sudo_adapter_impl_.OutDegree(adapter);
  }
  int InDegree(SudoVertexAdapter &adapter) {
    return sudo_adapter_impl_.InDegree(adapter);
  }
  SudoVertexAdapter NextSimplePathAdapter(SudoVertexAdapter &adapter) {
    return sudo_adapter_impl_.NextSimplePathAdapter(adapter);
  }
  SudoVertexAdapter PrevSimplePathAdapter(SudoVertexAdapter &adapter) {
    return sudo_adapter_impl_.PrevSimplePathAdapter(adapter);
  }

 private:
  /**
   * A wrapper for operating different types of adapters
   * @tparam AdapterType type of the vertex adapter
   */
  template <class AdapterType>
  class AdapterImpl {
   public:
    AdapterImpl(UnitigGraph *graph) : graph_(graph) {}

   public:
    AdapterType MakeVertexAdapter(size_type id, int strand = 0) {
      return {graph_->vertices_[id], strand, id};
    }
    int GetNextAdapters(AdapterType &adapter, AdapterType *out) {
      uint64_t next_starts[4];
      int degree = graph_->sdbg_->OutgoingEdges(adapter.e(), next_starts);
      if (out) {
        for (int i = 0; i < degree; ++i) {
          out[i] = MakeVertexAdapterWithSdbgId(next_starts[i]);
        }
      }
      return degree;
    }
    int GetPrevAdapters(AdapterType &adapter, AdapterType *out) {
      adapter.ReverseComplement();
      int degree = GetNextAdapters(adapter, out);
      if (out) {
        for (int i = 0; i < degree; ++i) {
          out[i].ReverseComplement();
        }
      }
      adapter.ReverseComplement();
      return degree;
    }
    int OutDegree(AdapterType &adapter) {
      return GetNextAdapters(adapter, nullptr);
    }
    int InDegree(AdapterType &adapter) {
      adapter.ReverseComplement();
      int degree = OutDegree(adapter);
      adapter.ReverseComplement();
      return degree;
    }
    AdapterType NextSimplePathAdapter(AdapterType &adapter) {
      uint64_t next_sdbg_id = graph_->sdbg_->NextSimplePathEdge(adapter.e());
      if (next_sdbg_id != SDBG::kNullID) {
        return MakeVertexAdapterWithSdbgId(next_sdbg_id);
      } else {
        return AdapterType{};
      }
    }
    AdapterType PrevSimplePathAdapter(AdapterType &adapter) {
      adapter.ReverseComplement();
      AdapterType ret = NextSimplePathAdapter(adapter);
      ret.ReverseComplement();
      adapter.ReverseComplement();
      return ret;
    }

   private:
    AdapterType MakeVertexAdapterWithSdbgId(uint64_t sdbg_id) {
      uint32_t id = graph_->id_map_.at(sdbg_id);
      AdapterType adapter(graph_->vertices_[id], 0, id);
      if (adapter.b() != sdbg_id) {
        adapter.ReverseComplement();
      }
      return adapter;
    }

   private:
    UnitigGraph *graph_;
  };

  void RefreshDisconnected();

 private:
  SDBG *sdbg_{};
  std::deque<UnitigGraphVertex> vertices_;
  phmap::flat_hash_map<uint64_t, size_type> id_map_;
  AdapterImpl<VertexAdapter> adapter_impl_;
  AdapterImpl<SudoVertexAdapter> sudo_adapter_impl_;
};

#endif  // MEGAHIT_UNITIG_GRAPH_H
