//
// Created by vout on 11/10/18.
//

#ifndef MEGAHIT_UNITIG_GRAPH_H
#define MEGAHIT_UNITIG_GRAPH_H

#include <limits>
#include <atomic>
#include <kmlib/kmbitvector.h>
#include "sparsepp/sparsepp/spp.h"

class SuccinctDBG;

class UnitigGraph {
 public:
  class Vertex {
    struct Status {
      int64_t depth: 60;
      bool is_marked: 1;
      bool is_dead: 1;
      bool is_loop: 1;
    };
    class Adapter {
    };
   private:
    std::atomic<Status> status_;
  };
 public:
  UnitigGraph() = default;

 private:
  const SuccinctDBG *sdbg_{};
  spp::sparse_hash_map<uint64_t, Vertex> start_node_map_;
  kmlib::AtomicBitVector<> locks_;
};

#endif //MEGAHIT_UNITIG_GRAPH_H
