/*
 *  MEGAHIT
 *  Copyright (C) 2014 The University of Hong Kong
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

#ifndef UNITIG_GRAPH_H_
#define UNITIG_GRAPH_H_

#include <vector>
#include <map>
#include <limits>
#include <assert.h>

#include "hash_map.h"

class SuccinctDBG;
struct UnitigGraphVertex {
    UnitigGraphVertex(int64_t start_node, int64_t end_node,
                      int64_t rev_start_node, int64_t rev_end_node, int64_t depth, uint32_t length):
        start_node(start_node), end_node(end_node), rev_start_node(rev_start_node),
        rev_end_node(rev_end_node), depth(depth), length(length) {
        is_deleted = false;
        is_changed = false;
        is_dead = false;
        is_loop = false;
        is_palindrome = false;
    }

    int64_t start_node, end_node;
    int64_t rev_start_node, rev_end_node;
    int64_t depth: 60; // if is_loop, depth is equal to average depth
    bool is_deleted: 1;
    bool is_changed: 1;
    bool is_dead: 1;
    bool is_loop: 1;
    uint32_t length: 31;
    bool is_palindrome: 1;
} __attribute__((packed));

class UnitigGraph {
  public:
    typedef uint32_t vertexID_t;

    UnitigGraph(SuccinctDBG *sdbg): sdbg_(sdbg) {}
    ~UnitigGraph() {}

    void InitFromSdBG();
    uint32_t size() {
        return vertices_.size();
    }
    bool RemoveLocalLowDepth(int min_depth, int min_len, int local_width, double local_ratio, int64_t &num_removed);

    // output
    void OutputInitUnitigs(FILE *contig_file, FILE *multi_file, std::map<int64_t, int> &histo);
    void OutputChangedUnitigs(FILE *addi_contig_file, FILE *addi_multi_file, std::map<int64_t, int> &histo);
    void OutputInitUnitigs(FILE *contig_file, FILE *multi_file, FILE *final_contig_file, std::map<int64_t, int> &histo, int min_final_contig_len);
    void OutputFinalUnitigs(FILE *final_contig_file, std::map<int64_t, int> &histo, int min_final_contig_len);

  private:
    // functions
    double LocalDepth_(UnitigGraphVertex &path, int local_width);
    void Refresh_();

  private:
    // data
    static const size_t kMaxNumVertices;// = std::numeric_limits<vertexID_t>::max();
    SuccinctDBG *sdbg_;
    HashMap<int64_t, vertexID_t> start_node_map_;
    std::vector<UnitigGraphVertex> vertices_;
};

#endif // UNITIG_GRAPH_H_