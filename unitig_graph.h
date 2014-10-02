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
#include "compact_sequence.h"

class SuccinctDBG;
struct UnitigGraphVertex {
    UnitigGraphVertex(int64_t start_node, int64_t end_node, 
        int64_t rev_start_node, int64_t rev_end_node, int64_t depth, const CompactSequence &label): 
            start_node(start_node), end_node(end_node), rev_start_node(rev_start_node), rev_end_node(rev_end_node), depth(depth), label(label) {

        is_deleted = false;
        is_changed = false;
        is_dead = false;
        is_loop = false;
    }

    int64_t start_node, end_node;
    int64_t rev_start_node, rev_end_node;
    int64_t depth: 60; // if is_loop, depth is equal to average depth
    bool is_deleted: 1;
    bool is_changed: 1;
    bool is_dead: 1;
    bool is_loop: 1;
    CompactSequence label;
};

class UnitigGraph {
public:
    UnitigGraph(SuccinctDBG *sdbg): sdbg_(sdbg) {}
    ~UnitigGraph() {}

    void InitFromSdBG();
    uint32_t size() { return vertices_.size(); }
    bool RemoveLocalLowDepth(int min_depth, int min_len, int local_width, double local_ratio, int64_t &num_removed);

    // output without final file, for test only
    void OutputInitUnitigs(FILE *contig_file, FILE *multi_file, std::map<int64_t, int> &histo);
    void OutputChangedUnitigs(FILE *addi_contig_file, FILE *addi_multi_file, std::map<int64_t, int> &histo);
    void OutputFinalUnitigs(FILE *contig_file, std::map<int64_t, int> &histo);

    // output with final file
    void OutputInitUnitigs(FILE *contig_file, FILE *multi_file, FILE *final_contig_file, std::map<int64_t, int> &histo, int min_final_contig_len);
    void OutputFinalUnitigs(FILE *final_contig_file, std::map<int64_t, int> &histo, int min_final_contig_len);

private:
    // functions
    double LocalDepth_(UnitigGraphVertex &path, int local_width);
    void Refresh_();

private:
    // data
    static const size_t kMaxNumVertices = uint32_t(4294967295ULL); // std::numeric_limits<uint32_t>::max();
    SuccinctDBG *sdbg_;
    HashMap<int64_t, uint32_t> start_node_map_;
    std::vector<UnitigGraphVertex> vertices_;
};

#endif // UNITIG_GRAPH_H_