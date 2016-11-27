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

#ifndef UNITIG_GRAPH_H_
#define UNITIG_GRAPH_H_

#include <vector>
#include <limits>
#include <assert.h>
#include <omp.h>

#include "hash_map.h"
#include "histgram.h"

class SuccinctDBG;
struct UnitigGraphVertex {
    UnitigGraphVertex(int64_t start_node, int64_t end_node,
                      int64_t rev_start_node, int64_t rev_end_node, int64_t depth, uint32_t length):
        start_node(start_node), end_node(end_node), rev_start_node(rev_start_node),
        rev_end_node(rev_end_node), depth(depth), length(length) {
        is_deleted = false;
        is_changed = false;
        is_marked = false;
        is_dead = false;
        is_loop = false;
        is_palindrome = false;
    }

    int64_t start_node, end_node;
    int64_t rev_start_node, rev_end_node;
    int64_t depth: 60; // if is_loop, depth is equal to average depth
    bool is_deleted: 1;
    bool is_changed: 1;
    bool is_marked: 1;
    bool is_dead: 1;
    bool is_loop: 1;
    uint32_t length: 30;
    bool is_palindrome: 1;

    UnitigGraphVertex ReverseComplement() {
        UnitigGraphVertex ret = *this;
        std::swap(ret.start_node, ret.rev_start_node);
        std::swap(ret.end_node, ret.rev_end_node);
        return ret;
    }

    int64_t Representation() {
        int64_t ret = std::max(start_node, end_node);
        ret = std::max(rev_start_node, ret);
        ret = std::max(rev_end_node, ret);
        return ret;
    }
};

class UnitigGraph {
  public:
    typedef uint32_t vertexID_t;

    UnitigGraph(SuccinctDBG *sdbg): sdbg_(sdbg) {}
    ~UnitigGraph() {
        for (vertexID_t i = 0; i < locks_.size(); ++i) {
            omp_destroy_lock(&locks_[i]);
        }
    }

    void InitFromSdBG();
    uint32_t size() {
        return vertices_.size();
    }
    int64_t RemoveLowDepth(double min_depth);
    bool RemoveLocalLowDepth(double min_depth, int min_len, int local_width, double local_ratio, int64_t &num_removed, bool permanent_rm = false);
    uint32_t MergeBubbles(bool permanent_rm, bool careful, FILE *bubble_file, Histgram<int64_t> &hist);
    uint32_t MergeComplexBubbles(double similarity, int merge_level, bool permanent_rm, bool careful, FILE *bubble_file, Histgram<int64_t> &hist);
    uint32_t DisconnectWeakLinks(double local_ratio);
    uint32_t RemoveTips(int max_tip_len);
    uint32_t MergeSuperBubbles(int max_len, bool permanent_rm, bool careful, FILE *bubble_file, Histgram<int64_t> &hist);
    int SearchAndMergeSuperBubble_(int64_t source, int max_len, bool careful, FILE *bubble_file, Histgram<int64_t> &hist);

    // output
    void OutputContigs(FILE *contig_file, FILE *final_file, Histgram<int64_t> &hist, bool change_only, int min_final_standalone);

  private:
    // functions
    double LocalDepth_(vertexID_t id, int local_width);
    void Refresh_(bool set_changed = true);

  private:
    // data
    static const size_t kMaxNumVertices;// = std::numeric_limits<vertexID_t>::max();
    SuccinctDBG *sdbg_;
    HashMap<int64_t, vertexID_t> start_node_map_;
    std::vector<UnitigGraphVertex> vertices_;
    std::vector<omp_lock_t> locks_;
};

#endif // UNITIG_GRAPH_H_