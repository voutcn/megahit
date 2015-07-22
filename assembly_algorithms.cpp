/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong
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

#include "assembly_algorithms.h"
#include <assert.h>
#include <stdio.h>
#include <omp.h>
#include <assert.h>
#include <vector>
#include <map>
#include <algorithm>
#include <unordered_set>
#include <parallel/algorithm>

#include "atomic_bit_vector.h"
#include "unitig_graph.h"
#include "utils.h"

using std::vector;
using std::string;
using std::map;

namespace assembly_algorithms {

static AtomicBitVector marked;

int64_t Trim(SuccinctDBG &dbg, int len, int min_final_standalone) {
    int64_t number_tips = 0;
    marked.reset(dbg.size);

    #pragma omp parallel for reduction(+:number_tips)
    for (int64_t edge_idx = 0; edge_idx < dbg.size; ++edge_idx) {
        if (dbg.IsValidEdge(edge_idx) && !marked.get(edge_idx) && dbg.EdgeOutdegree(edge_idx) == 0) {
            vector<int64_t> path = {edge_idx};
            int64_t prev_edge;
            int64_t cur_edge = edge_idx;
            bool is_tip = false;
            for (int i = 1; i < len; ++i) {
                prev_edge = dbg.UniquePrevEdge(cur_edge);
                if (prev_edge == -1) {
                    is_tip = dbg.EdgeIndegree(cur_edge) == 0; // && (i + dbg.kmer_k - 1 < min_final_standalone);
                    break;
                } else if (dbg.UniqueNextEdge(prev_edge) == -1) {
                    is_tip = true;
                    break;
                } else {
                    path.push_back(prev_edge);
                    cur_edge = prev_edge;
                }
            }

            if (is_tip) {
                for (unsigned i = 0; i < path.size(); ++i) {
                    marked.set(path[i]);
                }
                ++number_tips;
            }
        }
    }

    #pragma omp parallel for reduction(+:number_tips)
    for (int64_t edge_idx = 0; edge_idx < dbg.size; ++edge_idx) {
        if (dbg.IsValidEdge(edge_idx) && !marked.get(edge_idx) && dbg.EdgeIndegree(edge_idx) == 0) {
            vector<int64_t> path = {edge_idx};
            int64_t next_node;
            int64_t cur_edge = edge_idx;
            bool is_tip = false;
            for (int i = 1; i < len; ++i) {
                next_node = dbg.UniqueNextEdge(cur_edge);
                if (next_node == -1) {
                    is_tip = dbg.EdgeOutdegree(cur_edge) == 0; // && (i + dbg.kmer_k - 1 < min_final_standalone);
                    break;
                } else if (dbg.UniquePrevEdge(next_node) == -1) {
                    is_tip = true;
                } else {
                    path.push_back(next_node);
                    cur_edge = next_node;
                }
            }

            if (is_tip) {
                for (unsigned i = 0; i < path.size(); ++i) {
                    marked.set(path[i]);
                }
                ++number_tips;
            }
        }
    }

    #pragma omp parallel for
    for (int64_t edge_idx = 0; edge_idx < dbg.size; ++edge_idx) {
        if (marked.get(edge_idx)) {
            dbg.SetInvalidEdge(edge_idx);
        }
    }

    return number_tips;
}

int64_t RemoveTips(SuccinctDBG &dbg, int max_tip_len, int min_final_standalone) {
    int64_t number_tips = 0;
    xtimer_t timer;
    for (int len = 2; len < max_tip_len; len *= 2) {
        xlog("Removing tips with length less than %d; ", len);
        timer.reset();
        timer.start();
        number_tips += Trim(dbg, len, min_final_standalone);
        timer.stop();
        xlog_ext("Accumulated tips removed: %lld; time elapsed: %.4f\n", (long long)number_tips, timer.elapsed());
    }
    xlog("Removing tips with length less than %d; ", max_tip_len);
    timer.reset();
    timer.start();
    number_tips += Trim(dbg, max_tip_len, min_final_standalone);
    timer.stop();
    xlog_ext("Accumulated tips removed: %lld; time elapsed: %.4f\n", (long long)number_tips, timer.elapsed());

    {
        AtomicBitVector empty;
        marked.swap(empty);
    }

    return number_tips;
}

} // namespace assembly_algorithms