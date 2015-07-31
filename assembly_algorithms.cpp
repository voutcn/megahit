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

static AtomicBitVector removed_nodes;

int64_t Trim(SuccinctDBG &dbg, int len, int min_final_standalone) {
    int64_t number_tips = 0;

    #pragma omp parallel for reduction(+:number_tips)
    for (int64_t node_idx = 0; node_idx < dbg.size; ++node_idx) {
        if (dbg.IsLast(node_idx) && !removed_nodes.get(node_idx) && dbg.NodeOutdegreeZero(node_idx)) {
            vector<int64_t> path = {node_idx};
            int64_t prev_node;
            int64_t cur_node = node_idx;
            bool is_tip = false;
            for (int i = 1; i < len; ++i) {
                prev_node = dbg.UniquePrevNode(cur_node);
                if (prev_node == -1) {
                    is_tip = dbg.NodeIndegreeZero(cur_node); // && (i + dbg.kmer_k - 1 < min_final_standalone);
                    break;
                } else if (dbg.UniqueNextNode(prev_node) == -1) {
                    is_tip = true;
                    break;
                } else {
                    path.push_back(prev_node);
                    cur_node = prev_node;
                }
            }

            if (is_tip) {
                for (unsigned i = 0; i < path.size(); ++i) {
                    removed_nodes.set(path[i]);
                }
                ++number_tips;
            }
        }
    }

    #pragma omp parallel for reduction(+:number_tips)
    for (int64_t node_idx = 0; node_idx < dbg.size; ++node_idx) {
        if (dbg.IsLast(node_idx) && !removed_nodes.get(node_idx) && dbg.NodeIndegreeZero(node_idx)) {
            vector<int64_t> path = {node_idx};
            int64_t next_node;
            int64_t cur_node = node_idx;
            bool is_tip = false;
            for (int i = 1; i < len; ++i) {
                next_node = dbg.UniqueNextNode(cur_node);
                if (next_node == -1) {
                    is_tip = dbg.NodeOutdegreeZero(cur_node); // && (i + dbg.kmer_k - 1 < min_final_standalone);
                    break;
                } else if (dbg.UniquePrevNode(next_node) == -1) {
                    is_tip = true;
                } else {
                    path.push_back(next_node);
                    cur_node = next_node;
                }
            }

            if (is_tip) {
                for (unsigned i = 0; i < path.size(); ++i) {
                    removed_nodes.set(path[i]);
                }
                ++number_tips;
            }
        }
    }

    #pragma omp parallel for
    for (int64_t node_idx = 0; node_idx < dbg.size; ++node_idx) {
        if (removed_nodes.get(node_idx)) {
            dbg.DeleteAllEdges(node_idx);
        }
    }

    return number_tips;
}

int64_t RemoveTips(SuccinctDBG &dbg, int max_tip_len, int min_final_standalone) {
    int64_t number_tips = 0;
    xtimer_t timer;
    removed_nodes.reset(dbg.size);
    
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
        removed_nodes.swap(empty);
    }

    return number_tips;
}

} // namespace assembly_algorithms