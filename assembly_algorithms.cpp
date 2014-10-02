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

#include "branch_group.h"
#include "atomic_bit_vector.h"
#include "unitig_graph.h"
#include "timer.h"

using std::vector;
using std::string;
using std::map;

namespace assembly_algorithms {

static AtomicBitVector marked;
static map<int64_t, int> histogram;
static inline void MarkNode(SuccinctDBG &dbg, int64_t node_idx);

int64_t NextSimplePathNode(SuccinctDBG &dbg, int64_t cur_node) {
    int64_t next_node = dbg.UniqueOutgoing(cur_node);
    if (next_node != -1 && dbg.UniqueIncoming(next_node) != -1) {
        return next_node;
    } else {
        return -1;
    }
}

int64_t PrevSimplePathNode(SuccinctDBG &dbg, int64_t cur_node) {
    int64_t prev_node = dbg.UniqueIncoming(cur_node);
    if (prev_node != -1 && dbg.UniqueOutgoing(prev_node) != -1) {
        return prev_node;
    } else {
        return -1;
    }
}

int64_t Trim(SuccinctDBG &dbg, int len, int min_final_contig_len) {
    int64_t number_tips = 0;
    omp_lock_t path_lock;
    omp_init_lock(&path_lock);
    marked.reset(dbg.size);

#pragma omp parallel for reduction(+:number_tips)  
    for (int64_t node_idx = 0; node_idx < dbg.size; ++node_idx) {
        if (dbg.IsValidNode(node_idx) && !marked.get(node_idx) && dbg.IsLast(node_idx) && dbg.OutdegreeZero(node_idx)) {
            vector<int64_t> path = {node_idx};
            int64_t prev_node;
            int64_t cur_node = node_idx;
            bool is_tip = false;
            for (int i = 1; i < len; ++i) {
                prev_node = dbg.UniqueIncoming(cur_node);
                if (prev_node == -1) {
                    is_tip = dbg.IndegreeZero(cur_node) && (i + dbg.kmer_k - 1 < min_final_contig_len);
                    break;
                } else if (dbg.UniqueOutgoing(prev_node) == -1) {
                    is_tip = true;
                    break;
                } else {
                    path.push_back(prev_node);
                    cur_node = prev_node;
                }
            }

            if (is_tip) {
                for (unsigned i = 0; i < path.size(); ++i) {
                    MarkNode(dbg, path[i]);
                }
                ++number_tips;
            }
        }
    }

#pragma omp parallel for reduction(+:number_tips)
    for (int64_t node_idx = 0; node_idx < dbg.size; ++node_idx) {
        if (dbg.IsValidNode(node_idx) && dbg.IsLast(node_idx) && !marked.get(node_idx) && dbg.IndegreeZero(node_idx)) {
            vector<int64_t> path = {node_idx};
            int64_t next_node;
            int64_t cur_node = node_idx;
            bool is_tip = false;
            for (int i = 1; i < len; ++i) {
                next_node = dbg.UniqueOutgoing(cur_node);
                if (next_node == -1) {
                    is_tip = dbg.OutdegreeZero(cur_node) && (i + dbg.kmer_k - 1 < min_final_contig_len);
                    break;
                } else if (dbg.UniqueIncoming(next_node) == -1) {
                    is_tip = true;
                } else {
                    path.push_back(next_node);
                    cur_node = next_node;
                }
            }

            if (is_tip) {
                for (unsigned i = 0; i < path.size(); ++i) {
                    MarkNode(dbg, path[i]);
                }
                ++number_tips;
            }
        }
    }

#pragma omp parallel for
    for (int64_t node_idx = 0; node_idx < dbg.size; ++node_idx) {
        if (marked.get(node_idx)) {
            dbg.SetInvalid(node_idx);
        }
    }

    return number_tips;
}

int64_t RemoveTips(SuccinctDBG &dbg, int max_tip_len, int min_final_contig_len) {
    int64_t number_tips = 0;
    xtimer_t timer;
    for (int len = 2; len < max_tip_len; len *= 2) {
        printf("Removing tips with length less than %d\n", len);
        timer.reset();
        timer.start();
        number_tips += Trim(dbg, len, min_final_contig_len);
        timer.stop();
        printf("Accumulated tips removed: %ld; time elapsed: %.4f\n", number_tips, timer.elapsed());
    }
    printf("Removing tips with length less than %d\n", max_tip_len);
    timer.reset();
    timer.start();
    number_tips += Trim(dbg, max_tip_len, min_final_contig_len);
    timer.stop();
    printf("Accumulated tips removed: %ld; time elapsed: %.4f\n", number_tips, timer.elapsed());
    return number_tips;
}

int64_t PopBubbles(SuccinctDBG &dbg, int max_bubble_len, double low_depth_ratio) {
    omp_lock_t bubble_lock;
    omp_init_lock(&bubble_lock);
    const int kMaxBranchesPerGroup = 4;
    if (max_bubble_len <= 0) { max_bubble_len = dbg.kmer_k * 2 + 2; }
    vector<std::pair<int, int64_t> > bubble_candidates;
    int64_t num_bubbles = 0;

#pragma omp parallel for
    for (int64_t node_idx = 0; node_idx < dbg.size; ++node_idx) {
        if (dbg.IsValidNode(node_idx) && dbg.IsLast(node_idx) && dbg.Outdegree(node_idx) > 1) {
            BranchGroup bubble(&dbg, node_idx, kMaxBranchesPerGroup, max_bubble_len);
            if (bubble.Search()) {
                omp_set_lock(&bubble_lock);
                bubble_candidates.push_back(std::make_pair(bubble.length(), node_idx));
                omp_unset_lock(&bubble_lock);
            }
        }
    }

    for (unsigned i = 0; i < bubble_candidates.size(); ++i) {
        BranchGroup bubble(&dbg, bubble_candidates[i].second, kMaxBranchesPerGroup, max_bubble_len);
        if (bubble.Search() && bubble.RemoveErrorBranches(low_depth_ratio)) {
            ++num_bubbles;
        }
    }

    omp_destroy_lock(&bubble_lock);
    return num_bubbles;
}

void AssembleFromUnitigGraph(SuccinctDBG &dbg, FILE *contigs_file, FILE *multi_file, FILE *final_contig_file, int min_final_contig_len) {
    xtimer_t timer;
    timer.reset();
    timer.start();
    UnitigGraph unitig_graph(&dbg);
    unitig_graph.InitFromSdBG();
    timer.stop();
    printf("unitig graph size: %u, time for building: %lf\n", unitig_graph.size(), timer.elapsed());
    
    timer.reset();
    timer.start();
    histogram.clear();
    if (final_contig_file == NULL) {
        unitig_graph.OutputInitUnitigs(contigs_file, multi_file, histogram);
    } else {
        unitig_graph.OutputInitUnitigs(contigs_file, multi_file, final_contig_file, histogram, min_final_contig_len);
    }
    PrintStat();
    timer.stop();
    printf("Time to output: %lf\n", timer.elapsed());
}

void AssembleFinalFromUnitigGraph(SuccinctDBG &dbg, FILE *final_contig_file, int min_final_contig_len) {
    xtimer_t timer;
    timer.reset();
    timer.start();
    UnitigGraph unitig_graph(&dbg);
    unitig_graph.InitFromSdBG();
    timer.stop();
    printf("unitig graph size: %u, time for building: %lf\n", unitig_graph.size(), timer.elapsed());
    
    timer.reset();
    timer.start();
    histogram.clear();
    unitig_graph.OutputFinalUnitigs(final_contig_file, histogram, min_final_contig_len);
    PrintStat();
    timer.stop();
    printf("Time to output: %lf\n", timer.elapsed());
}

void RemoveLowLocalAndOutputChanged(SuccinctDBG &dbg, FILE *contigs_file, FILE *multi_file, FILE *final_contig_file, 
                                    FILE *addi_contig_file, FILE *addi_multi_file, 
                                    double min_depth, int min_len, double local_ratio, int min_final_contig_len) {
    xtimer_t timer;
    timer.reset();
    timer.start();
    UnitigGraph unitig_graph(&dbg);
    unitig_graph.InitFromSdBG();
    timer.stop();
    printf("Simple path graph size: %u, time for building: %lf\n", unitig_graph.size(), timer.elapsed());

    timer.reset();
    timer.start();
    histogram.clear();
    if (final_contig_file == NULL) {
        unitig_graph.OutputInitUnitigs(contigs_file, multi_file, histogram);
    } else {
        unitig_graph.OutputInitUnitigs(contigs_file, multi_file, final_contig_file, histogram, min_final_contig_len);
    }
    PrintStat();
    timer.stop();
    printf("Time to output: %lf\n", timer.elapsed());

    const double kMaxDepth = 65535;
    const int kLocalWidth = 1000;
    int64_t num_removed = 0;
    
    timer.reset();
    timer.start();
    while (min_depth < kMaxDepth) {
        // xtimer_t local_timer;
        // local_timer.reset();
        // local_timer.start();
        if (!unitig_graph.RemoveLocalLowDepth(min_depth, min_len, kLocalWidth, local_ratio, num_removed)) {
            break;
        }

        min_depth *= 1.1;
        // local_timer.stop();
        // printf("depth: %lf, num: %ld, time: %lf\n", min_depth, num_removed, local_timer.elapsed());
    }
    timer.stop();
    printf("Number of unitigs removed: %ld, time: %lf\n", num_removed, timer.elapsed());

    histogram.clear();
    unitig_graph.OutputChangedUnitigs(addi_contig_file, addi_multi_file, histogram);
    PrintStat();
}

void RemoveLowLocalAndOutputFinal(SuccinctDBG &dbg, FILE *final_contig_file, 
                                  double min_depth, int min_len, double local_ratio, int min_final_contig_len) {
    UnitigGraph unitig_graph(&dbg);
    unitig_graph.InitFromSdBG();
    printf("Simple path graph size: %u\n", unitig_graph.size());

    const double kMaxDepth = 65535;
    const int kLocalWidth = 1000;
    int64_t num_removed = 0;

    while (min_depth < kMaxDepth && 
           unitig_graph.RemoveLocalLowDepth(min_depth, min_len, kLocalWidth, local_ratio, num_removed)) {
        min_depth *= 1.1;
    }
    printf("Number of unitigs removed: %ld\n", num_removed);

    histogram.clear();
    unitig_graph.OutputFinalUnitigs(final_contig_file, histogram, min_final_contig_len);
    PrintStat();
}

void PrintStat(long long genome_size) {
    // total length
    int64_t total_length = 0;
    int64_t total_contigs = 0;
    int64_t average_length = 0;
    for (auto it = histogram.begin(); it != histogram.end(); ++it) {
        total_length += it->first * it->second;
        total_contigs += it->second;
    }
    if (genome_size == 0) { genome_size = total_length; }

    if (total_contigs > 0) {
        average_length = total_length / total_contigs;
    }

    // N50
    int64_t n50 = -1;
    int64_t acc_length = 0;
    for (auto it = histogram.rbegin(); it != histogram.rend(); ++it) {
        acc_length += it->first * it->second;
        if (n50 == -1 && acc_length * 2 >= genome_size) {
            n50 = it->first;
            break;
        }
    }

    printf("Total length: %ld, N50: %ld, Mean: %ld, number of contigs: %ld\n", total_length, n50, average_length, total_contigs);
    printf("Maximum length: %ld\n", histogram.size() > 0 ? histogram.rbegin()->first : 0);
}

static inline void MarkNode(SuccinctDBG &dbg, int64_t node_idx) {
    node_idx = dbg.GetLastIndex(node_idx);
    marked.set(node_idx);
}

} // namespace assembly_algorithms