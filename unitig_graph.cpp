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

#include "unitig_graph.h"

#include <omp.h>
#include <queue>
#include <set>
#include <vector>
#include <algorithm>
#include <string>
#include <utility>
#include <map>

#include "utils.h"
#include "definitions.h"
#include "succinct_dbg.h"
#include "assembly_algorithms.h"
#include "atomic_bit_vector.h"

// -- helper functions --
static inline char Complement(char c) {
    if (c >= 0 && c < 4) {
        return 3 - c;
    }

    switch (c) {
    case 'A': {
        return 'T';
    }

    case 'C': {
        return 'G';
    }

    case 'G': {
        return 'C';
    }

    case 'T': {
        return 'A';
    }

    default: {
        assert(false);
    }
    }
}

static inline void ReverseComplement(std::string &s) {
    int i, j;

    for (i = 0, j = s.length() - 1; i < j; ++i, --j) {
        std::swap(s[i], s[j]);
        s[i] = Complement(s[i]);
        s[j] = Complement(s[j]);
    }

    if (i == j) {
        s[i] = Complement(s[i]);
    }
}

std::string VertexToDNAString(SuccinctDBG *sdbg_, const UnitigGraphVertex &v) {
    std::string label;
    int64_t cur_edge = v.end_node;

    for (int i = 1; i < v.length; ++i) {
        int8_t cur_char = sdbg_->GetW(cur_edge);
        assert(1 <= cur_char && cur_char <= 8);
        label.append(1, "ACGT"[cur_char > 4 ? (cur_char - 5) : (cur_char - 1)]);

        cur_edge = sdbg_->PrevSimplePathEdge(cur_edge);
        // assert(cur_edge != -1);
        if (cur_edge == -1) {
            xerr_and_exit("%lld, %lld, %lld, %lld, (%lld, %lld), %d, %d\n", v.start_node, v.end_node, v.rev_start_node, v.rev_end_node, sdbg_->EdgeReverseComplement(v.end_node), sdbg_->EdgeReverseComplement(v.start_node), v.length, i);
        }
    }

    int8_t cur_char = sdbg_->GetW(cur_edge);
    label.append(1, "ACGT"[cur_char > 4 ? (cur_char - 5) : (cur_char - 1)]);

    if (cur_edge != v.start_node) {
        xerr("fwd: %lld, %lld, rev: %lld, %lld, (%lld, %lld) length: %d\n", v.start_node, v.end_node, v.rev_start_node, v.rev_end_node, sdbg_->EdgeReverseComplement(v.end_node), sdbg_->EdgeReverseComplement(v.start_node), v.length);
    }

    assert(cur_edge == v.start_node);

    uint8_t seq[sdbg_->kMaxKmerK];
    sdbg_->Label(v.start_node, seq);

    for (int i = sdbg_->kmer_k - 1; i >= 0; --i) {
        assert(seq[i] >= 1 && seq[i] <= 4);
        label.append(1, "ACGT"[seq[i] - 1]);
    }

    std::reverse(label.begin(), label.end());
    return label;
}

void FoldPalindrome(std::string &s, int kmer_k, bool is_loop) {
    if (is_loop) {
        for (unsigned i = 1; i + kmer_k <= s.length(); ++i) {
            std::string rc = s.substr(i, kmer_k);
            ReverseComplement(rc);

            if (rc == s.substr(i - 1, kmer_k)) {
                assert(i <= s.length() / 2);
                s = s.substr(i, s.length() / 2);
                break;
            }
        }
    }
    else {
        int num_edges = s.length() - kmer_k;
        assert(num_edges % 2 == 1);
        s.resize((num_edges - 1) / 2 + kmer_k + 1);
    }
}

void WriteContig(const std::string &label, int k_size, long long &id, int flag, double multiplicity, omp_lock_t *lock, FILE *file) {
    std::string rev_label(label);
    ReverseComplement(rev_label);

    omp_set_lock(lock);

    ++id;
    fprintf(file, ">k%d_%lld flag=%d multi=%.4lf len=%d\n%s\n",
            k_size,
            id,
            flag,
            multiplicity,
            (int)label.length(),
            label < rev_label ? label.c_str() : rev_label.c_str());

    omp_unset_lock(lock);
}

double GetSimilarity(std::string &a, std::string &b, double min_similar) {
    int n = a.length();
    int m = b.length();
    int max_indel = std::max(n, m) * (1 - min_similar);

    if (abs(n - m) > max_indel) {
        return 0;
    }

    if (max_indel < 1) {
        return 0;
    }

    std::vector<int> dp[2];

    for (int i = 0; i < 2; ++i) {
        dp[i].resize(max_indel * 2 + 1, 0);
    }

#define IDX(j, i) ((j) - (i) + max_indel)

    for (int j = 0; j <= max_indel; ++j) {
        dp[0][IDX(j, 0)] = j;
    }

    for (int i = 1; i <= n; ++i) {
        std::fill(dp[i & 1].begin(), dp[i & 1].end(), 99999999);

        if (i - max_indel <= 0) {
            dp[i & 1][IDX(0, i)] = i;
        }

        for (int j = std::max(i - max_indel, 1); j <= m && j <= i + max_indel; ++j) {
            // assert(IDX(j,i) >= 0 && IDX(j,i) < max_indel * 2 + 1);
            dp[i & 1][IDX(j, i)] = std::min(dp[i & 1][IDX(j, i)], dp[(i ^ 1) & 1][IDX(j - 1, i - 1)] + (a[i - 1] != b[j - 1]));

            if (j > i - max_indel) {
                // assert(IDX(j-1,i) >= 0 && IDX(j-1,i) < max_indel * 2 + 1);
                dp[i & 1][IDX(j, i)] = std::min(dp[i & 1][IDX(j, i)], dp[i & 1][IDX(j - 1, i)] + 1);
            }

            if (j < i + max_indel) {
                // assert(IDX(j,i-1) >= 0 && IDX(j,i-1) < max_indel * 2 + 1);
                dp[i & 1][IDX(j, i)] = std::min(dp[i & 1][IDX(j, i)], dp[(i ^ 1) & 1][IDX(j, i - 1)] + 1);
            }
        }
    }

    return 1 - dp[n & 1][IDX(m, n)] * 1.0 / std::max(n, m);
#undef IDX
}

// -- end of helper functions --

const size_t UnitigGraph::kMaxNumVertices = std::numeric_limits<UnitigGraph::vertexID_t>::max();

void UnitigGraph::InitFromSdBG() {
    start_node_map_.clear();
    vertices_.clear();

    omp_lock_t path_lock;
    omp_init_lock(&path_lock);
    AtomicBitVector marked(sdbg_->size);

    // assemble simple paths
#pragma omp parallel for

    for (int64_t edge_idx = 0; edge_idx < sdbg_->size; ++edge_idx) {
        if (sdbg_->IsValidEdge(edge_idx) && sdbg_->NextSimplePathEdge(edge_idx) == -1 && marked.try_lock(edge_idx)) {
            bool will_be_added = true;
            int64_t cur_edge = edge_idx, prev_edge;
            int64_t depth = sdbg_->EdgeMultiplicity(edge_idx);
            uint32_t length = 1;

            while ((prev_edge = sdbg_->PrevSimplePathEdge(cur_edge)) != -1) {
                cur_edge = prev_edge;

                if (!marked.try_lock(cur_edge)) {
                    will_be_added = false;
                    break;
                }

                depth += sdbg_->EdgeMultiplicity(cur_edge);
                ++length;
            }

            if (!will_be_added) {
                continue;
            }

            int64_t rc_start = sdbg_->EdgeReverseComplement(edge_idx);
            int64_t rc_end = -1;
            assert(rc_start != -1);

            if (!marked.try_lock(rc_start)) {
                rc_end = sdbg_->EdgeReverseComplement(cur_edge);

                if (std::max(edge_idx, cur_edge) < std::max(rc_start, rc_end)) {
                    will_be_added = false;
                }
            }
            else {
                // lock through the rc path
                int64_t rc_cur_edge = rc_start;
                rc_end = rc_cur_edge;
                bool extend_full = true;

                while ((rc_cur_edge = sdbg_->NextSimplePathEdge(rc_cur_edge)) != -1) {
                    rc_end = rc_cur_edge;

                    if (!marked.try_lock(rc_cur_edge)) {
                        extend_full = false;
                        break;
                    }
                }

                if (!extend_full) {
                    rc_end = sdbg_->EdgeReverseComplement(cur_edge);
                }
            }

            if (!will_be_added) {
                continue;
            }

            omp_set_lock(&path_lock);
            vertices_.push_back(UnitigGraphVertex(cur_edge, edge_idx, rc_start, rc_end, depth, length));
            omp_unset_lock(&path_lock);
        }
    }

    // assemble looped paths
#pragma omp parallel for

    for (int64_t edge_idx = 0; edge_idx < sdbg_->size; ++edge_idx) {
        if (!marked.get(edge_idx) && sdbg_->IsValidEdge(edge_idx)) {
            omp_set_lock(&path_lock);

            if (!marked.get(edge_idx)) {
                int64_t cur_edge = edge_idx;
                int64_t depth = sdbg_->EdgeMultiplicity(edge_idx);
                uint32_t length = 0;

                bool rc_marked = marked.get(sdbg_->EdgeReverseComplement(edge_idx)); // whether it is marked before entering the loop

                while (!marked.get(cur_edge)) {
                    marked.set(cur_edge);
                    depth += sdbg_->EdgeMultiplicity(cur_edge);
                    ++length;

                    cur_edge = sdbg_->PrevSimplePathEdge(cur_edge);
                    assert(cur_edge != -1);
                }

                assert(cur_edge == edge_idx);

                if (!rc_marked) {
                    int64_t start = sdbg_->NextSimplePathEdge(edge_idx);
                    int64_t end = edge_idx;
                    vertices_.push_back(UnitigGraphVertex(start, end, sdbg_->EdgeReverseComplement(end), sdbg_->EdgeReverseComplement(start), depth, length));
                    vertices_.back().is_loop = true;
                    vertices_.back().is_deleted = true; // loop path will not process to further steps, but still should be there for output

                    if (marked.get(sdbg_->EdgeReverseComplement(edge_idx))) {
                        // this loop is palindrome
                        vertices_.back().is_palindrome = true;
                    }
                }
            }

            omp_unset_lock(&path_lock);
        }
    }

    if (vertices_.size() >= kMaxNumVertices) {
        xerr_and_exit("[ERROR] Too many vertices in the unitig graph (%llu >= %llu)\n",
                      (unsigned long long)vertices_.size(), (unsigned long long)kMaxNumVertices);
    }

    // free memory for hash table construction
    sdbg_->FreeMul();
    {
        AtomicBitVector empty_abv;
        marked.swap(empty_abv);
    }

    start_node_map_.reserve(vertices_.size() * 2);

#pragma omp parallel for

    for (vertexID_t i = 0; i < vertices_.size(); ++i) {
        if (!vertices_[i].is_deleted) {
            start_node_map_[vertices_[i].start_node] = i;
            start_node_map_[vertices_[i].rev_start_node] = i;
        }
    }

    locks_.resize(vertices_.size());
#pragma omp parallel for

    for (vertexID_t i = 0; i < vertices_.size(); ++i) {
        omp_init_lock(&locks_[i]);
    }

    omp_destroy_lock(&path_lock);
}

uint32_t UnitigGraph::RemoveTips(int max_tip_len) {
    uint32_t num_removed = 0;
    for (int thre = 2; thre < max_tip_len; thre = std::min(thre * 2, max_tip_len)) {
#pragma omp parallel for schedule(static) reduction(+: num_removed)
        for (vertexID_t i = 0; i < vertices_.size(); ++i) {
            if (vertices_[i].is_deleted || vertices_[i].length >= thre) {
                continue;
            }

            int64_t outs[4], rev_outs[4];
            int outdegree = sdbg_->OutgoingEdges(vertices_[i].end_node, outs);
            int indegree = sdbg_->OutgoingEdges(vertices_[i].rev_end_node, rev_outs);

            if (indegree + outdegree == 0) {
                vertices_[i].is_dead = true;
            } else if (outdegree == 1 && indegree == 0) {
                if (sdbg_->PrevSimplePathEdge(outs[0]) == vertices_[i].end_node) {
                    vertices_[i].is_dead = true;
                }
            } else if (indegree == 1 && outdegree == 0) {
                if (sdbg_->PrevSimplePathEdge(rev_outs[0]) == vertices_[i].rev_end_node) {
                    vertices_[i].is_dead = true;
                }
            }

            num_removed += vertices_[i].is_dead;
        }

        Refresh_(false);

        if (thre >= max_tip_len) { break; }
    }

    return num_removed;
}

uint32_t UnitigGraph::MergeBubbles(bool permanent_rm, bool careful, FILE *bubble_file, Histgram<int64_t> &hist) {
    int max_bubble_len = sdbg_->kmer_k + 2; // allow 1 indel
    uint32_t num_removed = 0;
    omp_lock_t output_lock;
    omp_init_lock(&output_lock);
    long long output_id = 0;

    std::vector<std::tuple<double, int64_t, vertexID_t, int64_t> > branches; // depth, representative id, id, out_id

#pragma omp parallel for private(branches) reduction(+: num_removed)
    for (vertexID_t i = 0; i < vertices_.size(); ++i) {
        if (vertices_[i].is_deleted) {
            continue;
        }

        for (int strand = 0; strand < 2; ++strand) {
            int64_t outgoings[4];
            int outdegree = sdbg_->OutgoingEdges(strand == 0 ? vertices_[i].end_node : vertices_[i].rev_end_node, outgoings);

            if (outdegree <= 1) {
                continue;
            }

            branches.clear();
            bool converged = true;
            int max_len = -1, min_len = 99999999;
            int64_t next_outgoings[4];

            for (int j = 0; j < outdegree; ++j) {
                auto next_vertex_iter = start_node_map_.find(outgoings[j]);
                assert(next_vertex_iter != start_node_map_.end());
                UnitigGraphVertex &next_vertex = vertices_[next_vertex_iter->second];
                assert(!next_vertex.is_deleted);

                if (next_vertex.length > max_bubble_len) {
                    converged = false;
                    break;
                }

                if (next_vertex.start_node == outgoings[j] && sdbg_->EdgeOutdegree(next_vertex.rev_end_node) != 1) {
                    converged = false;
                    break;
                }

                if (next_vertex.rev_start_node == outgoings[j] && sdbg_->EdgeOutdegree(next_vertex.end_node) != 1) {
                    converged = false;
                    break;
                }

                if (sdbg_->OutgoingEdges(outgoings[j] == next_vertex.start_node ? next_vertex.end_node : next_vertex.rev_end_node, next_outgoings) != 1) {
                    converged = false;
                    break;
                }

                max_len = std::max(max_len, (int)next_vertex.length);
                min_len = std::min(min_len, (int)next_vertex.length);

                if (max_len - min_len > 2) {
                    converged = false;
                    break;
                }

                branches.push_back(std::make_tuple(-next_vertex.depth * 1.0 / next_vertex.length, next_vertex.Representation(),
                                                   next_vertex_iter->second, next_outgoings[0]));
            }

            for (int j = 1; converged && j < outdegree; ++j) {
                if (std::get<3>(branches[j]) != std::get<3>(branches[0])) {
                    converged = false;
                    break;
                }
            }

            if (!converged) {
                continue;
            }

            std::sort(branches.begin(), branches.end());
            bool careful_merged = false;

            for (int j = 1; j < outdegree; ++j) {
                UnitigGraphVertex &vj = vertices_[std::get<2>(branches[j])];
                vj.is_dead = true;
                ++num_removed;

                if (careful && -std::get<0>(branches[j]) >= -std::get<0>(branches[0]) * 0.2) {
                    careful_merged = true;
                    string label = VertexToDNAString(sdbg_, vj);
                    WriteContig(label, sdbg_->kmer_k, output_id, 0, -std::get<0>(branches[j]), &output_lock, bubble_file);
                    hist.insert(label.length());
                }
            }

            if (careful_merged) {
                UnitigGraphVertex &leftVertex = vertices_[i];
                UnitigGraphVertex &rightVertex = vertices_[start_node_map_[std::get<3>(branches[0])]];
                string left_label = VertexToDNAString(sdbg_, leftVertex);
                string right_label = VertexToDNAString(sdbg_, rightVertex);
                WriteContig(left_label, sdbg_->kmer_k, output_id, 0, leftVertex.depth * 1.0 / leftVertex.length, &output_lock, bubble_file);
                WriteContig(right_label, sdbg_->kmer_k, output_id, 0, rightVertex.depth * 1.0 / rightVertex.length, &output_lock, bubble_file);
                hist.insert(left_label.length());
                hist.insert(right_label.length());
            }
        }
    }
    omp_destroy_lock(&output_lock);

    Refresh_(!permanent_rm);
    return num_removed;
}

struct SuperBubbleNodeInfo {
    int dist : 30;
    int indegree : 3;
    UnitigGraph::vertexID_t from;
    double score; // sum of depths
    bool marked : 1;
    int max_dist : 30;

    SuperBubbleNodeInfo(int dist = 0, int indegree = 0, UnitigGraph::vertexID_t from = 0, int score = 0, int max_dist = 0):
        dist(dist), indegree(indegree), from(from), score(score), marked(false), max_dist(max_dist) {}
};

int UnitigGraph::SearchAndMergeSuperBubble_(int64_t source, int max_len, bool careful, FILE *bubble_file, Histgram<int64_t> &hist) {
    // see miniasm paper
    std::map<int64_t, SuperBubbleNodeInfo> visited;
    std::vector<int64_t> st;
    int num_inplay = 0;
    int num_too_long = 0;
    int num_tips = 0;
    vertexID_t vid = start_node_map_[source];

    if (sdbg_->EdgeOutdegree(vertices_[vid].start_node == source ? vertices_[vid].end_node : vertices_[vid].rev_end_node) <= 1) {
        return 0;
    }

    if (vertices_[vid].is_dead || !omp_test_lock(&locks_[vid])) {
        return -1;
    }

#define CLEAN_LOCK \
    do { \
        for (auto it = visited.begin(); it != visited.end(); ++it) { \
            omp_unset_lock(&locks_[start_node_map_[it->first]]); \
        } \
    } while (0); \

    visited[source].marked = true;
    st.push_back(source);

    while (!st.empty()) {
        int64_t cur = st.back();
        st.pop_back();
        SuperBubbleNodeInfo cur_node_info = visited[cur];

        auto &cur_vertex = vertices_[start_node_map_[cur]];

        int64_t outs[4];
        int outdegree = sdbg_->OutgoingEdges(cur_vertex.start_node == cur ? cur_vertex.end_node : cur_vertex.rev_end_node, outs);

        for (int j = 0; j < outdegree; ++j) {
            assert(start_node_map_.find(outs[j]) != start_node_map_.end());
            vertexID_t next = start_node_map_[outs[j]];
            if (outs[j] == source) {
                CLEAN_LOCK
                return 0;
            }

            auto iter = visited.find(outs[j]);
            if (iter == visited.end()) {
                // a new node
                if (vertices_[next].is_dead || !omp_test_lock(&locks_[next])) {
                    CLEAN_LOCK
                    return -1;
                }

                visited[outs[j]] = SuperBubbleNodeInfo(cur_node_info.dist + vertices_[next].length,
                                                 sdbg_->EdgeIndegree(outs[j]),
                                                 cur,
                                                 cur_node_info.score + vertices_[next].depth,
                                                 cur_node_info.max_dist + vertices_[next].length);

                if (visited.size() > max_len / 2) {
                    CLEAN_LOCK
                    return 0; // too complex
                }

                iter = visited.find(outs[j]);
                num_inplay++;
            } else {
                if (iter->second.score / iter->second.dist < (vertices_[next].depth + cur_node_info.score) / (cur_node_info.dist + vertices_[next].length)) {
                    iter->second.score = vertices_[next].depth + cur_node_info.score;
                    iter->second.dist = cur_node_info.dist + vertices_[next].length;
                    iter->second.from = cur;
                }

                if (iter->second.max_dist < cur_node_info.max_dist + vertices_[next].length) {
                    iter->second.max_dist = cur_node_info.max_dist + vertices_[next].length;
                }
            }

            if (--iter->second.indegree == 0) {
                if (!sdbg_->EdgeOutdegreeZero(outs[j] == vertices_[next].start_node ? vertices_[next].end_node : vertices_[next].rev_end_node)) {
                    st.push_back(outs[j]);
                } else {
                    // is a tip, we don't allow
                    CLEAN_LOCK
                    return 0;
                }

                if (iter->second.max_dist > max_len) {
                    if (++num_too_long > 1) {
                        CLEAN_LOCK
                        return 0;
                    }
                }

                num_inplay--;
            }
        }

        if (st.size() == 1 && num_inplay == 0) {
            int64_t sink = st.back();
            int64_t sink_len = vertices_[start_node_map_[sink]].length;
            visited[sink].dist -= sink_len;
            visited[sink].max_dist -= sink_len;

            if ((num_too_long > 0 && visited[sink].dist + sink_len <= max_len) ||
                visited[sink].dist < std::min(visited[sink].max_dist * 0.9, visited[sink].max_dist - 3.0)) {
                CLEAN_LOCK
                return -2;
            }

            double total_depth_remain = 0;
            double total_length_remain = 0;
            double max_avg_depth_removed = -.5;

            while (sink != source) {
                visited[sink].marked = true;
                if (careful && sink != st.back()) {
                    vertexID_t vid = start_node_map_[sink];
                    total_depth_remain += vertices_[vid].depth;
                    total_length_remain += vertices_[vid].length;
                }
                sink = visited[sink].from;
            }

            int num_removed = 0;

            if (careful) {
                for (auto iter = visited.begin(); iter != visited.end(); ++iter) {
                    if (!iter->second.marked) {
                        vertexID_t vid = start_node_map_[iter->first];
                        max_avg_depth_removed = std::max(max_avg_depth_removed, (double)vertices_[vid].depth / vertices_[vid].length);
                    }
                }

                assert(total_length_remain > 0);
                assert(max_avg_depth_removed >= 0);
                if (max_avg_depth_removed * total_length_remain > 0.2 * total_depth_remain) {
                    CLEAN_LOCK
                    return 0;
                }
            }

            for (auto iter = visited.begin(); iter != visited.end(); ++iter) {
                if (!iter->second.marked) {
                    vertices_[start_node_map_[iter->first]].is_dead = true;
                    num_removed++;
                }
            }
            CLEAN_LOCK
            return num_removed;
        }
    }

    CLEAN_LOCK
    return 0;

#undef CLEAN_LOCK
}

uint32_t UnitigGraph::MergeSuperBubbles(int max_len, bool permanent_rm, bool careful, FILE *bubble_file, Histgram<int64_t> &hist) {
    uint32_t num_merged = 0;
    uint32_t not_similar = 0;
    std::vector<vertexID_t> to_retry;
    omp_lock_t vec_lock;
    omp_init_lock(&vec_lock);

#pragma omp parallel for schedule(dynamic, 1) reduction(+: num_merged, not_similar)
    for (vertexID_t i = 0; i < vertices_.size(); ++i) {
        if (vertices_[i].is_deleted) continue;
        int search_res = SearchAndMergeSuperBubble_(vertices_[i].start_node, max_len, careful, bubble_file, hist);
        if (search_res > 0) {
            num_merged += search_res;
        } else if (search_res == -1) {
            omp_set_lock(&vec_lock);
            to_retry.push_back(i);
            omp_unset_lock(&vec_lock);
        } else {
            if (search_res == -2) not_similar++;
            int rev_res = SearchAndMergeSuperBubble_(vertices_[i].rev_start_node, max_len, careful, bubble_file, hist);
            if (rev_res == -1) {
                omp_set_lock(&vec_lock);
                to_retry.push_back(i);
                omp_unset_lock(&vec_lock);
            } else if (rev_res == -2) {
                not_similar++;
            } else {
                num_merged += rev_res;
            }
        }
    }

    Refresh_(!permanent_rm);

    for (auto it = to_retry.begin(); it != to_retry.end(); ++it) {
        vertexID_t i = *it;
        if (vertices_[i].is_deleted) continue;
        int search_res = SearchAndMergeSuperBubble_(vertices_[i].start_node, max_len, careful, bubble_file, hist);
        if (search_res > 0) {
            num_merged += search_res;
        } else if (search_res == 0) {
            int rev_res = SearchAndMergeSuperBubble_(vertices_[i].rev_start_node, max_len, careful, bubble_file, hist);
            if (rev_res > 0) {
                num_merged += rev_res;
            }
        }
    }

    Refresh_(!permanent_rm);

    xlog("Super bubble rejected b/c not similar: %u\n", not_similar);
    return num_merged;
}

uint32_t UnitigGraph::MergeComplexBubbles(double similarity, int merge_level, bool permanent_rm, bool careful, FILE *bubble_file, Histgram<int64_t> &hist) {
    int max_bubble_len = sdbg_->kmer_k * merge_level / similarity + 0.5;

    if (max_bubble_len * (1 - similarity) < 1) {
        return 0;
    }

    uint32_t num_removed = 0;

    std::vector<std::tuple<double, int64_t, vertexID_t, std::vector<int64_t>, bool> > branches; // depth, representative id, id, in_out_ids, strand
    std::vector<std::string> vertex_labels;

    long long output_id = 0;
    omp_lock_t output_lock;
    omp_init_lock(&output_lock);

#pragma omp parallel for private(branches, vertex_labels) reduction(+: num_removed)

    for (vertexID_t i = 0; i < vertices_.size(); ++i) {
        if (vertices_[i].is_deleted || vertices_[i].is_dead) {
            continue;
        }

        for (int strand = 0; strand < 2; ++strand) {
            int64_t outgoings[4];
            int outdegree = sdbg_->OutgoingEdges(strand == 0 ? vertices_[i].end_node : vertices_[i].rev_end_node, outgoings);

            if (outdegree <= 1) {
                continue;
            }

            branches.clear();
            vertex_labels.resize(outdegree);
            std::fill(vertex_labels.begin(), vertex_labels.end(), "");

            for (int j = 0; j < outdegree; ++j) {
                auto next_vertex_iter = start_node_map_.find(outgoings[j]);
                assert(next_vertex_iter != start_node_map_.end());
                UnitigGraphVertex &next_vertex = vertices_[next_vertex_iter->second];
                assert(!next_vertex.is_deleted);

                vector<int64_t> next_outgoings(8, -1);
                sdbg_->OutgoingEdges(outgoings[j] == next_vertex.start_node ? next_vertex.end_node : next_vertex.rev_end_node, &next_outgoings[0]);
                sdbg_->OutgoingEdges(outgoings[j] == next_vertex.start_node ? next_vertex.rev_end_node : next_vertex.end_node, &next_outgoings[4]);

                branches.push_back(std::make_tuple(-next_vertex.depth * 1.0 / next_vertex.length, next_vertex.Representation(),
                                                   next_vertex_iter->second, next_outgoings, outgoings[j] == next_vertex.start_node));
            }

            std::sort(branches.begin(), branches.end());
            std::vector<vertexID_t> left_or_right;

            bool careful_merged = false;
            for (int j = 0; j < outdegree; ++j) {
                UnitigGraphVertex &vj = vertices_[std::get<2>(branches[j])];

                if (vj.is_dead) {
                    continue;
                }

                if (vj.length > max_bubble_len) {
                    continue;
                }

                for (int k = j + 1; k < outdegree; ++k) {
                    UnitigGraphVertex &vk = vertices_[std::get<2>(branches[k])];

                    if (vk.is_dead) {
                        continue;
                    }

                    if (vk.length > max_bubble_len) {
                        continue;
                    }

                    if (std::get<3>(branches[j]) != std::get<3>(branches[k])) {
                        continue;
                    }

                    if ((vk.length + sdbg_->kmer_k - 1) * similarity <= (vj.length + sdbg_->kmer_k - 1) &&
                            (vj.length + sdbg_->kmer_k - 1) * similarity <= (vk.length + sdbg_->kmer_k - 1)) {
                        if (vertex_labels[j] == "") {
                            vertex_labels[j] = VertexToDNAString(sdbg_, std::get<4>(branches[j]) ? vj.ReverseComplement() : vj);
                        }

                        if (vertex_labels[k] == "") {
                            vertex_labels[k] = VertexToDNAString(sdbg_, std::get<4>(branches[k]) ? vk.ReverseComplement() : vk);
                        }

                        if (GetSimilarity(vertex_labels[j], vertex_labels[k], similarity) >= similarity) {
                            num_removed++;
                            vk.is_dead = true;

                            if (careful && -std::get<0>(branches[k]) >= -std::get<0>(branches[j]) * 0.2) {
                                careful_merged = true;
                                WriteContig(vertex_labels[k], sdbg_->kmer_k, output_id, 0, -std::get<0>(branches[j]), &output_lock, bubble_file);
                                hist.insert(vertex_labels[k].length());
                                for (int ni = 0; ni < 8; ++ni) {
                                    if (std::get<3>(branches[k])[ni] != -1) {
                                        left_or_right.push_back(start_node_map_[std::get<3>(branches[k])[ni]]);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if (careful_merged) {
                std::sort(left_or_right.begin(), left_or_right.end());
                for (int j = 0, sz = unique(left_or_right.begin(), left_or_right.end()) - left_or_right.begin();
                     j < sz; ++j) {
                    UnitigGraphVertex &vertex = vertices_[left_or_right[j]];
                    string label = VertexToDNAString(sdbg_, vertex);
                    WriteContig(label, sdbg_->kmer_k, output_id, 0, vertex.depth * 1.0 / vertex.length, &output_lock, bubble_file);
                    hist.insert(label.length());
                }
            }
        }
    }
    omp_destroy_lock(&output_lock);

    Refresh_(!permanent_rm);
    return num_removed;
}

int64_t UnitigGraph::RemoveLowDepth(double min_depth) {
    int64_t num_removed = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+: num_removed)

    for (vertexID_t i = 0; i < vertices_.size(); ++i) {
        if (!vertices_[i].is_deleted && vertices_[i].depth < min_depth) {
            vertices_[i].is_dead = true;
            ++num_removed;
        }
    }

    Refresh_(false);
    return num_removed;
}

bool UnitigGraph::RemoveLocalLowDepth(double min_depth, int min_len, int local_width, double local_ratio, int64_t &num_removed, bool permanent_rm) {
    bool is_changed = false;
    bool need_refresh = false;
    int64_t num_removed_ = 0;

#pragma omp parallel for schedule(dynamic, 1) reduction(+: num_removed_)

    for (vertexID_t i = 0; i < vertices_.size(); ++i) {
        if (vertices_[i].is_deleted || vertices_[i].length >= min_len) {
            continue;
        }

        assert(vertices_[i].length > 0);

        int indegree = sdbg_->EdgeIndegree(vertices_[i].start_node);
        int outdegree = sdbg_->EdgeOutdegree(vertices_[i].end_node);

        if (indegree + outdegree == 0) {
            continue;
        }

        if ((indegree <= 1 && outdegree <= 1) || indegree == 0 || outdegree == 0) {
            double depth = (double)vertices_[i].depth / vertices_[i].length;

            if (is_changed && depth > min_depth)
                continue;

            double mean = LocalDepth_(i, local_width);
            double threshold = min_depth;

            if (min_depth < mean * local_ratio)
                is_changed = true;
            else
                threshold = mean * local_ratio;

            if (depth < threshold) {
                is_changed = true;
                need_refresh = true;
                vertices_[i].is_dead = true;
                ++num_removed_;
            }
        }
    }

    if (need_refresh) {
        bool set_changed = !permanent_rm;
        Refresh_(set_changed);
    }

    num_removed += num_removed_;

    return is_changed;
}

uint32_t UnitigGraph::DisconnectWeakLinks(double local_ratio = 0.1) {
    // see metaspades paper
    vector<int64_t> to_remove;
    omp_lock_t lock;
    omp_init_lock(&lock);

#pragma omp parallel for schedule(dynamic, 1)
    for (vertexID_t i = 0; i < vertices_.size(); ++i) {
        if (vertices_[i].is_deleted) {
            continue;
        }

        for (int strand = 0; strand < 2; ++strand) {
            int64_t outgoings[4];
            double next_depths[4];
            double total_depth = 0;
            int outdegree = sdbg_->OutgoingEdges(strand == 0 ? vertices_[i].end_node : vertices_[i].rev_end_node, outgoings);

            if (outdegree <= 1) {
                continue;
            }


            for (int j = 0; j < outdegree; ++j) {
                auto next_vertex_iter = start_node_map_.find(outgoings[j]);
                assert(next_vertex_iter != start_node_map_.end());
                UnitigGraphVertex &next_vertex = vertices_[next_vertex_iter->second];
                assert(!next_vertex.is_deleted);

                next_depths[j] = next_vertex.depth * 1.0 / next_vertex.length;
                total_depth += next_depths[j];
            }

            for (int j = 0; j < outdegree; ++j) {
                if (next_depths[j] <= total_depth * local_ratio) {
                    omp_set_lock(&lock);
                    to_remove.push_back(outgoings[j]);
                    omp_unset_lock(&lock);
                }
            }
        }
    }

    omp_destroy_lock(&lock);

    uint32_t num_removed = to_remove.size();

#pragma omp parallel for schedule(dynamic, 1)
    for (vertexID_t i = 0; i < to_remove.size(); ++i) {
        int64_t id = to_remove[i];

        auto iter = start_node_map_.find(id);
        if (iter == start_node_map_.end()) { continue; }
        vertexID_t vid = iter->second;

        if (!omp_test_lock(&locks_[vid])) {
            continue;
        }

        if (vertices_[vid].length == 1 || vertices_[vid].start_node == vertices_[vid].rev_start_node) {
            vertices_[vid].is_dead = true;
        } else {
            int64_t rm1, rm2;
            start_node_map_.remove(id);

            if (vertices_[vid].start_node == id) {
                rm1 = vertices_[vid].start_node;
                rm2 = vertices_[vid].rev_end_node;

                vertices_[vid].start_node = sdbg_->NextSimplePathEdge(id);
                assert(vertices_[vid].start_node != -1);
                start_node_map_[vertices_[vid].start_node] = vid;
                vertices_[vid].rev_end_node = sdbg_->PrevSimplePathEdge(vertices_[vid].rev_end_node);

                assert(vertices_[vid].rev_end_node != -1);
                assert(rm1 != vertices_[vid].start_node);
                assert(rm2 != vertices_[vid].rev_end_node);

            } else {
                assert(vertices_[vid].rev_start_node == id);

                rm1 = vertices_[vid].rev_start_node;
                rm2 = vertices_[vid].end_node;

                vertices_[vid].rev_start_node = sdbg_->NextSimplePathEdge(id);
                assert(vertices_[vid].rev_start_node != -1);

                start_node_map_[vertices_[vid].rev_start_node] = vid;

                vertices_[vid].end_node = sdbg_->PrevSimplePathEdge(vertices_[vid].end_node);

                assert(vertices_[vid].end_node != -1);
                assert(rm1 != vertices_[vid].rev_start_node);
                assert(rm2 != vertices_[vid].end_node);
            }

            vertices_[vid].depth = (double)vertices_[vid].depth / vertices_[vid].length * (vertices_[vid].length - 1);
            vertices_[vid].length--;

            sdbg_->SetInvalidEdge(rm1);
            sdbg_->SetInvalidEdge(rm2);
        }
    }

#pragma omp parallel for
    for (vertexID_t i = 0; i < vertices_.size(); ++i) {
        omp_unset_lock(&locks_[i]);
    }

    Refresh_(false);
    return num_removed;
}

double UnitigGraph::LocalDepth_(vertexID_t id, int local_width) {
    double total_depth = 0;
    double num_added_edges = 0;

    for (int dir = 0; dir < 2; ++dir) {
        int64_t outgoings[4];
        int outdegree = sdbg_->OutgoingEdges(dir == 1 ? vertices_[id].rev_end_node : vertices_[id].end_node, outgoings);

        for (int i = 0; i < outdegree; ++i) {
            auto next_vertex_iter = start_node_map_.find(outgoings[i]);
            assert(next_vertex_iter != start_node_map_.end());
            UnitigGraphVertex &next_vertex = vertices_[next_vertex_iter->second];
            assert(!next_vertex.is_deleted);

            if (next_vertex.length <= local_width) {
                num_added_edges += next_vertex.length;
                total_depth += next_vertex.depth;
            }
            else {
                num_added_edges += local_width;
                total_depth += (double)next_vertex.depth * local_width / next_vertex.length;
            }
        }
    }

    if (num_added_edges == 0) {
        return 0;
    }
    else {
        return total_depth / num_added_edges;
    }
}

void UnitigGraph::Refresh_(bool set_changed) {
    omp_lock_t reassemble_lock;
    omp_init_lock(&reassemble_lock);

    // update the sdbg
#pragma omp parallel for

    for (vertexID_t i = 0; i < vertices_.size(); ++i) {
        if (vertices_[i].is_dead && !vertices_[i].is_deleted) {
            int64_t cur_edge = vertices_[i].end_node, prev_edge;

            while (cur_edge != vertices_[i].start_node) {
                prev_edge = sdbg_->UniquePrevEdge(cur_edge);
                sdbg_->SetInvalidEdge(cur_edge);
                cur_edge = prev_edge;
                assert(cur_edge != -1);
            }

            sdbg_->SetInvalidEdge(cur_edge);

            if (vertices_[i].rev_end_node != vertices_[i].end_node) {
                cur_edge = vertices_[i].rev_end_node;

                while (cur_edge != vertices_[i].rev_start_node) {
                    prev_edge = sdbg_->UniquePrevEdge(cur_edge);
                    sdbg_->SetInvalidEdge(cur_edge);
                    cur_edge = prev_edge;
                    assert(cur_edge != -1);
                }

                sdbg_->SetInvalidEdge(cur_edge);
            }

            vertices_[i].is_deleted = true;
        }
    }

#pragma omp parallel for

    for (vertexID_t i = 0; i < vertices_.size(); ++i) {
        if (vertices_[i].is_deleted) {
            continue;
        }

        int dir;

        if (sdbg_->PrevSimplePathEdge(vertices_[i].start_node) == -1) {
            dir = 0;
        }
        else if (sdbg_->PrevSimplePathEdge(vertices_[i].rev_start_node) == -1) {
            dir = 1;
        }
        else {
            continue;
        }

        if (!omp_test_lock(&locks_[i])) {
            continue;
        }

        std::vector<std::pair<vertexID_t, bool> > linear_path; // first: vertex_id, second: is_rc
        int64_t cur_end = dir == 0 ? vertices_[i].end_node : vertices_[i].rev_end_node;
        int64_t new_start = dir == 0 ? vertices_[i].start_node : vertices_[i].rev_start_node;
        int64_t new_rc_end = dir == 0 ? vertices_[i].rev_end_node : vertices_[i].end_node;

        while (true) {
            int64_t next_start = sdbg_->NextSimplePathEdge(cur_end);

            if (next_start == -1) {
                break;
            }

            auto next_vertex_iter = start_node_map_.find(next_start);
            assert(next_vertex_iter != start_node_map_.end());
            UnitigGraphVertex &next_vertex = vertices_[next_vertex_iter->second];
            assert(!next_vertex.is_deleted);

            bool is_rc = next_vertex.start_node != next_start;
            linear_path.push_back(std::make_pair(next_vertex_iter->second, is_rc));

            cur_end = is_rc ? next_vertex.rev_end_node : next_vertex.end_node;
        }

        if (linear_path.empty()) {
            vertices_[i].is_marked = true;
            continue;
        }

        if (i != linear_path.back().first && !omp_test_lock(&locks_[linear_path.back().first])) { // if i == linear_path.back().first, it is a palindrome self loop
            if (linear_path.back().first > i) {
                omp_unset_lock(&locks_[i]);
                continue;
            }
            else {
                omp_set_lock(&locks_[linear_path.back().first]);
            }
        }

        vertices_[i].is_marked = true;

        // assemble the linear path

        int64_t depth = vertices_[i].depth;
        int64_t length = vertices_[i].length;

        for (unsigned j = 0; j < linear_path.size(); ++j) {
            UnitigGraphVertex &next_vertex = vertices_[linear_path[j].first];
            length += next_vertex.length;
            depth += next_vertex.depth;
            next_vertex.is_deleted = true;
        }

        vertices_[i].length = length;
        vertices_[i].depth = depth;

        int64_t new_end;
        int64_t new_rc_start;

        if (linear_path.back().second) {
            new_end = vertices_[linear_path.back().first].rev_end_node;
            new_rc_start = vertices_[linear_path.back().first].start_node;
        }
        else {
            new_end = vertices_[linear_path.back().first].end_node;
            new_rc_start = vertices_[linear_path.back().first].rev_start_node;
        }

        vertices_[i].start_node = new_start;
        vertices_[i].end_node = new_end;
        vertices_[i].rev_start_node = new_rc_start;
        vertices_[i].rev_end_node = new_rc_end;
        vertices_[i].is_changed |= set_changed;

        if (i == linear_path.back().first) {
            vertices_[i].is_deleted = false;
        }
    }

    // looped path
#pragma omp parallel for

    for (vertexID_t i = 0; i < vertices_.size(); ++i) {
        if (!vertices_[i].is_deleted && !vertices_[i].is_marked) {
            omp_set_lock(&reassemble_lock);

            if (!vertices_[i].is_deleted) {
                uint32_t length = vertices_[i].length;
                int64_t depth = vertices_[i].depth;

                vertices_[i].is_changed |= set_changed;
                vertices_[i].is_loop = true;
                vertices_[i].is_deleted = true;
                bool is_palindrome = false;

                int64_t cur_end = vertices_[i].end_node;

                while (true) {
                    int64_t next_start = sdbg_->NextSimplePathEdge(cur_end);
                    assert(next_start != -1);

                    if (next_start == vertices_[i].start_node) {
                        break;
                    }

                    auto next_vertex_iter = start_node_map_.find(next_start);
                    assert(next_vertex_iter != start_node_map_.end());
                    UnitigGraphVertex &next_vertex = vertices_[next_vertex_iter->second];

                    if (next_vertex.is_deleted) {
                        // that means the loop has alrealy gone through its rc
                        is_palindrome = true;
                    }

                    length += next_vertex.length;
                    depth += next_vertex.depth;
                    next_vertex.is_deleted = true;

                    cur_end = (next_vertex.start_node == next_start) ? next_vertex.end_node : next_vertex.rev_end_node;
                }

                vertices_[i].depth = depth;
                vertices_[i].length = length;
                vertices_[i].is_palindrome = is_palindrome;
                vertices_[i].end_node = sdbg_->PrevSimplePathEdge(vertices_[i].start_node);
                vertices_[i].rev_start_node = sdbg_->EdgeReverseComplement(vertices_[i].end_node);
                vertices_[i].rev_end_node = sdbg_->EdgeReverseComplement(vertices_[i].start_node);
            }

            omp_unset_lock(&reassemble_lock);
        }
    }

#pragma omp parallel for

    for (vertexID_t i = 0; i < vertices_.size(); ++i) {
        if (!vertices_[i].is_deleted) {
            vertices_[i].is_marked = false;
            start_node_map_[vertices_[i].rev_start_node] = i;
        }
        else {
            assert(!vertices_[i].is_marked);
        }
    }

#pragma omp parallel for
    for (vertexID_t i = 0; i < vertices_.size(); ++i) {
        omp_unset_lock(&locks_[i]);
    }

    omp_destroy_lock(&reassemble_lock);
}

void UnitigGraph::OutputContigs(FILE *contig_file, FILE *final_file, Histgram<int64_t> &histo,
                                bool change_only, int min_final_standalone) {
    omp_lock_t output_lock;
    long long output_id = 0;

    omp_init_lock(&output_lock);
    histo.clear();

    assert(!(change_only && final_file != NULL)); // if output changed contigs, must not output final contigs

#pragma omp parallel for

    for (vertexID_t i = 0; i < vertices_.size(); ++i) {
        if (vertices_[i].is_deleted && !vertices_[i].is_loop) {
            continue;
        }

        double multi = std::min((double)kMaxMulti_t, (double)vertices_[i].depth / vertices_[i].length);

        if (change_only) {
            multi = 1;
        }

        std::string label = VertexToDNAString(sdbg_, vertices_[i]);

        if (vertices_[i].is_palindrome) {
            FoldPalindrome(label, sdbg_->kmer_k, vertices_[i].is_loop);
        }

        histo.insert(label.length());

        if (change_only && !vertices_[i].is_changed) {
            continue;
        }

        if (vertices_[i].is_loop) {
            int flag = contig_flag::kLoop | contig_flag::kIsolated;
            FILE *out_file = contig_file;

            if (vertices_[i].is_palindrome) {
                flag = contig_flag::kIsolated;
            }

            if (final_file != NULL) {
                if (label.length() < (unsigned)min_final_standalone) {
                    continue;
                }
                else {
                    out_file = final_file;
                }
            }

            WriteContig(label, sdbg_->kmer_k, output_id, flag, multi, &output_lock, out_file);

        }
        else {
            FILE *out_file = contig_file;
            int flag = 0;

            int indegree = sdbg_->EdgeIndegree(vertices_[i].start_node);
            int outdegree = sdbg_->EdgeOutdegree(vertices_[i].end_node);

            if (indegree == 0 && outdegree == 0) {
                vertices_[i].is_deleted = true;

                if (vertices_[i].start_node == vertices_[i].rev_start_node) {
                    FoldPalindrome(label, sdbg_->kmer_k, vertices_[i].is_loop);
                }

                flag = contig_flag::kIsolated;

                if (final_file != NULL) {
                    if (label.length() < (unsigned)min_final_standalone) {
                        continue;
                    }
                    else {
                        out_file = final_file;
                    }
                }
            }

            WriteContig(label, sdbg_->kmer_k, output_id, flag, multi, &output_lock, out_file);
        }
    }

    omp_destroy_lock(&output_lock);
}