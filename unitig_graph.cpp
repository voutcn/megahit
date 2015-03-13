/*
 *  unitig_graph.cpp
 *  This is a part of MEGAHIT
 *  
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

#include "unitig_graph.h"

#include <omp.h>
#include <queue>
#include <set>
#include <vector>
#include <algorithm>
#include <string>
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
    if (i == j) { s[i] = Complement(s[i]); }
}

std::string VertexToDNAString(SuccinctDBG *sdbg_, UnitigGraphVertex &v) {
    static char acgt[] = "ACGT";
    std::string label;
    int64_t start_node = v.start_node;
    int64_t end_node = v.end_node;

    int64_t cur_node = end_node;
    while (cur_node != start_node) {
        int64_t prev_node =  assembly_algorithms::PrevSimplePathNode(*sdbg_, cur_node);
        cur_node = sdbg_->GetLastIndex(prev_node);
        int8_t cur_char = sdbg_->GetW(prev_node);
        assert(1 <= cur_char && cur_char <= 8);
        label.append(1, acgt[cur_char > 4 ? (cur_char - 5) : (cur_char - 1)]);
    }

    uint8_t seq[sdbg_->kMaxKmerK];
    sdbg_->Label(start_node, seq);
    for (int i = sdbg_->kmer_k - 1; i >= 0; --i) {
        assert(seq[i] >= 1 && seq[i] <= 4);
        label.append(1, acgt[seq[i] - 1]);
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
        assert(false);
    } else {
        int num_kmer = s.length() - kmer_k + 1;
        assert(num_kmer % 2 == 0);
        s.resize(num_kmer / 2 + (kmer_k - 1));
    }
}
// -- end of helper functions --

void UnitigGraph::InitFromSdBG() {
    start_node_map_.clear();
    vertices_.clear();

    omp_lock_t path_lock;
    omp_init_lock(&path_lock);
    AtomicBitVector marked(sdbg_->size);

    // assemble simple paths
#pragma omp parallel for
    for (int64_t node_idx = 0; node_idx < sdbg_->size; ++node_idx) {
        if (sdbg_->IsValidNode(node_idx) && 
            sdbg_->IsLast(node_idx) && 
            assembly_algorithms::NextSimplePathNode(*sdbg_, node_idx) == -1 &&
            marked.lock(node_idx)) {

            bool will_be_added = true;
            int64_t cur_node = node_idx, prev_node;
            int64_t depth = sdbg_->NodeMultiplicity(cur_node);
            uint32_t length = 1;
            while ((prev_node = assembly_algorithms::PrevSimplePathNode(*sdbg_, cur_node)) != -1) {
                cur_node = sdbg_->GetLastIndex(prev_node);
                if (!marked.lock(cur_node)) {
                    will_be_added = false;
                    break;
                }

                depth += sdbg_->NodeMultiplicity(cur_node);
                ++length;
            }

            if (!will_be_added) { continue; }

            int64_t rc_start = sdbg_->ReverseComplement(node_idx);
            if (rc_start == -1) {
                fprintf(stderr, "Node: %lld\n", (long long)node_idx); 
                fprintf(stderr, "Graph is incorrect!\n");
                exit(1);
            }
            int64_t rc_end = -1;

            if (!marked.lock(rc_start)) {
                rc_end = sdbg_->ReverseComplement(cur_node);
                // compare whose id is bigger
                if (std::max(node_idx, cur_node) < std::max(rc_start, rc_end)) {
                    will_be_added = false;
                }
            } else {
                // lock through the rc path
                int64_t rc_cur_node = rc_start;
                rc_end = rc_cur_node;
                bool extend_full = true;
                while ((rc_cur_node = assembly_algorithms::NextSimplePathNode(*sdbg_, rc_cur_node)) != -1) {
                    rc_cur_node = sdbg_->GetLastIndex(rc_cur_node);
                    rc_end = rc_cur_node;
                    if (!marked.lock(rc_cur_node)) {
                        extend_full = false;
                        break;
                    }
                }

                if (!extend_full) {
                    rc_end = sdbg_->ReverseComplement(cur_node);
                }
            }

            if (!will_be_added) { continue; }

            omp_set_lock(&path_lock);
            vertices_.push_back(UnitigGraphVertex(cur_node, node_idx, rc_start, rc_end, depth, length));
            omp_unset_lock(&path_lock);
        } // end if
    } // end for

    // assemble looped paths
#pragma omp parallel for
    for (int64_t node_idx = 0; node_idx < sdbg_->size; ++node_idx) {
        if (!marked.get(node_idx) && sdbg_->IsValidNode(node_idx) && sdbg_->IsLast(node_idx)) {
            omp_set_lock(&path_lock);
            if (!marked.get(node_idx)) {
                int64_t cur_node = node_idx;
                int64_t depth = sdbg_->NodeMultiplicity(node_idx);
                uint32_t length = 1;

                bool rc_marked = marked.get(sdbg_->ReverseComplement(node_idx)); // whether it is marked before entering the loop

                while (!marked.get(cur_node)) {
                    marked.set(cur_node);
                    int64_t prev_node = assembly_algorithms::PrevSimplePathNode(*sdbg_, cur_node);
                    assert(prev_node != -1);
                    cur_node = sdbg_->GetLastIndex(prev_node);
                    depth += sdbg_->NodeMultiplicity(cur_node);
                    ++length;
                }

                if (!rc_marked) {
                    vertices_.push_back(UnitigGraphVertex(cur_node, node_idx, 0, 0, depth, length));
                    vertices_.back().is_loop = true;
                    vertices_.back().is_deleted = true;
                    if (marked.get(sdbg_->ReverseComplement(node_idx))) {
                        // this loop is palindrome
                        vertices_.back().is_palindrome = true;
                        vertices_.back().length;
                    }
                }
            }
            omp_unset_lock(&path_lock);
        }
    }

    if (vertices_.size() >= kMaxNumVertices) {
        fprintf(stderr, "[ERROR] Too many vertices in the unitig graph (%llu >= %llu)\n", 
                (unsigned long long)vertices_.size(), (unsigned long long)kMaxNumVertices);
        exit(1);
    }

    // free memory for hash table construction
    sdbg_->FreeMul();
    {
        AtomicBitVector empty_abv;
        marked.swap(empty_abv);
    }

    start_node_map_.reserve(vertices_.size() * 2);

#pragma omp parallel for
    for (uint32_t i = 0; i < vertices_.size(); ++i) {
        if (!vertices_[i].is_deleted) {
            start_node_map_[vertices_[i].start_node] = i;
            start_node_map_[vertices_[i].rev_start_node] = i;        
        }
    }

    omp_destroy_lock(&path_lock);
}

bool UnitigGraph::RemoveLocalLowDepth(int min_depth, int min_len, int local_width, double local_ratio, int64_t &num_removed) {
    bool is_changed = false;
    bool need_refresh = false;

#pragma omp parallel for schedule(static, 1)
    for (uint32_t i = 0; i < vertices_.size(); ++i) {
        if (vertices_[i].is_deleted || vertices_[i].length >= min_len) { continue; }
        assert(vertices_[i].length > 0);

        int indegree = sdbg_->Indegree(vertices_[i].start_node);
        int outdegree = sdbg_->Outdegree(vertices_[i].end_node);

        if (indegree + outdegree == 0) { continue; }

        if ((indegree <= 1 && outdegree <= 1) || indegree == 0 || outdegree == 0) {
            double depth = (double)vertices_[i].depth / vertices_[i].length;
            if (is_changed && depth > min_depth)
                continue;
            
            double mean = LocalDepth_(vertices_[i], local_width);
            double threshold = min_depth;
            if (min_depth < mean * local_ratio)
                is_changed = true;
            else
                threshold = mean * local_ratio;

            if (depth < threshold) {
                is_changed = true;
                need_refresh = true;
                vertices_[i].is_dead = true;
#pragma omp atomic
                ++num_removed;
            }
        }
    }

    if (need_refresh) {
        Refresh_();
    }

    return is_changed;
}

double UnitigGraph::LocalDepth_(UnitigGraphVertex &vertex, int local_width) {
    double total_depth = 0;
    double num_added_kmer = 0;

    for (int dir = 0; dir < 2; ++dir) {
        int64_t outgoings[4];
        int outdegree = sdbg_->Outgoings(dir == 1 ? vertex.rev_end_node : vertex.end_node, outgoings);
        for (int i = 0; i < outdegree; ++i) {
            auto next_vertex_iter = start_node_map_.find(outgoings[i]);
            assert(next_vertex_iter != start_node_map_.end());
            UnitigGraphVertex &next_vertex = vertices_[next_vertex_iter->second];
            assert(!next_vertex.is_deleted);

            if (next_vertex.length <= local_width) {
                num_added_kmer += next_vertex.length;
                total_depth += next_vertex.depth;
            } else {
                num_added_kmer += local_width;
                total_depth += (double)next_vertex.depth * local_width / next_vertex.length;
            }
        }
    }

    if (num_added_kmer == 0) { return 0; }
    else { return total_depth / num_added_kmer; }
}

void UnitigGraph::Refresh_() {
    omp_lock_t reassemble_lock;
    omp_init_lock(&reassemble_lock);
    static AtomicBitVector marked;
    marked.reset(vertices_.size());

    // update the sdbg
#pragma omp parallel for
    for (uint32_t i = 0; i < vertices_.size(); ++i) {
        if (vertices_[i].is_dead && !vertices_[i].is_deleted) {
            int64_t cur_node = vertices_[i].end_node;
            while (cur_node != vertices_[i].start_node) {
                sdbg_->SetInvalid(cur_node);
                cur_node = sdbg_->UniqueIncoming(cur_node);
                assert(cur_node != -1);
                cur_node = sdbg_->GetLastIndex(cur_node);
            }
            sdbg_->SetInvalid(cur_node);

            if (vertices_[i].rev_end_node != vertices_[i].end_node) {
                cur_node = vertices_[i].rev_end_node;
                while (cur_node != vertices_[i].rev_start_node) {
                    sdbg_->SetInvalid(cur_node);
                    cur_node = sdbg_->UniqueIncoming(cur_node);
                    assert(cur_node != -1);
                    cur_node = sdbg_->GetLastIndex(cur_node);
                }
                sdbg_->SetInvalid(cur_node);
            }

            vertices_[i].is_deleted = true;
        }
    }

#pragma omp parallel for
    for (uint32_t i = 0; i < vertices_.size(); ++i) {
        if (vertices_[i].is_deleted) { continue; }
        int dir;
        if (assembly_algorithms::PrevSimplePathNode(*sdbg_, vertices_[i].start_node) == -1) {
            dir = 0;
        } else if (assembly_algorithms::PrevSimplePathNode(*sdbg_, vertices_[i].rev_start_node) == -1) {
            dir = 1;
        } else {
            continue;
        }

        if (!marked.lock(i)) { continue; }

        std::vector<std::pair<uint32_t, bool> > linear_path; // first: vertex_id, second: is_rc
        int64_t cur_end = dir == 0 ? vertices_[i].end_node : vertices_[i].rev_end_node;
        int64_t new_start = dir == 0 ? vertices_[i].start_node : vertices_[i].rev_start_node;
        int64_t new_rc_end = dir == 0 ? vertices_[i].rev_end_node : vertices_[i].end_node;

        while (true) {
            int64_t next_start = assembly_algorithms::NextSimplePathNode(*sdbg_, cur_end);
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

        if (linear_path.empty()) { continue; }

        if (!marked.lock(linear_path.back().first)) {
            if (linear_path.back().first > i) {
                marked.unset(i);
                continue;
            } else {
                while (!marked.lock(linear_path.back().first)) {
                    // wait for the other thread release the lock
                }
            }
        }

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
        } else {
            new_end = vertices_[linear_path.back().first].end_node;
            new_rc_start = vertices_[linear_path.back().first].rev_start_node;
        }

        vertices_[i].start_node = new_start;
        vertices_[i].end_node = new_end;
        vertices_[i].rev_start_node = new_rc_start;
        vertices_[i].rev_end_node = new_rc_end;
        vertices_[i].is_changed = true;
        if (i == linear_path.back().first) {
            vertices_[i].is_deleted = false;
        }
    }

    // looped path
#pragma omp parallel for
    for (uint32_t i = 0; i < vertices_.size(); ++i) {
        if (!vertices_[i].is_deleted && !marked.get(i)) {
            omp_set_lock(&reassemble_lock);
            if (!vertices_[i].is_deleted && !marked.get(i)) {
                uint32_t length = vertices_[i].length;
                int64_t depth = vertices_[i].depth;

                vertices_[i].is_changed = true;
                vertices_[i].is_loop = true;
                vertices_[i].is_deleted = true;

                for (int dir = 0; dir < 2; ++dir) {
                    int64_t cur_end = vertices_[i].end_node;
                    if (dir == 1) { cur_end = vertices_[i].rev_end_node; }
                    while (true) {
                        int64_t next_start = assembly_algorithms::NextSimplePathNode(*sdbg_, cur_end);
                        assert(next_start != -1);
                        auto next_vertex_iter = start_node_map_.find(next_start);
                        assert(next_vertex_iter != start_node_map_.end());
                        UnitigGraphVertex &next_vertex = vertices_[next_vertex_iter->second];
                        if (next_vertex.is_deleted) { break; }

                        length += next_vertex.length;
                        depth += next_vertex.depth;
                        next_vertex.is_deleted = true;
                    }
                }

                vertices_[i].depth = depth;
                vertices_[i].length = length;
            }
            omp_unset_lock(&reassemble_lock);
        }
    }

#pragma omp parallel for
    for (uint32_t i = 0; i < vertices_.size(); ++i) {
        if (!vertices_[i].is_deleted) {
            start_node_map_[vertices_[i].rev_start_node] = i;
        }
    }

    omp_destroy_lock(&reassemble_lock);
}

void UnitigGraph::OutputInitUnitigs(FILE *contig_file, FILE *multi_file, std::map<int64_t, int> &histo) {
    uint32_t output_id = 0;
    omp_lock_t output_lock;
    omp_init_lock(&output_lock);
    histo.clear();

#pragma omp parallel for
    for (unsigned i = 0; i < vertices_.size(); ++i) {
        uint16_t multi = std::min(kMaxMulti_t, int(vertices_[i].depth / vertices_[i].length + 0.5));
        std::string label = VertexToDNAString(sdbg_, vertices_[i]);

        if (vertices_[i].is_loop) {
            omp_set_lock(&output_lock);
            fprintf(contig_file, ">contig%d_length_%ld_multi_%d_loop\n%s\n", 
                                 output_id,
                                 label.length(),
                                 multi,
                                 label.c_str());
            fwrite(&multi, sizeof(uint16_t), 1, multi_file);
            ++output_id;
            ++histo[label.length()];
            omp_unset_lock(&output_lock);
        } else {
            int indegree = sdbg_->Indegree(vertices_[i].start_node);
            int outdegree = sdbg_->Outdegree(vertices_[i].end_node);
            if (indegree == 0 && outdegree == 0) {
                vertices_[i].is_deleted = true;
            }
            omp_set_lock(&output_lock);
            fprintf(contig_file, ">contig%d_length_%ld_multi_%d_in_%d_out_%d\n%s\n", 
                                 output_id, 
                                 label.length(), 
                                 multi, 
                                 indegree, 
                                 outdegree, 
                                 label.c_str());
            fwrite(&multi, sizeof(uint16_t), 1, multi_file);
            ++output_id;
            ++histo[label.length()];
            omp_unset_lock(&output_lock);
        }
    }

    omp_destroy_lock(&output_lock);
}

void UnitigGraph::OutputInitUnitigs(FILE *contig_file, 
                                    FILE *multi_file, 
                                    FILE *final_contig_file, 
                                    std::map<int64_t, int> &histo,
                                    int min_final_contig_length) {
    uint32_t output_id = 0;
    omp_lock_t output_lock;
    omp_init_lock(&output_lock);
    histo.clear();

#pragma omp parallel for
    for (unsigned i = 0; i < vertices_.size(); ++i) {
        uint16_t multi = std::min(kMaxMulti_t, int(vertices_[i].depth / vertices_[i].length + 0.5));
        std::string label = VertexToDNAString(sdbg_, vertices_[i]);

        if (vertices_[i].is_loop) {
            if (vertices_[i].is_palindrome) {
                FoldPalindrome(label, sdbg_->kmer_k, vertices_[i].is_loop);
            }
            if (label.length() < (unsigned)min_final_contig_length) {
                continue;
            }
            omp_set_lock(&output_lock);
            fprintf(final_contig_file, ">contig_%d_%d_length_%ld_multi_%d_loop\n%s\n", 
                                 sdbg_->kmer_k,
                                 output_id, 
                                 label.length(),
                                 multi,
                                 label.c_str());
            ++output_id;
            ++histo[label.length()];
            omp_unset_lock(&output_lock);
        } else {
            int indegree = sdbg_->Indegree(vertices_[i].start_node);
            int outdegree = sdbg_->Outdegree(vertices_[i].end_node);
            FILE *out_file = contig_file;
            if (indegree == 0 && outdegree == 0) {
                vertices_[i].is_deleted = true;
                if (vertices_[i].start_node == vertices_[i].rev_start_node) {
                    FoldPalindrome(label, sdbg_->kmer_k, vertices_[i].is_loop);
                }

                if (label.length() >= (unsigned)min_final_contig_length) {
                    out_file = final_contig_file;
                } else {
                    continue;
                }
            }

            omp_set_lock(&output_lock);
            fprintf(out_file, ">contig_%d_%d_length_%ld_multi_%d_in_%d_out_%d\n%s\n", 
                                 sdbg_->kmer_k,
                                 output_id, 
                                 label.length(), 
                                 multi, 
                                 indegree, 
                                 outdegree, 
                                 label.c_str());
            if (out_file == contig_file) {
                fwrite(&multi, sizeof(uint16_t), 1, multi_file);
            }
            ++output_id;
            ++histo[label.length()];
            omp_unset_lock(&output_lock);
        }
    }

    omp_destroy_lock(&output_lock);
}

void UnitigGraph::OutputChangedUnitigs(FILE *add_contig_file, FILE *addi_multi_file, std::map<int64_t, int> &histo) {
    uint32_t output_id = 0;
    omp_lock_t output_lock;
    omp_lock_t histo_lock;
    omp_init_lock(&output_lock);
    omp_init_lock(&histo_lock);
    histo.clear();

#pragma omp parallel for
    for (unsigned i = 0; i < vertices_.size(); ++i) {
        if ((vertices_[i].is_deleted && !vertices_[i].is_loop)) {
            continue;
        }

        uint16_t multi = std::min(kMaxMulti_t, int(vertices_[i].depth / vertices_[i].length + 0.5));
        std::string label = VertexToDNAString(sdbg_, vertices_[i]);

        omp_set_lock(&histo_lock);
        ++histo[label.length()];
        omp_unset_lock(&histo_lock);

        if (!vertices_[i].is_changed) { continue; }

        if (vertices_[i].is_loop) {
            if (label.length() >= (unsigned)sdbg_->kmer_k && 
                label.substr(label.length() - sdbg_->kmer_k + 1) == label.substr(0, sdbg_->kmer_k - 1)) {
                int num_vertex = label.length() - sdbg_->kmer_k + 1;
                if (num_vertex < sdbg_->kmer_k + 1) {
                    continue;
                }
                // WARN: hard code 28: the maximum step
                unsigned max_next_k = 28 + sdbg_->kmer_k;
                int j = sdbg_->kmer_k - 1;
                while (label.length() <= max_next_k + 1 ||
                       label.substr(0, max_next_k + 1) != label.substr(label.length() - max_next_k - 1)) {
                    label.push_back(label[j]);
                    ++j;
                }
            }
            omp_set_lock(&output_lock);
            fprintf(add_contig_file, ">addi%d_length_%ld_multi_%d_loop\n%s\n", 
                                 output_id, 
                                 label.length(),
                                 multi, 
                                 label.c_str());
            fwrite(&multi, sizeof(uint16_t), 1, addi_multi_file);
            ++output_id;
            omp_unset_lock(&output_lock);
        } else {
            int indegree = sdbg_->Indegree(vertices_[i].start_node);
            int outdegree = sdbg_->Outdegree(vertices_[i].end_node);
            omp_set_lock(&output_lock);
            fprintf(add_contig_file, ">addi%d_length_%ld_multi_%d_in_%d_out_%d\n%s\n", 
                                 output_id, 
                                 label.length(), 
                                 multi,
                                 indegree,
                                 outdegree,
                                 label.c_str());
            fwrite(&multi, sizeof(uint16_t), 1, addi_multi_file);
            ++output_id;
            omp_unset_lock(&output_lock);
        }

    }

    omp_destroy_lock(&histo_lock);
    omp_destroy_lock(&output_lock);
}


void UnitigGraph::OutputFinalUnitigs(FILE *final_contig_file,
                                     std::map<int64_t, int> &histo,
                                     int min_final_contig_length) {
    uint32_t output_id = 0;
    omp_lock_t output_lock;
    omp_init_lock(&output_lock);
    histo.clear();

#pragma omp parallel for
    for (unsigned i = 0; i < vertices_.size(); ++i) {
        if ((vertices_[i].is_deleted && !vertices_[i].is_loop)) {
            continue;
        }
        
        uint16_t multi = std::min(kMaxMulti_t, int(vertices_[i].depth / vertices_[i].length + 0.5));
        std::string label = VertexToDNAString(sdbg_, vertices_[i]);

        if (vertices_[i].is_loop) {
            if (vertices_[i].is_palindrome) {
                FoldPalindrome(label, sdbg_->kmer_k, vertices_[i].is_loop);
            }
            if (label.length() < (unsigned)min_final_contig_length) {
                continue;
            }
            omp_set_lock(&output_lock);
            fprintf(final_contig_file, ">contig_%d_%d_length_%ld_multi_%d_loop\n%s\n", 
                                 sdbg_->kmer_k,
                                 output_id, 
                                 label.length(),
                                 multi, 
                                 label.c_str());
            ++histo[label.length()];
            ++output_id;
            omp_unset_lock(&output_lock);
        } else if (vertices_[i].start_node == vertices_[i].rev_start_node) {
            // it is a palindrome
            FoldPalindrome(label, sdbg_->kmer_k, vertices_[i].is_loop);

            if (label.length() < (unsigned)min_final_contig_length) {
                continue;
            }
            omp_set_lock(&output_lock);
            fprintf(final_contig_file, ">contig_%d_%d_length_%ld_multi_%d_palindrome\n%s\n", 
                                 sdbg_->kmer_k,
                                 output_id, 
                                 label.length(),
                                 multi, 
                                 label.c_str());
            ++output_id;
            ++histo[label.length()];
            omp_unset_lock(&output_lock);
        } else {
            int indegree = sdbg_->Indegree(vertices_[i].start_node);
            int outdegree = sdbg_->Outdegree(vertices_[i].end_node);
            if (label.length() < (unsigned)min_final_contig_length) {
                continue;
            }
            omp_set_lock(&output_lock);
            fprintf(final_contig_file, ">contig_%d_%d_length_%ld_multi_%d_in_%d_out_%d\n%s\n", 
                                 sdbg_->kmer_k,
                                 output_id, 
                                 label.length(), 
                                 multi, 
                                 indegree, 
                                 outdegree, 
                                 label.c_str());
            ++output_id;
            ++histo[label.length()];
            omp_unset_lock(&output_lock);
        }
    }

    omp_destroy_lock(&output_lock);
}