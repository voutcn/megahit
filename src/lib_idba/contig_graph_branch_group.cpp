/**
 * @file contig_graph_branch_group.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.9
 * @date 2011-12-27
 */

#include "lib_idba/contig_graph_branch_group.h"

#include <algorithm>
#include <cstdio>
#include <iostream>

using namespace std;

bool ContigGraphBranchGroup::Search()
{
    int kmer_size = contig_graph_->kmer_size();
    branches_.reserve(max_branches_);

    ContigGraphPath path;
    path.Append(begin_, 0);
    branches_.push_back(path);

    if ((int)begin_.out_edges().size() <= 1 
            || (int)begin_.out_edges().size() > max_branches_ 
            || (int)begin_.contig_size() == kmer_size)
        return false;

    bool is_converge = false;
    for (int k = 1; k < max_length_; ++k)
    {
        int num_branches = branches_.size();
        bool is_extend = false;
        for (int i = 0; i < num_branches; ++i)
        {
            if ((int)branches_[i].internal_size(kmer_size) >= max_length_)
                continue;

            ContigGraphVertexAdaptor current = branches_[i].back();

            if (current.out_edges().size() == 0)
                return false;

            bool is_first = true;
            ContigGraphPath path = branches_[i];
            for (int x = 0; x < 4; ++x)
            {
                if (current.out_edges()[x])
                {
                    ContigGraphVertexAdaptor next = contig_graph_->GetNeighbor(current, x);

                    if (next.status().IsDead())
                        return false;

                    if (is_first)
                    {
                        branches_[i].Append(next, -kmer_size + 1);
                        is_first = false;
                    }
                    else
                    {
                        if ((int)branches_.size() == max_branches_)
                            return false;

                        path.Append(next, -kmer_size + 1);
                        branches_.push_back(path);
                        path.Pop();
                    }

                    is_extend = true;
                }
            }
        }

        end_ = branches_[0].back();

        if ((int)end_.contig_size() > kmer_size)
        {
            is_converge = true;
            for (unsigned i = 0; i < branches_.size(); ++i)
            {
                if (branches_[i].back() != end_ 
                        || (int)branches_[i].internal_size(kmer_size) != max_length_)
                {
                    is_converge = false;
                    break;
                }
            }

            if (is_converge)
                break;
        }

        if (!is_extend)
            break;
    }

    return is_converge && begin_ != end_;
}


void ContigGraphBranchGroup::Merge()
{
    unsigned best = 0;
    for (unsigned i = 1; i < branches_.size(); ++i)
    {
        if (branches_[i].kmer_count() > branches_[best].kmer_count())
            best = i;
    }

    for (unsigned i = 0; i < branches_.size(); ++i)
    {
        ContigGraphPath &path = branches_[i];
        path.front().out_edges() = 0;
        path.back().in_edges() = 0;
        for (unsigned j = 1; j+1 < path.num_nodes(); ++j)
        {
            path[j].in_edges() = 0;
            path[j].out_edges() = 0;
            path[j].status().SetDeadFlag();
        }
    }

    ContigGraphPath &path = branches_[best];
    for (unsigned j = 1; j+1 < path.num_nodes(); ++j)
        path[j].status().ResetDeadFlag();

    for (unsigned j = 0; j+1 < path.num_nodes(); ++j)
        contig_graph_->AddEdge(path[j], path[j+1]);
}

