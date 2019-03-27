/**
 * @file contig_graph.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-26
 */

#include "contig_graph.h"

#include <deque>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <sstream>

#include "contig_graph_branch_group.h"
#include "sequence.h"


using namespace std;


void ContigGraph::Initialize(const deque<Sequence> &contigs, const deque<ContigInfo> &contig_infos)
{
    vertices_.clear();
    vertices_.resize(contigs.size());

    for (int64_t i = 0; i < (int64_t)contigs.size(); ++i)
    {
        vertices_[i].clear();
        vertices_[i].set_contig(contigs[i]);
        vertices_[i].set_contig_info(contig_infos[i]);
        vertices_[i].set_id(i);
    }
    RefreshEdges();
}

void ContigGraph::Refresh()
{
    RefreshVertices();
    RefreshEdges();
}

void ContigGraph::RefreshVertices()
{
    uint64_t index = 0;
    for (unsigned i = 0; i < vertices_.size(); ++i)
    {
        if (!vertices_[i].status().IsDead())
        {
            vertices_[index].swap(vertices_[i]);
            vertices_[index].set_id(index);
            ++index;
        }
    }
    vertices_.resize(index);
}

void ContigGraph::RefreshEdges()
{
    BuildBeginIdbaKmerMap();

    uint64_t total_degree = 0;

    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        for (int strand = 0; strand < 2; ++strand)
        {
            ContigGraphVertexAdaptor current(&vertices_[i], strand);

            for (int x = 0; x < 4; ++x)
            {
                if (current.out_edges()[x])
                {
                    IdbaKmer kmer = current.end_kmer(kmer_size_);
                    kmer.ShiftAppend(x);
                    if (FindVertexAdaptorByBeginIdbaKmer(kmer).is_null())
                        current.out_edges().Remove(x);
                }
            }

            total_degree += current.out_edges().size();
        }

        if (vertices_[i].contig().size() == kmer_size_ 
                && vertices_[i].contig().IsPalindrome())
        {
            vertices_[i].in_edges() = vertices_[i].out_edges() | vertices_[i].out_edges();
            vertices_[i].out_edges() = vertices_[i].in_edges();
        }
    }

    num_edges_ = total_degree / 2;
}

void ContigGraph::ClearStatus()
{
    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
        vertices_[i].status().clear();
}

int64_t ContigGraph::Trim(int min_length)
{
    uint64_t old_num_vertices = vertices_.size();

    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        if (vertices_[i].contig().size() == kmer_size_
                && vertices_[i].contig().IsPalindrome())
            continue;

        if ((vertices_[i].in_edges().empty() || vertices_[i].out_edges().empty())
                && vertices_[i].contig().size() < min_length + kmer_size_ - 1
                && (vertices_[i].in_edges().size() + vertices_[i].out_edges().size() <= 1)
           )
        {
            vertices_[i].status().SetDeadFlag();
        }
    }
    Refresh();
    MergeSimplePaths();

    return old_num_vertices - vertices_.size();
}

int64_t ContigGraph::Trim(int min_length, double min_cover)
{
    uint64_t old_num_vertices = vertices_.size();

    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        if (vertices_[i].contig().size() == kmer_size_
                && vertices_[i].contig().IsPalindrome())
            continue;

        if ((vertices_[i].in_edges().empty() || vertices_[i].out_edges().empty())
                && vertices_[i].contig().size() < min_length + kmer_size_ - 1
                && (vertices_[i].in_edges().size() + vertices_[i].out_edges().size() <= 1
                && vertices_[i].coverage() < min_cover)
           )
        {
            vertices_[i].status().SetDeadFlag();
        }
    }
    Refresh();
    MergeSimplePaths();

    return old_num_vertices - vertices_.size();
}

int64_t ContigGraph::RemoveDeadEnd(int min_length)
{
    uint64_t num_deadend = 0;
    int l = 1;
    while (true)
    {
        l = min(2*l, min_length);
        num_deadend += Trim(l);

        if (l == min_length)
            break;
    }
    num_deadend += Trim(min_length);
    return num_deadend;
}

int64_t ContigGraph::RemoveBubble()
{
    deque<ContigGraphVertexAdaptor> candidates;

    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        for (int strand = 0; strand < 2; ++strand)
        {
            ContigGraphVertexAdaptor current(&vertices_[i], strand);

            if (current.out_edges().size() > 1 && current.contig_size() > kmer_size_)
            {
                ContigGraphBranchGroup branch_group(this, current, 4, kmer_size_ + 2);

                if (branch_group.Search())
                {
                    ContigGraphVertexAdaptor begin = branch_group.begin();
                    ContigGraphVertexAdaptor end = branch_group.end();

                    begin.ReverseComplement();
                    end.ReverseComplement();
                    std::swap(begin, end);
                    ContigGraphBranchGroup rev_branch_group(this, begin, 4, kmer_size_ + 2);

                    if (rev_branch_group.Search() && rev_branch_group.end() == end)
                    {
                        candidates.push_back(current);
                    }
                }
            }
        }
    }

    int64_t bubble = 0;
    for (unsigned i = 0; i < candidates.size(); ++i)
    {
        ContigGraphVertexAdaptor current = candidates[i];

        if (current.out_edges().size() > 1)
        {
            ContigGraphBranchGroup branch_group(this, current, 4, kmer_size_ + 2);

            if (branch_group.Search())
            {
                ContigGraphVertexAdaptor begin = branch_group.begin();
                ContigGraphVertexAdaptor end = branch_group.end();

                begin.ReverseComplement();
                end.ReverseComplement();
                std::swap(begin, end);
                ContigGraphBranchGroup rev_branch_group(this, begin, 4, kmer_size_ + 2);

                if (rev_branch_group.Search() && rev_branch_group.end() == end)
                {
                    branch_group.Merge();
                    ++bubble;
                }
            }
        }
    }

    Refresh();
    MergeSimplePaths();

    return bubble;
}

double ContigGraph::IterateCoverage(int min_length, double min_cover, double max_cover, double factor)
{
    min_cover = min(min_cover, max_cover);
    while (true)
    {
        RemoveLowCoverage(min_cover, min_length);
        min_cover *= factor;
        if (min_cover >= max_cover)
            break;
    }
    return min_cover;
}

bool ContigGraph::RemoveLowCoverage(double min_cover, int min_length)
{
    bool is_changed = false;

    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        ContigGraphVertexAdaptor current(&vertices_[i]);

        if (current.contig_size() < min_length + kmer_size_ - 1
                && ((current.in_edges().size() <= 1 && current.out_edges().size() <= 1)
                        || current.in_edges().size() == 0 || current.out_edges().size() == 0)
           )
        {
            if (current.coverage() < min_cover)
            {
                is_changed = true;
                current.status().SetDeadFlag();
            }
        }
    }

    Refresh();
    //Trim(min_length);
    MergeSimplePaths();

    return is_changed;
}

void ContigGraph::MergeSimplePaths()
{
    deque<Sequence> contigs;
    deque<ContigInfo> contig_infos;
    Assemble(contigs, contig_infos);
    Initialize(contigs, contig_infos);
}


int64_t ContigGraph::Assemble(deque<Sequence> &contigs, deque<ContigInfo> &contig_infos)
{
    contigs.clear();
    contig_infos.clear();

    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        if (vertices_[i].contig().size() == kmer_size_
                && vertices_[i].contig().IsPalindrome())
        {
            vertices_[i].status().Lock(1);

            Sequence contig = vertices_[i].contig();
            //ContigInfo contig_info(vertices_[i].kmer_count(), vertices_[i].in_edges(), vertices_[i].out_edges());
            ContigInfo contig_info;
            contig_info.set_kmer_count(vertices_[i].kmer_count());
            contig_info.in_edges() = vertices_[i].in_edges();
            contig_info.out_edges() = vertices_[i].out_edges();

            contigs.push_back(contig);
            contig_infos.push_back(contig_info);
        }
    }

    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        if (!vertices_[i].status().Lock(0))
            continue;

        ContigGraphPath path;
        path.Append(ContigGraphVertexAdaptor(&vertices_[i]), 0);

        Sequence contig;
        ContigInfo contig_info;
        for (int strand = 0; strand < 2; ++strand)
        {
            while (true)
            {
                ContigGraphVertexAdaptor current = path.back();
                ContigGraphVertexAdaptor next;

                if (!GetNextVertexAdaptor(current, next))
                    break;

                if (IsPalindromeLoop(path, next))
                    break;

                if (IsLoop(path, next))
                    goto FAIL;

                if (!next.status().LockPreempt(0))
                    goto FAIL;

                path.Append(next, -kmer_size_ + 1);
            }

            path.ReverseComplement();
        }

        path.Assemble(contig, contig_info);
        contigs.push_back(contig);
        contig_infos.push_back(contig_info);
FAIL:
        ;
    }

    ClearStatus();

    return contigs.size();
}

struct SearchNode
{
    ContigGraphVertexAdaptor node;
    int distance;
    int label;
};

bool ContigGraph::IsConverged(ContigGraphVertexAdaptor current)
{
    int TimeLimit = 1000;
    int DistanceLimit = 300;
    map<ContigGraphVertexAdaptor, int> reachable;
    queue<SearchNode> qu;

    for (int x = 0; x < 4; ++x)
    {
        if (current.out_edges()[x])
        {
            SearchNode search_node;
            search_node.node = GetNeighbor(current, x);
            search_node.distance = -(int)kmer_size_ + 1;
            search_node.label = x;

            //if (!search_node.node.status().IsDead())
                qu.push(search_node);
        }
    }

    int time = 0;
    while (!qu.empty())
    {
        if (time++ == TimeLimit)
            break;

        SearchNode search_node = qu.front();
        qu.pop();

        reachable[search_node.node] |= (1 << search_node.label);
        
        if (reachable[search_node.node] == (int)current.out_edges())
        {
            return true;
        }

        if (search_node.distance + (int)search_node.node.contig_size() - (int)kmer_size_ + 1 > DistanceLimit)
            continue;

        for (int x = 0; x < 4; ++x)
        {
            if (search_node.node.out_edges()[x])
            {
                ContigGraphVertexAdaptor next = GetNeighbor(search_node.node, x);

                SearchNode new_search_node;
                new_search_node.node = next;
                new_search_node.distance = search_node.distance + (int)search_node.node.contig_size() - (int)kmer_size_ + 1;
                new_search_node.label = search_node.label;

//                if (new_search_node.node == current)
//                    continue;

                if (reachable[new_search_node.node] & (1 << new_search_node.label))
                    continue;

                //if (!new_search_node.node.status().IsDead())
                    qu.push(new_search_node);
            }
        }
    }

    return false;
}

int64_t ContigGraph::SplitBranches()
{
    //cout << num_vertices() << " " << num_edges() << endl;

    deque<ContigGraphVertexAdaptor> branches;

    int64_t count = 0;
    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        ContigGraphVertexAdaptor current(&vertices_[i]);
        for (int strand = 0; strand < 2; ++strand)
        {
            if (!IsConverged(current))
            {
                ++count;

                branches.push_back(current);
            }

            current.ReverseComplement();
        }
    }

    set<ContigGraphVertexAdaptor> sources;

    for (unsigned i = 0; i < branches.size(); ++i)
        sources.insert(branches[i]);

    for (unsigned i = 0; i < branches.size(); ++i)
    {
        ContigGraphVertexAdaptor u = branches[i];

        for (int x = 0; x < 4; ++x)
        {
            if (u.out_edges()[x])
            {
                ContigGraphVertexAdaptor v = GetNeighbor(u, x);
                v.ReverseComplement();
                if (sources.find(v) == sources.end())
                {
                    sources.insert(v);
                    branches.push_back(v);
                }
                //RemoveEdge(u, x);
            }
        }
    }

    for (unsigned i = 0; i < branches.size(); ++i)
    {
        ContigGraphVertexAdaptor u = branches[i];

        for (int x = 0; x < 4; ++x)
        {
            if (u.out_edges()[x])
                RemoveEdge(u, x);
        }
    }

    RefreshEdges();

    return count;
}

void ContigGraph::GetComponents(deque<deque<ContigGraphVertexAdaptor> > &components, deque<string> &component_strings)
{
    components.clear();
    component_strings.clear();

    for (unsigned i = 0; i < vertices().size(); ++i)
    {
        if (vertices()[i].status().IsUsed())
            continue;

        deque<ContigGraphVertexAdaptor> qu;
        qu.push_back(ContigGraphVertexAdaptor(&vertices()[i], 0));
        vertices()[i].status().SetUsedFlag();

        stringstream ss;
        for (int index = 0; index < (int)qu.size(); ++index)
        {
            ContigGraphVertexAdaptor current = qu[index];

            for (int strand = 0; strand < 2; ++strand)
            {
                //for (connection_list_iterator p = connections()[current].begin(); p != connections()[current].end(); ++p)
                for (int x = 0; x < 4; ++x)
                {
                    if (current.out_edges()[x])
                    {
                        ContigGraphVertexAdaptor next = GetNeighbor(current, x);

                        if (strand == 0)
                        {
                            ss << current.id() << "_" << current.is_reverse() << "_" << current.contig_size() << "_" << current.kmer_count() << " " 
                                << next.id() << "_" << next.is_reverse() << "_" << next.contig_size() << "_" << next.kmer_count() << endl;

                            if (!next.status().IsUsed())
                                qu.push_back(next);
                        }
                        else
                        {
                            ss << next.id() << "_" << next.is_reverse() << "_" << next.contig_size() << "_" << next.kmer_count() << " "
                                << current.id() << "_" << current.is_reverse() << "_" << current.contig_size() << "_" << current.kmer_count() << endl;

                            if (!next.status().IsUsed())
                                qu.push_back(next.ReverseComplement());
                        }

                        next.status().SetUsedFlag();
                    }
                }

                current.ReverseComplement();
            }
        }

        components.push_back(qu);
        component_strings.push_back(ss.str());
    }

    ClearStatus();
}

double ContigGraph::GetSimilarity(const Sequence &a, const Sequence &b)
{
    vector<vector<int> > table;
    table.resize(a.size() + 1);
    for (unsigned i = 0; i < table.size(); ++i)
        table[i].resize(b.size() + 1);

    for (int i = 0; i <= (int)a.size(); ++i)
        table[i][0] = i;

    for (int j = 0; j <= (int)b.size(); ++j)
        table[0][j] = j;

    for (int i = 1; i <= (int)a.size(); ++i)
    {
        for (int j = 1; j <= (int)b.size(); ++j)
        {
            table[i][j] = 1000000000;
            if (table[i-1][j] + 1 < table[i][j])
                table[i][j] = table[i-1][j] + 1;
            if (table[i][j-1] + 1 < table[i][j])
                table[i][j] = table[i][j-1] + 1;
            if (table[i-1][j-1] + (a[i-1] != b[j-1]) < table[i][j])
                table[i][j] = table[i-1][j-1] + (a[i-1] != b[j-1]);
        }
    }

    return 1.0 - 1.0 * table[a.size()][b.size()] / max(a.size(), b.size());
}

void ContigGraph::BuildBeginIdbaKmerMap()
{
    begin_kmer_map_.clear();
    begin_kmer_map_.reserve(vertices_.size()*2);

    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        for (int strand = 0; strand < 2; ++strand)
        {
            ContigGraphVertexAdaptor current(&vertices_[i], strand);
            IdbaKmer kmer = current.begin_kmer(kmer_size_);

            IdbaKmer key = kmer.unique_format();
            begin_kmer_map_[key] = i;
        }
    }
}

bool ContigGraph::CycleDetect(ContigGraphVertexAdaptor current, map<int, int> &status)
{
    if (status[current.id()] == 0)
    {
        bool flag = false;
        status[current.id()] = 1;
        for (int x = 0; x < 4; ++x)
        {
            if (current.out_edges()[x])
            {
                if (CycleDetect(GetNeighbor(current, x), status))
                    flag = true;
            }
        }
        status[current.id()] = 2;
        return flag;
    }
    else if (status[current.id()] == 1)
        return true;
    else
        return false;
}

void ContigGraph::TopSortDFS(deque<ContigGraphVertexAdaptor> &order, ContigGraphVertexAdaptor current, map<int, int> &status)
{
    if (status[current.id()] == 0)
    {
        status[current.id()] = 1;
        for (int x = 0; x < 4; ++x)
        {
            if (current.out_edges()[x])
                TopSortDFS(order, GetNeighbor(current, x), status);
        }
        order.push_back(current);
    }
}

int ContigGraph::GetDepth(ContigGraphVertexAdaptor current, int depth, int &maximum, int min_length)
{
    if (depth > maximum)
        maximum = depth;

    if (maximum >= min_length)
        return min_length;

    deque<ContigGraphVertexAdaptor> neighbors;
    GetNeighbors(current, neighbors);
    for (unsigned i = 0; i < neighbors.size(); ++i)
    {
        if (neighbors[i].status().IsDead())
            continue;

        GetDepth(neighbors[i], depth - kmer_size_ + 1 + neighbors[i].contig_size(), maximum, min_length);
    }

    return min(maximum, min_length);
}

double ContigGraph::FindSimilarPath(ContigGraphVertexAdaptor target, ContigGraphPath &path, int &time)
{
    if (++time > 100)
        return 0;

    ContigGraphVertexAdaptor current = path.back();
    if (path.size() > 1.1 * target.contig_size())
        return 0;
    else if (current.end_kmer(kmer_size_-1) == target.end_kmer(kmer_size_-1)
            && current.out_edges() == target.out_edges())
    {
        Sequence contig;
        ContigInfo contig_info;
        path.Assemble(contig, contig_info);
        return GetSimilarity(target.contig(), contig);
    }
    else
    {
        double maximum = 0;
        deque<ContigGraphVertexAdaptor> neighbors;
        GetNeighbors(current, neighbors);
        for (unsigned i = 0; i < neighbors.size(); ++i)
        {
            path.Append(neighbors[i], -kmer_size_+1);
            double tmp = FindSimilarPath(target, path, time);
            path.Pop();
            if (tmp > maximum)
                maximum = tmp;
        }
        return maximum;
    }
}


