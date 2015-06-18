/**
 * @file contig_graph.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-26
 */

#include "lib_idba/contig_graph.h"

#include <deque>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <sstream>

#include "lib_idba/contig_graph_branch_group.h"
#include "lib_idba/sequence.h"


using namespace std;


double ContigGraph::Binormial(int n, int m)
{
    double product = 1;
    for (int i = 1; i <= n; ++i)
        product *= i;
    for (int i = 1; i <= m; ++i)
        product /= i;
    return product;
}


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

void ContigGraph::AddAllEdges()
{
    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        vertices_[i].in_edges() = 15;
        vertices_[i].out_edges() = 15;
    }
    RefreshEdges();
}

void ContigGraph::RemoveAllEdges()
{
    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        vertices_[i].in_edges() = 0;
        vertices_[i].out_edges() = 0;
    }
    RefreshEdges();
}

void ContigGraph::ClearStatus()
{
    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
        vertices_[i].status().clear();
}

void ContigGraph::MergeSimilarPath()
{
    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        for (int strand = 0; strand < 2; ++strand)
        {
            ContigGraphVertexAdaptor current(&vertices_[i], strand);

            if (current.status().IsDead())
                continue;

            if (current.out_edges().size() > 1)
            {
                deque<ContigGraphVertexAdaptor> neighbors;
                GetNeighbors(current, neighbors);
                sort(neighbors.begin(), neighbors.end(), CompareContigCoverage);

                for (unsigned j = 0; j < neighbors.size(); ++j)
                {
                    if (neighbors[j].status().IsDead())
                        continue;

                    for (unsigned k = j+1; k < neighbors.size(); ++k)
                    {
                        if (!neighbors[k].status().IsDead()
                                && neighbors[j].in_edges() == neighbors[k].in_edges() 
                                && neighbors[j].out_edges() == neighbors[k].out_edges()
                                && neighbors[j].begin_kmer(kmer_size_-1) == neighbors[k].begin_kmer(kmer_size_-1)
                                && neighbors[j].end_kmer(kmer_size_-1) == neighbors[k].end_kmer(kmer_size_-1)
                                && GetSimilarity(neighbors[j], neighbors[k]) > 0.98)
                        {
                            neighbors[k].status().SetDeadFlag();
                        }
                    }
                }
            }
        }
    }
    Refresh();
    MergeSimplePaths();

}

int64_t ContigGraph::Prune(int min_length)
{
    uint64_t old_num_vertices = vertices_.size();

    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        for (int strand = 0; strand < 2; ++strand)
        {
            ContigGraphVertexAdaptor current(&vertices_[i], strand);

            if (current.status().IsDead())
                continue;

            if (current.out_edges().size() <= 1)
                continue;

            int maximum = 0;
            int depth = GetDepth(current, kmer_size_ - 1, maximum, min_length + kmer_size_ - 1);
            if (depth > min_length + (int)kmer_size_ - 1)
                depth = min_length + (int)kmer_size_ - 1;

            deque<ContigGraphVertexAdaptor> neighbors;
            GetNeighbors(current, neighbors);
            for (unsigned j = 0; j < neighbors.size(); ++j)
            {
                if (neighbors[j].in_edges().size() == 1 
                        && neighbors[j].out_edges().size() == 0
                        && (int)neighbors[j].contig_size() < depth)
                    neighbors[j].status().SetDeadFlag();
            }
        }
    }
    Refresh();
    MergeSimplePaths();

    return old_num_vertices - vertices_.size();
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

int64_t ContigGraph::RemoveStandAlone(int min_length)
{
    uint64_t old_num_vertices = vertices_.size();

    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        if (vertices_[i].contig().size() == kmer_size_
                && vertices_[i].contig().IsPalindrome())
            continue;

        if ((vertices_[i].in_edges().empty() && vertices_[i].out_edges().empty())
                && vertices_[i].contig().size() < min_length + kmer_size_ - 1
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

int64_t ContigGraph::RemoveDeadEnd(int min_length, double min_cover)
{
    uint64_t num_deadend = 0;
    int l = 1;
    while (true)
    {
        l = min(2*l, min_length);
        num_deadend += Trim(l, min_cover);

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

double ContigGraph::IterateLocalCoverage(int min_length, double ratio, double min_cover, double max_cover, double factor)
{
    in_kmer_count_table_.reserve(vertices_.size());

    min_cover = min(min_cover, max_cover);
    while (true)
    {
        bool is_changed = RemoveLocalLowCoverage(min_cover, min_length, ratio);

        if (!is_changed)
            break;

        if (min_cover >= max_cover)
            break;

        min_cover *= factor;
    }
    return min_cover;
}

double ContigGraph::IterateComponentCoverage(int min_length, double ratio, double min_cover, double max_cover, double factor, int max_component_size)
{
    in_kmer_count_table_.reserve(vertices_.size());

    min_cover = min(min_cover, max_cover);
    while (true)
    {
        bool is_changed = RemoveComponentLowCoverage(min_cover, min_length, ratio, max_component_size);

        if (!is_changed)
            break;

        if (min_cover >= max_cover)
            break;

        min_cover *= factor;
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

bool ContigGraph::RemoveLocalLowCoverage(double min_cover, int min_length, double ratio)
{
    int region_length = 1000;
    //int region_length = 100;
    bool is_changed = false;

    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        ContigGraphVertexAdaptor current(&vertices_[i]);

        if (current.contig_size() < min_length + kmer_size_ - 1
                && ((current.in_edges().size() <= 1 && current.out_edges().size() <= 1)
                        || current.in_edges().size() == 0 || current.out_edges().size() == 0)
           )
        {
            if (is_changed && current.coverage() > min_cover)
                continue;

            double mean = LocalCoverage(current, region_length);
            double threshold = min_cover;
            if (min_cover < mean * ratio)
                is_changed = true;
            else
                threshold = mean * ratio;

            if (current.coverage() < threshold)
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

bool ContigGraph::RemoveComponentLowCoverage(double min_cover, int min_length, double ratio, int max_component_size)
{
    int region_length = 300;
    deque<deque<ContigGraphVertexAdaptor> > components;
    deque<string> component_strings;
    GetComponents(components, component_strings);

    deque<double> average_coverage(components.size());
    deque<int> component_id_table(vertices_.size());

    for (int64_t i = 0; i < (int64_t)components.size(); ++i)
    {
        double total_kmer_count = 0;
        double total = 0;

        for (unsigned j = 0; j < components[i].size(); ++j)
        {
            total_kmer_count += components[i][j].kmer_count();
            total += components[i][j].contig_size() - kmer_size_ + 1;
            component_id_table[components[i][j].id()] = i;
        }

        average_coverage[i] = total_kmer_count / total;
    }

    bool is_changed = false;
    //int max_component_size = 30;

    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        ContigGraphVertexAdaptor current(&vertices_[i]);
        int id = component_id_table[current.id()];

        if (components[id].size() <= 10)
            continue;

        if (current.contig_size() < min_length + kmer_size_ - 1
                && (current.in_edges().size() <= 1 && current.out_edges().size() <= 1)
                        //|| current.in_edges().size() == 0 || current.out_edges().size() == 0)
           )
        {
            if (is_changed && current.coverage() > min_cover)
                continue;

            double threshold = min_cover;
            double mean = LocalCoverage(current, region_length);
            //double mean = average_coverage[id];
            if (min_cover < ratio * mean || ((int)components[id].size() > max_component_size && min_cover < average_coverage[id]))
                is_changed = true;
            else
                threshold = ratio * mean;

            if (current.coverage() < threshold 
                    || ((int)components[id].size() > max_component_size && current.coverage() < average_coverage[id]))
            {
                is_changed = true;
                current.status().SetDeadFlag();
            }
        }
    }

    Refresh();
    MergeSimplePaths();

    return is_changed;
}

double ContigGraph::LocalCoverage(ContigGraphVertexAdaptor current, int region_length)
{
    double num_count = 0;
    int num_kmer = 0;
    LocalCoverageSingle(current, region_length, num_count, num_kmer);
    LocalCoverageSingle(current.ReverseComplement(), region_length, num_count, num_kmer);

    if (num_kmer == 0)
        //return 1e100;
        return 0;
    else
        return num_count / num_kmer;
}

double ContigGraph::LocalCoverageSingle(ContigGraphVertexAdaptor current, int region_length, double &total_count, int &total_kmer)
{
    map<int, int> visited;
    deque<ContigGraphVertexAdaptor> qu;
    qu.push_back(current);
    visited[current.id()] = 0;

    int index = 0;
    int num_added = 0;
    int num_count = 0;
    int num_kmer = 0;
    while (index < (int)qu.size())
    {
        current = qu[index++];

        if (num_added >= 4 * region_length)
            break;

        if (visited.size() > 32)
            break;

        if (visited[current.id()] >= region_length)
            continue;

        int dist = visited[current.id()];

        for (int x = 0; x < 4; ++x)
        {
            if (current.out_edges()[x])
            {
                ContigGraphVertexAdaptor next = GetNeighbor(current, x);
                if (visited.find(next.id()) == visited.end())
                {
                    visited[next.id()] = dist + next.num_kmer();
                    qu.push_back(next);

                    if ((int)next.num_kmer() + dist > region_length)
                    {
                        if ((int)next.num_kmer() < region_length)
                        {
                            num_count += (int64_t)next.kmer_count() * (region_length - dist) / next.num_kmer();
                            num_kmer += region_length - dist;
                            num_added += region_length - dist;
                        }
                        else
                        {
                            IdbaKmer begin = next.begin_kmer(kmer_size_);
                            if (in_kmer_count_table_.find(begin) == in_kmer_count_table_.end())
                            {
                                int in_kmer_count = 0;
                                for (int i = 0; i < region_length; ++i)
                                    in_kmer_count += next.get_count(i);
                                in_kmer_count_table_[begin] = in_kmer_count;
                            }

                            num_count += (int64_t)in_kmer_count_table_[begin] * (region_length - dist) / region_length;
                            num_kmer += region_length - dist;
                            num_added += region_length - dist;
                        }
                    }
                    else
                    {
                        num_count += next.kmer_count();
                        num_kmer += next.num_kmer();
                        num_added += next.num_kmer();
                    }
                }
            }
        }
    }

    total_count += num_count;
    total_kmer += num_kmer;

    if (num_kmer == 0)
        return 0;
    else
        return num_count * 1.0 / num_kmer;
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

    //cout << "palindrome " << contigs.size() << endl;

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
        
        //cout << (reachable[search_node.node] == current.out_edges()) << " " << reachable[search_node.node] << " " << (int)current.out_edges() << endl;
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

void ContigGraph::Decomposite()
{
    int64_t last = 0;
    for (int i = 0; i < 100; ++i)
    {
        int64_t split = SplitBranches();
        //cout << split << " " << 2*vertices_.size() << endl;

        if (last == split)
            break;
        last = split;
    }
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

void ContigGraph::GetConsensus(deque<Sequence> &contigs)
{
    deque<deque<ContigGraphVertexAdaptor> > components;
    deque<string> component_strings;

    GetComponents(components, component_strings);
    for (unsigned i = 0; i < components.size(); ++i)
    {
        ContigGraphVertexAdaptor begin = GetBeginVertexAdaptor(components[i]);
        ContigGraphVertexAdaptor end = GetEndVertexAdaptor(components[i]);

        if (begin.is_null() || end.is_null() || !IsValid(components[i]))
        {
            for (unsigned j = 0; j < components[i].size(); ++j)
                contigs.push_back(components[i][j].contig());
        }
        else
        {
            ContigGraphPath path;
            FindLongestPath(components[i], path);
            Sequence contig;
            ContigInfo contig_info;
            path.Assemble(contig, contig_info);
            contigs.push_back(contig);
        }
    }
}

bool ContigGraph::FindPath(ContigGraphVertexAdaptor from, ContigGraphVertexAdaptor to, ContigGraphPath &path)
{
    path.clear();
    map<int, int> is_used;
    map<ContigGraphVertexAdaptor, ContigGraphVertexAdaptor> prev;
    deque<ContigGraphVertexAdaptor> qu;
    qu.push_back(from);
    prev[from] = ContigGraphVertexAdaptor(NULL);
    is_used[from.id()] = true;

    int time = 0;
    while (!qu.empty())
    {
        if (++time >= 100)
            break;

        if (prev.find(to) != prev.end())
            break;

        ContigGraphVertexAdaptor current = qu.front();
        qu.pop_front();

        deque<ContigGraphVertexAdaptor> neighbors;
        GetNeighbors(current, neighbors);
        for (unsigned i = 0; i < neighbors.size(); ++i)
        {
            ContigGraphVertexAdaptor next = neighbors[i];
            //if (prev.find(next) == prev.end())
            if (!is_used[next.id()])
            {
                is_used[next.id()] = true;
                prev[next] = current;
                qu.push_back(next);
            }
        }
    }

    if (prev.find(to) != prev.end())
    {
        deque<ContigGraphVertexAdaptor> tmp;
        tmp.push_back(to);
        while (!prev[tmp.back()].is_null())
            tmp.push_back(prev[tmp.back()]);
        reverse(tmp.begin(), tmp.end());
        for (unsigned i = 0; i < tmp.size(); ++i)
            path.Append(tmp[i], -kmer_size_ + 1);
        return true;
    }
    else
        return false;
}

void ContigGraph::GetContigs(deque<Sequence> &contigs, deque<ContigInfo> &contig_infos)
{
    contigs.resize(vertices_.size());
    contig_infos.resize(vertices_.size());

    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        contigs[i] = vertices_[i].contig();
        contig_infos[i] = vertices_[i].contig_info();
    }
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

bool ContigGraph::IsValid(deque<ContigGraphVertexAdaptor> &component)
{
    ContigGraphVertexAdaptor begin = GetBeginVertexAdaptor(component);
    ContigGraphVertexAdaptor end = GetEndVertexAdaptor(component);

    map<int, int> status;
    if (CycleDetect(begin, status))
        return false;

    if (status.size() != component.size())
        return false;

    status.clear();
    end.ReverseComplement();
    if (CycleDetect(end, status))
        return false;

    if (status.size() != component.size())
        return false;

    return true;
}

void ContigGraph::FindLongestPath(deque<ContigGraphVertexAdaptor> &component, ContigGraphPath &path)
{
    ContigGraphVertexAdaptor begin = GetBeginVertexAdaptor(component);
    ContigGraphVertexAdaptor end = GetEndVertexAdaptor(component);

    deque<ContigGraphVertexAdaptor> order;
    TopSort(component, order);

    map<ContigGraphVertexAdaptor, int> dist;
    map<ContigGraphVertexAdaptor, ContigGraphVertexAdaptor> prev;
    dist[begin] = 0;
    prev[begin] = ContigGraphVertexAdaptor(NULL);

    for (unsigned i = 0; i < order.size(); ++i)
    {
        ContigGraphVertexAdaptor current = order[i];
        for (int x = 0; x < 4; ++x)
        {
            if (current.out_edges()[x])
            {
                ContigGraphVertexAdaptor next = GetNeighbor(current, x);
                int tmp = dist[current] + (int)current.contig_size() - (int)kmer_size_ + 1;
                if (current.id() != next.id() && tmp > dist[next])
                {
                    dist[next] = tmp;
                    prev[next] = current;
                }
            }
        }
    }

    deque<ContigGraphVertexAdaptor> v;
    v.push_back(end);

    while (!prev[v.back()].is_null())
        v.push_back(prev[v.back()]);
    reverse(v.begin(), v.end());

    path.clear();
    for (unsigned i = 0; i < v.size(); ++i)
        path.Append(v[i], -(int)kmer_size_ + 1);
}

void ContigGraph::TopSort(deque<ContigGraphVertexAdaptor> &component, deque<ContigGraphVertexAdaptor> &order)
{
    ContigGraphVertexAdaptor begin = GetBeginVertexAdaptor(component);

    map<int, int> status;
    TopSortDFS(order, begin, status);
    reverse(order.begin(), order.end());
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

double ContigGraph::FindSimilarPath(ContigGraphVertexAdaptor target, ContigGraphVertexAdaptor start)
{
    if (start.status().IsDead()
            || target.begin_kmer(kmer_size_-1) != start.begin_kmer(kmer_size_-1)
            || target.in_edges() != start.in_edges())
        return 0;

    ContigGraphPath path;
    path.Append(start, 0);

    int time = 0;
    return FindSimilarPath(target, path, time);
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


