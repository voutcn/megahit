/**
 * @file hash_graph_path.h
 * @brief HashGraphPath Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.4
 * @date 2011-09-21
 */

#ifndef __GRAPH_HASH_GRAPH_PATH_H_

#define __GRAPH_HASH_GRAPH_PATH_H_

#include "lib_idba/kmer.h"
#include "lib_idba/hash_graph.h"

#include <deque>


/**
 * @brief It is a path of k-mers in de Bruijn graph (HashGraph).
 */
class HashGraphPath
{
public:
    HashGraphPath() {}
    HashGraphPath(const HashGraphPath &path)
        : vertices_(path.vertices_) {}

    const HashGraphPath &operator =(const HashGraphPath &path)
    { vertices_ = path.vertices_; return *this; }

    HashGraphVertexAdaptor &operator [](uint32_t index)
    { return vertices_[index]; }

    const HashGraphVertexAdaptor &operator [](uint32_t index) const
    { return vertices_[index]; }

    void Append(const HashGraphVertexAdaptor &vertex)
    { vertices_.push_back(vertex); }

    void Pop()
    { vertices_.pop_back(); }

    const HashGraphPath &ReverseComplement()
    {
        std::reverse(vertices_.begin(), vertices_.end());
        for (unsigned i = 0; i < vertices_.size(); ++i)
            vertices_[i].ReverseComplement();
        return *this;
    }

    bool IsSimplePath() const
    {
        for (unsigned i = 1; i+1 < vertices_.size(); ++i)
        {
            if (vertices_[i].out_edges().size() != 1)
                return false;

            if (vertices_[i].in_edges().size() != 1)
                return false;
        }

        return true;
    }

    void swap(HashGraphPath &path)
    {
        if (this != &path)
            vertices_.swap(path.vertices_);
    }

    uint64_t kmer_count()
    {
        uint64_t sum = 0;
        for (unsigned i = 0; i < vertices_.size(); ++i)
            sum += vertices_[i].count();
        return sum;
    }

    HashGraphVertexAdaptor &front() { return vertices_.front(); }
    const HashGraphVertexAdaptor &front() const { return vertices_.front(); }

    HashGraphVertexAdaptor &back() { return vertices_.back(); }
    const HashGraphVertexAdaptor &back() const { return vertices_.back(); }

    uint32_t size() const
    {
        if (vertices_.empty())
            return 0;
        else
            return vertices_[0].kmer().size() + vertices_.size() - 1;
    }

    uint32_t num_nodes() const
    { return vertices_.size(); }

    void clear()
    { vertices_.clear(); }

private:
    std::deque<HashGraphVertexAdaptor> vertices_;
};

#endif

