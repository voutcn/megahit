/**
 * @file contig_graph_path.h
 * @brief ContigGraph Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-16
 */

#ifndef __GRAPH_CONTIG_GRAPH_PATH_H_

#define __GRAPH_CONTIG_GRAPH_PATH_H_

#include <stdint.h>

#include <deque>

#include "lib_idba/contig_builder.h"
#include "lib_idba/contig_graph_vertex.h"
#include "lib_idba/contig_info.h"


/**
 * @brief It is a path of contigs in ContigGraph.
 */
class ContigGraphPath
{
public:
    ContigGraphPath() {}
    ContigGraphPath(const ContigGraphPath &path)
        : vertices_(path.vertices_), distances_(path.distances_) {}

    const ContigGraphPath &operator =(const ContigGraphPath &path)
    { vertices_ = path.vertices_; distances_ = path.distances_; return *this; }

    bool operator <(const ContigGraphPath &path) const
    {
		for (unsigned i = 0; i < num_nodes() && i < path.num_nodes(); ++i)
		{
			if ((*this)[i] != path[i])
				return (*this)[i] < path[i];
		}

		return num_nodes() < path.num_nodes();
    }

    ContigGraphVertexAdaptor &operator [](uint32_t index) 
    { return vertices_[index]; }
    const ContigGraphVertexAdaptor &operator [](uint32_t index) const 
    { return vertices_[index]; }

    void Append(const ContigGraphVertexAdaptor &vertex, int d)
    {
        vertices_.push_back(vertex);
        if (vertices_.size() > 1)
            distances_.push_back(d);
    }

    void Append(const ContigGraphPath &path, int d)
    {
        for (unsigned i = 0; i < path.num_nodes(); ++i)
        {
            if (i == 0)
                Append(path[i], d);
            else
                Append(path[i], path.distances()[i-1]);
        }
    }

    void Pop()
    {
        vertices_.pop_back();
        if (!distances_.empty())
            distances_.pop_back();
    }

    const ContigGraphPath &ReverseComplement()
    {
        std::reverse(vertices_.begin(), vertices_.end());
        for (unsigned i = 0; i < vertices_.size(); ++i)
            vertices_[i].ReverseComplement();
        std::reverse(distances_.begin(), distances_.end());
        return *this;
    }

    void Assemble(Sequence &contig, ContigInfo &contig_info)
    {
        ContigBuilder contig_builder;

        if (vertices_.size() > 0)
        {
            contig_builder.Append(vertices_[0], 0);
            for (unsigned i = 1; i < vertices_.size(); ++i)
                contig_builder.Append(vertices_[i], distances_[i-1]);
        }

        contig = contig_builder.contig();
        contig_info = contig_builder.contig_info();
    }

    void swap(ContigGraphPath &path) 
    {
        if (this != &path)
        {
            vertices_.swap(path.vertices_); 
            distances_.swap(path.distances_); 
        }
    }

    ContigGraphVertexAdaptor &front() { return vertices_.front(); }
    const ContigGraphVertexAdaptor &front() const { return vertices_.front(); }

    ContigGraphVertexAdaptor &back() { return vertices_.back(); }
    const ContigGraphVertexAdaptor &back() const { return vertices_.back(); }

    uint64_t kmer_count() const
    {
        uint64_t sum = 0;
        for (unsigned i = 0; i < vertices_.size(); ++i)
            sum += vertices_[i].kmer_count();
        return sum;
    }

    uint32_t size() const
    {
        uint32_t size = 0;
        for (unsigned i = 0; i < vertices_.size(); ++i)
            size += vertices_[i].contig_size();
        for (unsigned i = 0; i < distances_.size(); ++i)
            size += distances_[i];
        return size;
    }

    uint32_t internal_size(int kmer_size) const
    {
        if (vertices_.size() <= 1)
            return vertices_.size();

        uint32_t size = kmer_size + 1;
        for (unsigned i = 1; i+1 < vertices_.size(); ++i)
            size += vertices_[i].contig_size();
        for (unsigned i = 0; i < distances_.size(); ++i)
            size += distances_[i];
        return size;
    }

    uint32_t num_nodes() const { return vertices_.size(); }
    void clear() { vertices_.clear(); distances_.clear(); }

    std::deque<int> &distances() { return distances_; }
    const std::deque<int> &distances() const { return distances_; }

private:
    std::deque<ContigGraphVertexAdaptor> vertices_;
    std::deque<int> distances_; 
};

namespace std
{
template <> inline void swap(ContigGraphPath &x, ContigGraphPath &y) { x.swap(y); }
}

#endif

