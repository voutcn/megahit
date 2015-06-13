/**
 * @file hash_graph_vertex.h
 * @brief HashGraphVertex Class and HashGraphVertexAdaptor Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-05
 */

#ifndef __GRAPH_HASH_GRAPH_VERTEX_H_

#define __GRAPH_HASH_GRAPH_VERTEX_H_

#include <algorithm>

#include "bit_operation.h"
#include "lib_idba/kmer.h"
#include "lib_idba/bit_edges.h"
#include "lib_idba/vertex_status.h"


/**
 * @brief It is the vertex class used in HashGraph.
 */
class HashGraphVertex
{
public:
    explicit HashGraphVertex(const IdbaKmer &kmer = IdbaKmer()): kmer_(kmer), count_(0) {}
    HashGraphVertex(const HashGraphVertex &x)
        : kmer_(x.kmer_), count_(x.count_), status_(x.status_), in_edges_(x.in_edges_), out_edges_(x.out_edges_) {}

    const HashGraphVertex &operator =(const HashGraphVertex &x)
    { 
        kmer_ = x.kmer_; 
        count_ = x.count_;
        status_ = x.status_; 
        in_edges_ = x.in_edges_;
        out_edges_ = x.out_edges_; 
        return *this; 
    }

    void FixPalindromeEdges()
    { if (kmer_.IsPalindrome()) out_edges_ = in_edges_ = (in_edges_ | out_edges_); }

    const IdbaKmer &key() const { return kmer_; }
    void set_key(const IdbaKmer &key) { kmer_ = key; }

    const IdbaKmer &kmer() const { return kmer_; }
    void set_kmer(const IdbaKmer &kmer) { kmer_ = kmer; }

    int32_t &count() { return count_; }
    const int32_t &count() const { return count_; }

    VertexStatus &status() { return status_; }
    const VertexStatus &status() const { return status_; }

    BitEdges &in_edges() { return in_edges_; }
    const BitEdges &in_edges() const { return in_edges_; }

    BitEdges &out_edges() { return out_edges_; }
    const BitEdges &out_edges() const { return out_edges_; }

    void swap(HashGraphVertex &x)
    { 
        if (this != &x)
        {
            kmer_.swap(x.kmer_); 
            std::swap(count_, x.count_); 
            status_.swap(x.status_); 
            in_edges_.swap(x.in_edges_); 
            out_edges_.swap(x.out_edges_); 
        }
    }

    uint32_t kmer_size() const { return kmer_.size(); }

    void clear() { in_edges_.clear(); out_edges_.clear() ; status_.clear(); count_ = 0; }

private:
    IdbaKmer kmer_;

    int32_t count_;
    VertexStatus status_;
    BitEdges in_edges_;
    BitEdges out_edges_;
};


/**
 * @brief It is adaptor class used for accessing HashGraphVertex. Because 
 * a k-mer and its reverse complemtn share the same vertex, using adaptor
 * makes sure the access to vertex consistant.
 */
class HashGraphVertexAdaptor
{
public:
    explicit HashGraphVertexAdaptor(HashGraphVertex *vertex = NULL, bool is_reverse = false)
    { vertex_ = vertex; is_reverse_ = is_reverse; }
    HashGraphVertexAdaptor(const HashGraphVertexAdaptor &x)
    { vertex_ = x.vertex_; is_reverse_ = x.is_reverse_; }

    const HashGraphVertexAdaptor &operator =(const HashGraphVertexAdaptor &x) 
    { vertex_ = x.vertex_; is_reverse_ = x.is_reverse_; return *this; }

    bool operator <(const HashGraphVertexAdaptor &x) const
    { return (vertex_ != x.vertex_) ? (vertex_ < x.vertex_) : (is_reverse_ < x.is_reverse_); }
    bool operator >(const HashGraphVertexAdaptor &x) const
    { return (vertex_ != x.vertex_) ? (vertex_ > x.vertex_) : (is_reverse_ > x.is_reverse_); }

    bool operator ==(const HashGraphVertexAdaptor &x) const
    { return vertex_ == x.vertex_ && is_reverse_ == x.is_reverse_; }
    bool operator !=(const HashGraphVertexAdaptor &x) const
    { return vertex_ != x.vertex_ || is_reverse_ != x.is_reverse_; }

    const HashGraphVertexAdaptor &ReverseComplement() { is_reverse_ = !is_reverse_; return *this; }

    IdbaKmer kmer() const 
    { 
        IdbaKmer kmer = vertex_->kmer(); 
        return !is_reverse_ ? kmer : kmer.ReverseComplement(); 
    }

    HashGraphVertex &vertex() { return *vertex_; }
    const HashGraphVertex &vertex() const { return *vertex_; }
    void set_vertex(HashGraphVertex *vertex, bool is_reverse = false)
    { vertex_ = vertex; is_reverse_ = is_reverse; }

    int32_t &count() { return vertex_->count(); }
    const int32_t &count() const { return vertex_->count(); }

    VertexStatus &status() { return vertex_->status(); }
    const VertexStatus &status() const { return vertex_->status(); }

    BitEdges &in_edges() { return !is_reverse_ ? vertex_->in_edges() : vertex_->out_edges(); }
    const BitEdges &in_edges() const { return !is_reverse_ ? vertex_->in_edges() : vertex_->out_edges(); }

    BitEdges &out_edges() { return !is_reverse_ ? vertex_->out_edges() : vertex_->in_edges(); }
    const BitEdges &out_edges() const { return !is_reverse_ ? vertex_->out_edges() : vertex_->in_edges(); }

    void swap(HashGraphVertexAdaptor &x)
    { 
        if (this != &x)
        {
            std::swap(vertex_, x.vertex_); 
            std::swap(is_reverse_, x.is_reverse_); 
        }
    }

    bool is_null() const { return vertex_ == NULL; }

    uint32_t kmer_size() const { return vertex_->kmer_size(); }

    void clear() { vertex_->clear(); }

private:
    HashGraphVertex *vertex_;
    bool is_reverse_;
};

namespace std
{
template <> inline void swap(HashGraphVertex &x, HashGraphVertex &y) { x.swap(y); }
template <> inline void swap(HashGraphVertexAdaptor &x, HashGraphVertexAdaptor &y) { x.swap(y); }
}

#endif


