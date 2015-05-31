/**
 * @file contig_graph_vertex.h
 * @brief ContigGraphVertex Class and ContigGraphVertexAdaptor Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-16
 */

#ifndef __GRAPH_CONTIG_GRAPH_VERTEX_H_

#define __GRAPH_CONTIG_GRAPH_VERTEX_H_

#include <stdint.h>

#include <algorithm>
#include <string>

#include "lib_idba/kmer.h"
#include "lib_idba/bit_edges.h"
#include "lib_idba/contig_info.h"
#include "lib_idba/vertex_status.h"
#include "lib_idba/sequence.h"


/**
 * @brief It is the vertex class used in ContigGraph class.
 */
class ContigGraphVertex
{
public:
    explicit ContigGraphVertex(const Sequence &contig = Sequence(), const ContigInfo &contig_info = ContigInfo())
        : contig_(contig), contig_info_(contig_info) {}

    ContigGraphVertex(const ContigGraphVertex &x)
        : contig_(x.contig_), id_(x.id_), status_(x.status_), contig_info_(x.contig_info_) {}

    const ContigGraphVertex &operator =(const ContigGraphVertex &x)
    {
        if (this != &x)
        {
            contig_ = x.contig_;
            id_ = x.id_;
            contig_info_ = x.contig_info_;
        }
        return *this;
    }

    const Sequence &contig() const { return contig_; }
    void set_contig(const Sequence &contig) { contig_ = contig; }

    uint32_t contig_size() const { return contig_.size(); }
    uint32_t num_kmer() const { return contig_.size() - kmer_size() + 1; }

    const ContigInfo &contig_info() const { return contig_info_; }
    void set_contig_info(const ContigInfo &contig_info) { contig_info_ = contig_info; }

    uint64_t kmer_count() const { return contig_info_.kmer_count(); }
    void set_kmer_count(uint64_t kmer_count) { contig_info_.set_kmer_count(kmer_count); }

    uint32_t id() const { return id_; }
    void set_id(uint32_t id) { id_ = id; }

    uint32_t kmer_size() const { return contig_info_.kmer_size(); }
    void set_kmer_size(uint32_t kmer_size) { contig_info_.set_kmer_size(kmer_size); }

    VertexStatus &status() { return status_; }
    const VertexStatus &status() const { return status_; }

    BitEdges &in_edges() { return contig_info_.in_edges(); }
    const BitEdges &in_edges() const { return contig_info_.in_edges(); }

    BitEdges &out_edges() { return contig_info_.out_edges(); }
    const BitEdges &out_edges() const { return contig_info_.out_edges(); }

    IdbaKmer begin_kmer(int kmer_size) const { return contig_.GetIdbaKmer(0, kmer_size); }
    IdbaKmer end_kmer(int kmer_size) const { return contig_.GetIdbaKmer(contig_.size() - kmer_size, kmer_size); }

    double coverage() const { return 1.0 * contig_info_.kmer_count() / (contig_size() - kmer_size() + 1); }

    const SequenceCount &counts() const { return contig_info_.counts(); }
    void set_counts(const SequenceCount &counts) { contig_info_.set_counts(counts); }

    char get_base(uint32_t index) const { return contig_[index]; }

    SequenceCountUnitType get_count(uint32_t index) const { return contig_info_.counts()[index]; }

    void swap(ContigGraphVertex &x)
    { 
        if (this != &x)
        {
            contig_.swap(x.contig_); 
            std::swap(id_, x.id_); 
            status_.swap(x.status_); 
            contig_info_.swap(x.contig_info_);
        }
    }

    void clear() 
    { 
        contig_.clear(); 
        id_ = 0; 
        status_.clear(); 
        contig_info_.clear();
    }

private:
    Sequence contig_;

    uint32_t id_;
    VertexStatus status_;
    ContigInfo contig_info_;
};

/**
 * @brief It is a adaptor class used to access ContigGraphVertex. Becase a contig and its
 * reverse complement share the same vertex, using adaptor makes sure that modification to 
 * the vertex consistant.
 */
class ContigGraphVertexAdaptor
{
public:
    explicit ContigGraphVertexAdaptor(ContigGraphVertex *vertex = NULL, bool is_reverse = false)
    { vertex_ = vertex; is_reverse_ = is_reverse; }
    ContigGraphVertexAdaptor(const ContigGraphVertexAdaptor &x)
    { vertex_ = x.vertex_, is_reverse_ = x.is_reverse_; }

    const ContigGraphVertexAdaptor &operator =(const ContigGraphVertexAdaptor &x)
    { vertex_ = x.vertex_; is_reverse_ = x.is_reverse_; return *this; }
        
    bool operator <(const ContigGraphVertexAdaptor &x) const
    { return (vertex_ != x.vertex_) ? (vertex_ < x.vertex_) : (is_reverse_ < x.is_reverse_); }
    bool operator >(const ContigGraphVertexAdaptor &x) const
    { return (vertex_ != x.vertex_) ? (vertex_ > x.vertex_) : (is_reverse_ > x.is_reverse_); }

    bool operator ==(const ContigGraphVertexAdaptor &x) const
    { return vertex_ == x.vertex_ && is_reverse_ == x.is_reverse_; }
    bool operator !=(const ContigGraphVertexAdaptor &x) const
    { return vertex_ != x.vertex_ || is_reverse_ != x.is_reverse_; }

    const ContigGraphVertexAdaptor &ReverseComplement() { is_reverse_ = !is_reverse_; return *this; }

    Sequence contig() const
    {
        Sequence contig = vertex_->contig();
        return !is_reverse_ ? contig : contig.ReverseComplement();
    }

    uint32_t contig_size() const { return vertex_->contig().size(); }
    uint32_t num_kmer() const { return vertex_->num_kmer(); }

    void set_vertex(ContigGraphVertex *vertex, bool is_reverse)
    { vertex_ = vertex; is_reverse_ = is_reverse; }

    ContigInfo contig_info() const
    {
        ContigInfo contig_info = vertex_->contig_info();
        return (!is_reverse_ ? contig_info : contig_info.ReverseComplement()); 
    }

    uint64_t kmer_size() const { return vertex_->kmer_size(); }
    void set_kmer_size(uint64_t kmer_size) { vertex_->set_kmer_size(kmer_size); }

    uint64_t kmer_count() const { return vertex_->kmer_count(); }
    void set_kmer_count(uint64_t kmer_count) { vertex_->set_kmer_count(kmer_count); }

    uint32_t id() const { return vertex_->id(); }
    void set_id(uint32_t id) { vertex_->set_id(id); }

    VertexStatus &status() { return vertex_->status(); }
    const VertexStatus &status() const { return vertex_->status(); }

    BitEdges &in_edges() { return !is_reverse_ ? vertex_->in_edges() : vertex_->out_edges(); }
    const BitEdges &in_edges() const { return !is_reverse_ ? vertex_->in_edges() : vertex_->out_edges(); }

    BitEdges &out_edges() { return !is_reverse_ ? vertex_->out_edges() : vertex_->in_edges(); }
    const BitEdges &out_edges() const { return !is_reverse_ ? vertex_->out_edges() : vertex_->in_edges(); }

    SequenceCount counts()
    {
        if (!is_reverse_) return vertex_->counts();
        else { SequenceCount counts = vertex_->counts(); std::reverse(counts.begin(), counts.end()); return counts; }
    }

    char get_base(uint32_t index) const
    { return (!is_reverse_) ? vertex_->get_base(index) : 3 - vertex_->get_base(contig_size() - 1 - index); }

    SequenceCountUnitType get_count(uint32_t index) const 
    { return (!is_reverse_) ? vertex_->get_count(index) : vertex_->get_count(vertex_->counts().size() - 1 - index); }

    IdbaKmer begin_kmer(int kmer_size) const 
    { return !is_reverse_ ? vertex_->begin_kmer(kmer_size) : vertex_->end_kmer(kmer_size).ReverseComplement(); }

    IdbaKmer end_kmer(int kmer_size) const
    { return !is_reverse_ ? vertex_->end_kmer(kmer_size) : vertex_->begin_kmer(kmer_size).ReverseComplement(); }

    double coverage() const 
    { return vertex_->coverage(); }

    bool is_reverse() const { return is_reverse_; }

    void swap(ContigGraphVertexAdaptor &x)
    { 
        if (this != &x)
        {
            std::swap(vertex_, x.vertex_); 
            std::swap(is_reverse_, x.is_reverse_); 
        }
    }

    bool is_null() const { return vertex_ == NULL; }
    void clear() { vertex_->clear(); }

private:
    ContigGraphVertex *vertex_;
    bool is_reverse_;
};

namespace std
{
template <> inline void swap(ContigGraphVertex &x, ContigGraphVertex &y) { x.swap(y); }
template <> inline void swap(ContigGraphVertexAdaptor &x, ContigGraphVertexAdaptor &y) { x.swap(y); }
}

#endif

