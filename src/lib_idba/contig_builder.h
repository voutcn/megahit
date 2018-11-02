/**
 * @file contig_builder.h
 * @brief Contig Build Class which builds contig and related contig info.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.9
 * @date 2011-12-27
 */

#ifndef __GRAPH_CONTIG_BUILDER_H_

#define __GRAPH_CONTIG_BUILDER_H_

#include "lib_idba/contig_graph_vertex.h"
#include "lib_idba/contig_info.h"
#include "lib_idba/hash_graph_vertex.h"
#include "lib_idba/sequence.h"


/**
 * @brief It is a builder class for building contigs.
 */
class ContigBuilder
{
public:
    ContigBuilder()
    {}

    explicit ContigBuilder(HashGraphVertexAdaptor x)
    { Append(x); }

    explicit ContigBuilder(ContigGraphVertexAdaptor x)
    { Append(x, 0); }

    void Append(HashGraphVertexAdaptor x)
    {
        if (contig_.size() == 0)
        {
            contig_.Assign(x.kmer());
            contig_info_.in_edges_ = x.in_edges();
            contig_info_.out_edges_ = x.out_edges();
            contig_info_.kmer_size_ = x.kmer().size();
            contig_info_.kmer_count_ = x.count();
            contig_info_.counts_.resize(1);
            contig_info_.counts_[0] = x.count();
        }
        else
        {
            contig_ += x.kmer()[x.kmer().size() - 1];
            contig_info_.out_edges_ = x.out_edges();
            contig_info_.kmer_count_ += x.count();
            contig_info_.counts_ += x.count();
        }
    }

    void Append(ContigGraphVertexAdaptor x, int d)
    {
        if (contig_.size() == 0)
        {
            contig_ = x.contig();
            contig_info_.in_edges_ = x.in_edges();
            contig_info_.out_edges_ = x.out_edges();
            contig_info_.kmer_size_ = x.kmer_size(); 
            contig_info_.kmer_count_ = x.kmer_count();
            contig_info_.counts_ = x.counts();
        }
        else
        {
            if (d <= 0)
            {
                contig_.Append(x.contig(), std::min(-d, (int)x.contig_size()));
                contig_info_.out_edges_ = x.out_edges();
                contig_info_.kmer_count_ += x.kmer_count();
                SequenceCount counts = x.counts();
                contig_info_.counts_ += counts.substr(std::min(-d - contig_info_.kmer_size_ + 1, (int)counts.size()));
            }
            else
            {
                contig_.Append(d, 4);
                contig_.Append(x.contig());
                contig_info_.out_edges_ = x.out_edges();
                contig_info_.kmer_count_ += x.kmer_count();
                contig_info_.counts_.append(d, 0);
                contig_info_.counts_ += x.counts();
            }
        }
    }

    const ContigBuilder &ReverseComplement()
    { contig_.ReverseComplement(); contig_info_.ReverseComplement(); return *this; }

    const Sequence &contig() const 
    { return contig_; }

    const ContigInfo &contig_info() const
    { return contig_info_; }

    void clear()
    { contig_.clear(); contig_info_.clear(); }

private:
    Sequence contig_;
    ContigInfo contig_info_;
};

#endif

