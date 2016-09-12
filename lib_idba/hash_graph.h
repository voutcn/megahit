/**
 * @file hash_graph.h
 * @brief HashGraph Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-05
 */

#ifndef __GRAPH_HASH_GRAPH_H_

#define __GRAPH_HASH_GRAPH_H_

#include <deque>
#include <istream>
#include <ostream>

#include "bit_operation.h"
#include "histgram.h"
#include "lib_idba/hash_table.h"
#include "lib_idba/kmer.h"
#include "lib_idba/contig_info.h"
#include "lib_idba/hash_graph_vertex.h"
#include "lib_idba/hash_graph_path.h"
#include "lib_idba/sequence.h"

class IdbaKmer;
class Sequence;
class ShortSequence;
class CompactSequence;


/**
 * @brief It is a hash table based de Bruijn graph implementation.
 */
class HashGraph
{
    class RefreshVerticesFunc;
    class RefreshEdgesFunc;

public:
    friend std::istream &operator >>(std::istream &is, HashGraph &hash_graph);
    friend std::ostream &operator <<(std::ostream &os, HashGraph &hash_graph);

    typedef HashTableST<HashGraphVertex, IdbaKmer> vertex_table_type;
    typedef vertex_table_type::iterator iterator;

    explicit HashGraph(uint32_t kmer_size = 0) { set_kmer_size(kmer_size); num_edges_ = 0; }
    ~HashGraph() {}

    iterator begin() { return vertex_table_.begin(); }
    iterator end() { return vertex_table_.end(); }

    HashGraphVertex *InsertVertex(const IdbaKmer &kmer, int count = 1)
    { 
        IdbaKmer key = kmer.unique_format();
        HashGraphVertex &vertex = vertex_table_.find_or_insert(HashGraphVertex(key));
        vertex.count() += count;
        return &vertex;
    }

    HashGraphVertex *InsertVertex(const HashGraphVertex &vertex)
    { return &vertex_table_.find_or_insert(vertex); }

    HashGraphVertex *FindVertex(const IdbaKmer &kmer)
    {
        IdbaKmer key = kmer.unique_format();
        vertex_table_type::iterator p = vertex_table_.find(key);
        return (p != vertex_table_.end()) ? &*p : NULL;
    }

    const HashGraphVertex *FindVertex(const IdbaKmer &kmer) const
    {
        IdbaKmer key = kmer.unique_format();
        vertex_table_type::const_iterator p = vertex_table_.find(key);
        return (p != vertex_table_.end()) ? &*p : NULL;
    }

    HashGraphVertexAdaptor FindVertexAdaptor(const IdbaKmer &kmer)
    { 
        IdbaKmer key = kmer.unique_format();
        vertex_table_type::iterator p = vertex_table_.find(key);
        return ((p != vertex_table_.end()) ? 
            HashGraphVertexAdaptor(&*p, kmer != key) : HashGraphVertexAdaptor(NULL));
    }

    HashGraphVertexAdaptor GetNeighbor(const HashGraphVertexAdaptor &current, int x)
    {
        IdbaKmer kmer = current.kmer();
        kmer.ShiftAppend(x);
        return FindVertexAdaptor(kmer);
    }

    int64_t InsertKmers(const Sequence &seq) { return InsertKmersWithPrefix(seq, 0, 0); }
    int64_t InsertKmersWithPrefix(const Sequence &seq, uint64_t prefix, uint64_t umask);
    int64_t InsertUncountKmers(const Sequence &seq);
    int64_t InsertInternalKmers(const Sequence &seq, int min_count = 0);
    int64_t InsertEdges(const Sequence &seq);
    int64_t InsertExistKmers(const Sequence &seq);

    int64_t RemoveKmers(const Sequence &seq);
    void RemoveEdge(HashGraphVertexAdaptor &node, int x)
    {
        node.out_edges().Remove(x);
        IdbaKmer kmer = node.kmer();
        kmer.ShiftAppend(x);
        HashGraphVertexAdaptor next = FindVertexAdaptor(kmer);
        if (!next.is_null())
            next.in_edges().Remove(3 - node.kmer()[0]);
    }

    void AddEdge(HashGraphVertexAdaptor &node, int x)
    {
        node.out_edges().Add(x);
        IdbaKmer kmer = node.kmer();
        kmer.ShiftAppend(x);
        HashGraphVertexAdaptor next = FindVertexAdaptor(kmer);
        if (!next.is_null())
            next.in_edges().Add(3 - node.kmer()[0]);
    }

    void BackupEdges()
    { BackupEdgesFunc func; vertex_table_.for_each(func); }

    void RestoreAndMergeEdges()
    { RestoreAndMergeEdgesFunc func; vertex_table_.for_each(func); }

    void AddAllEdges()
    { AddAllEdgesFunc func; vertex_table_.for_each(func); RefreshEdges(); }
    void ClearEdges()
    { ClearEdgesFunc func; vertex_table_.for_each(func); }
    void ClearStatus()
    { ClearStatusFunc func; vertex_table_.for_each(func); }
    void ClearCount()
    { ClearCountFunc func; vertex_table_.for_each(func); }

    void SetCountCap(int cap)
    { SetCountCapFunc func(cap); vertex_table_.for_each(func); }

    void Refresh(int min_count = 0) 
    { RefreshVertices(min_count); RefreshEdges(); }
    int64_t RefreshVertices(int min_count = 0)
    { RefreshVerticesFunc func(min_count); return vertex_table_.remove_if(func); }
    void RefreshEdges()
    { RefreshEdgesFunc func(this); vertex_table_.for_each(func); num_edges_ = func.num_edges(); }

    int64_t ErodeEnd(int min_cover);
    int64_t Trim(int min_length);
    int64_t RemoveDeadEnd(int min_length);
    int64_t RemoveLowCoverage(double min_cover, int min_contig = (1 << 20));
    int64_t RemoveBubble();

    int64_t Assemble(std::deque<Sequence> &contigs);
    int64_t Assemble(std::deque<Sequence> &contigs, std::deque<ContigInfo> &contig_infos); 
    
//    int64_t TrimSequentially(int min_length);
//    int64_t RemoveDeadEndSequentially(int min_length);
//    int64_t RemoveLowCoverageSequentially(double min_cover);
//    int64_t AssembleSequentially(std::deque<Sequence> &contigs);
//    int64_t AssembleSequentially(std::deque<Sequence> &contigs, std::deque<ContigInfo> &contig_infos); 

    void reserve(uint64_t capacity) { vertex_table_.reserve(capacity); }

    uint32_t kmer_size() const { return kmer_size_; }
    void set_kmer_size(uint32_t kmer_size) { kmer_size_ = kmer_size; }

    Histgram<int> coverage_histgram() const
    {
        CoverageHistgramFunc func;
        vertex_table_.for_each(func);
        return func.histgram();
    }

    void swap(HashGraph &hash_graph)
    {
        if (this != &hash_graph)
        {
            vertex_table_.swap(hash_graph.vertex_table_);
            std::swap(kmer_size_, hash_graph.kmer_size_);
            std::swap(num_edges_, hash_graph.num_edges_);
        }
    }

    uint64_t num_bucket() const { return vertex_table_.bucket_count(); }
    uint64_t num_vertices() const { return vertex_table_.size(); }
    uint64_t num_edges() const { return num_edges_; }
    void clear() { vertex_table_.clear(); num_edges_ = 0; }

private:
#if __cplusplus >= 201103L
    HashGraph(const HashGraph &) = delete;
    const HashGraph &operator =(const HashGraph &) = delete;
#else
    HashGraph(const HashGraph &);
    const HashGraph &operator =(const HashGraph &);
#endif

    bool GetNextVertexAdaptor(const HashGraphVertexAdaptor &current, HashGraphVertexAdaptor &next)
    {
        if (current.out_edges().size() != 1)
            return false;

        IdbaKmer kmer = current.kmer();
        kmer.ShiftAppend(bit_operation::BitToIndex(current.out_edges()));
        next = FindVertexAdaptor(kmer);

        return !kmer.IsPalindrome() && next.in_edges().size() == 1;
    }

    bool IsLoop(const Sequence &contig, HashGraphVertexAdaptor &next)
    {
        IdbaKmer kmer = next.kmer();
        IdbaKmer rev_comp = kmer;
        rev_comp.ReverseComplement();

        return contig.GetIdbaKmer(0, kmer_size_) == kmer;
    }

    bool IsPalindromeLoop(const Sequence &contig, HashGraphVertexAdaptor &next)
    {
        IdbaKmer kmer = next.kmer();
        IdbaKmer rev_comp = kmer;
        rev_comp.ReverseComplement();

        return contig.GetIdbaKmer(contig.size() - kmer_size_, kmer_size_) == rev_comp;   
    }

    class BackupEdgesFunc
    {
    public:
        BackupEdgesFunc() {}

        void operator() (HashGraphVertex &vertex)
        {
            vertex.in_edges() = (vertex.in_edges() << 4) | (vertex.in_edges() & 15);
            vertex.out_edges() = (vertex.out_edges() << 4) | (vertex.out_edges() & 15);
        }
    };

    class RestoreAndMergeEdgesFunc
    {
    public:
        RestoreAndMergeEdgesFunc() {}

        void operator() (HashGraphVertex &vertex)
        {
            vertex.in_edges() = ((unsigned)vertex.in_edges() >> 4) | (vertex.in_edges() & 15);
            vertex.out_edges() = ((unsigned)vertex.out_edges() >> 4) | (vertex.out_edges() & 15);
        }
    };

    class AddAllEdgesFunc
    {
    public:
        AddAllEdgesFunc() {}

        void operator ()(HashGraphVertex &vertex)
        { vertex.in_edges() = 15; vertex.out_edges() = 15; }
    };

    class ClearEdgesFunc
    {
    public:
        ClearEdgesFunc() {}

        void operator ()(HashGraphVertex &vertex)
        { vertex.in_edges() = 0; vertex.out_edges() = 0; }
    };

    class ClearStatusFunc
    {
    public:
        ClearStatusFunc() {}

        void operator ()(HashGraphVertex &vertex)
        { vertex.status().clear(); }
    };

    class ClearCountFunc
    {
    public:
        ClearCountFunc() {}

        void operator ()(HashGraphVertex &vertex)
        { vertex.count() = 0; }
    };

    class SetCountCapFunc
    {
    public:
        SetCountCapFunc(int cap): cap_(cap) { }

        void operator ()(HashGraphVertex &vertex)
        { if (vertex.count() > cap_) vertex.count() = cap_; }

    private:
        int cap_;
    };

    class RefreshVerticesFunc
    {
    public:
        explicit RefreshVerticesFunc(int min_count) : min_count_(min_count) {}

        bool operator ()(HashGraphVertex &vertex) const
        {
            if (vertex.count() < min_count_ || vertex.status().IsDead())
                return true;
            return false;
        }

    private:
        int min_count_;
    };

    class RefreshEdgesFunc
    {
    public:
        explicit RefreshEdgesFunc(HashGraph *hash_graph) { hash_graph_ = hash_graph; total_degree_ = 0; }

        void operator ()(HashGraphVertex &vertex)
        {
            HashGraphVertexAdaptor adaptor(&vertex);
            for (int strand = 0; strand < 2; ++strand)
            {
                IdbaKmer kmer = adaptor.kmer();
                for (int i = 0; i < 4; ++i)
                {
                    if (adaptor.out_edges()[i])
                    {
                        IdbaKmer next = kmer;
                        next.ShiftAppend(i);
                        if (hash_graph_->FindVertex(next) == NULL)
                            adaptor.out_edges().Remove(i);
                        else
                            total_degree_ += 1;
                    }
                }
                adaptor.ReverseComplement();
            }

            if ((vertex.kmer().size() & 1) == 0)
                vertex.FixPalindromeEdges();
        }

        uint64_t num_edges() { return total_degree_ / 2; }

    private:
        HashGraph *hash_graph_;
        uint64_t total_degree_;
    };

    class ErodeFunc
    {
    public:
        ErodeFunc(HashGraph *hash_graph, int min_cover)
        { hash_graph_ = hash_graph; min_cover_ = min_cover; }

        void operator ()(HashGraphVertex &vertex);

    private:
        HashGraph *hash_graph_;
        int min_cover_;
    };

    class TrimFunc
    {
    public:
        TrimFunc(HashGraph *hash_graph, int min_length)
        { hash_graph_ = hash_graph; min_length_ = min_length; }

        void operator ()(HashGraphVertex &vertex);

    private:
        HashGraph *hash_graph_;
        int min_length_;
    };

    class BubbleFunc
    {
    public:
        BubbleFunc(HashGraph *hash_graph)
        { hash_graph_ = hash_graph; }
        ~BubbleFunc()
        { }

        void operator ()(HashGraphVertex &vertex);

        std::deque<HashGraphVertexAdaptor> &candidates() { return candidates_; }

    private:
        HashGraph *hash_graph_;
        std::deque<HashGraphVertexAdaptor> candidates_;
    };

    class AssembleFunc
    {
    public:
        AssembleFunc(HashGraph *hash_graph)
            : hash_graph_(hash_graph)
        { }
        ~AssembleFunc()
        { }
            
        void operator ()(HashGraphVertex &vertex);

        std::deque<Sequence> &contigs() { return contigs_; }
        std::deque<ContigInfo> &contig_infos() { return contig_infos_; }

    private:
        HashGraph *hash_graph_;
        std::deque<Sequence> contigs_;
        std::deque<ContigInfo> contig_infos_;
    };

    class CoverageHistgramFunc
    {
    public:
        void operator ()(HashGraphVertex &vertex)
        { histgram_.insert(vertex.count()); }

        const Histgram<int> &histgram() { return histgram_; }

    private:
        Histgram<int> histgram_;
    };

    HashTableST<HashGraphVertex, IdbaKmer> vertex_table_;
    uint32_t kmer_size_;
    uint64_t num_edges_;
};

inline std::istream &operator >>(std::istream &is, HashGraph &hash_graph)
{ return is >> hash_graph.vertex_table_; }

inline std::ostream &operator <<(std::ostream &os, HashGraph &hash_graph)
{ os << hash_graph.vertex_table_; hash_graph.RefreshEdges(); return os; }

namespace std
{
inline void swap(HashGraph &x, HashGraph &y)
{ x.swap(y); }
}

#endif

