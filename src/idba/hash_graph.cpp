/**
 * @file hash_graph.cpp
 * @brief
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-05
 */

#include "idba/hash_graph.h"

#include <algorithm>
#include <deque>
#include <stdexcept>

#include "bit_operation.h"
#include "idba/contig_builder.h"
#include "idba/contig_info.h"
#include "idba/hash_graph_vertex.h"
#include "idba/kmer.h"
#include "idba/sequence.h"
#include "utils/histgram.h"

using namespace std;

#include <iostream>

int64_t HashGraph::InsertKmers(const Sequence &seq) {
  if (seq.size() < kmer_size_) return 0;

  IdbaKmer kmer(kmer_size_);
  int length = 0;
  int64_t num_kmers = 0;
  for (uint64_t i = 0; i < seq.size(); ++i) {
    kmer.ShiftAppend(seq[i]);
    length = (seq[i] < 4) ? length + 1 : 0;

    if (length < (int)kmer_size_) continue;

    IdbaKmer key = kmer.unique_format();
    HashGraphVertex &vertex =
        vertex_table_.find_or_insert(HashGraphVertex(key));
    vertex.count() += 1;
    HashGraphVertexAdaptor adaptor(&vertex, kmer != key);

    if (length > (int)kmer_size_ && seq[i - kmer_size_] < 4)
      adaptor.in_edges().Add(3 - seq[i - kmer_size_]);
    if (i + 1 < seq.size() && seq[i + 1] < 4)
      adaptor.out_edges().Add(seq[i + 1]);

    ++num_kmers;
  }

  return num_kmers;
}

int64_t HashGraph::InsertUncountKmers(const Sequence &seq) {
  if (seq.size() < kmer_size_) return 0;

  IdbaKmer kmer(kmer_size_);
  int length = 0;
  int64_t num_kmers = 0;
  for (uint64_t i = 0; i < seq.size(); ++i) {
    kmer.ShiftAppend(seq[i]);
    length = (seq[i] < 4) ? length + 1 : 0;

    if (length < (int)kmer_size_) continue;

    IdbaKmer key = kmer.unique_format();

    HashGraphVertex &vertex =
        vertex_table_.find_or_insert(HashGraphVertex(key));
    HashGraphVertexAdaptor adaptor(&vertex, kmer != key);

    if (length > (int)kmer_size_ && seq[i - kmer_size_] < 4)
      adaptor.in_edges().Add(3 - seq[i - kmer_size_]);
    if (i + 1 < seq.size() && seq[i + 1] < 4)
      adaptor.out_edges().Add(seq[i + 1]);

    ++num_kmers;
  }

  return num_kmers;
}

int64_t HashGraph::Assemble(std::deque<Sequence> &contigs,
                            std::deque<ContigInfo> &contig_infos) {
  contigs.clear();
  contig_infos.clear();
  AssembleFunc func(this);
  vertex_table_.for_each(func);
  contigs.swap(func.contigs());
  contig_infos.swap(func.contig_infos());
  ClearStatus();
  return contigs.size();
}

void HashGraph::AssembleFunc::operator()(HashGraphVertex &vertex) {
  if (!vertex.status().Lock(0)) return;

  ContigBuilder contig_builder;
  contig_builder.Append(HashGraphVertexAdaptor(&vertex));

  if (!vertex.kmer().IsPalindrome()) {
    for (int strand = 0; strand < 2; ++strand) {
      HashGraphVertexAdaptor current(&vertex, strand);

      while (true) {
        HashGraphVertexAdaptor next;
        if (!hash_graph_->GetNextVertexAdaptor(current, next)) break;

        if (hash_graph_->IsPalindromeLoop(contig_builder.contig(), next)) break;

        if (hash_graph_->IsLoop(contig_builder.contig(), next)) return;

        if (!next.status().LockPreempt(0)) return;

        contig_builder.Append(next);
        current = next;
      }

      contig_builder.ReverseComplement();
    }
  }

  contigs_.push_back(contig_builder.contig());
  contig_infos_.push_back(contig_builder.contig_info());
}
