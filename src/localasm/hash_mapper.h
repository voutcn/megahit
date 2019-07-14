//
// Created by vout on 7/10/19.
//

#ifndef MEGAHIT_LOCALASM_HASH_MAPPER_H
#define MEGAHIT_LOCALASM_HASH_MAPPER_H

#include <cstdint>

#include "parallel_hashmap/phmap.h"
#include "sequence/kmer_plus.h"
#include "sequence/sequence_package.h"

struct MappingRecord {
  uint32_t contig_id;
  int32_t contig_from;
  int32_t contig_to;
  uint64_t query_id;
  int32_t query_from;
  int32_t query_to;
  uint32_t mismatch;
  uint8_t strand;
  bool valid;

  bool operator<(const MappingRecord &rhs) const {
    if (contig_id != rhs.contig_id) return contig_id < rhs.contig_id;
    if (contig_from != rhs.contig_from) return contig_from < rhs.contig_from;
    if (contig_to != rhs.contig_to) return contig_to < rhs.contig_to;
    if (query_id != rhs.query_id) return query_id < rhs.query_id;
    if (query_from != rhs.query_from) return query_from < rhs.query_from;
    if (query_to != rhs.query_to) return query_to < rhs.query_to;
    return strand < rhs.strand;
  }

  bool operator==(const MappingRecord &rhs) const {
    return contig_id == rhs.contig_id && contig_from == rhs.contig_from &&
           contig_to == rhs.contig_to && query_id == rhs.query_id &&
           query_from == rhs.query_from && query_to == rhs.query_to &&
           strand == rhs.strand;
  };
};

class HashMapper {
 public:
  using TKmer = Kmer<2, uint32_t>;
  using TMapper = phmap::flat_hash_map<TKmer, uint64_t, KmerHash>;

  void LoadAndBuild(const std::string &contig_file, int32_t min_len,
                    int32_t seed_kmer_size, int32_t sparsity);

  void SetMappingThreshold(int32_t mapping_len, double similarity) {
    min_mapped_len_ = mapping_len;
    similarity_ = similarity;
  }

  MappingRecord TryMap(const SeqPackage ::SeqView &seq_view) const;

  const SeqPackage &refseq() const { return refseq_; }

 private:
  int32_t Match(const SeqPackage::SeqView &seq_view, int query_from,
                int query_to, size_t contig_id, int ref_from, int ref_to,
                bool strand) const;

 private:
  TMapper index_;
  SeqPackage refseq_;

  int32_t seed_kmer_size_{31};
  int32_t min_mapped_len_{50};
  double similarity_{0.95};
};

#endif  // MEGAHIT_LOCALASM_HASH_MAPPER_H
