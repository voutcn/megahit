//
// Created by vout on 7/13/19.
//

#ifndef MEGAHIT_LOCALASM_MAPPING_RESULT_H
#define MEGAHIT_LOCALASM_MAPPING_RESULT_H

#include <cstdint>
#include <deque>
#include <vector>
#include "hash_mapper.h"
#include "kmlib/kmbitvector.h"

class MappingResultCollector {
 public:
  explicit MappingResultCollector(size_t n_ref)
      : fwd_mappings_(n_ref), bwd_mappings_(n_ref), locks_(n_ref) {}

  unsigned AddSingle(const MappingRecord &rec, int32_t contig_len,
                     int32_t read_len, int32_t local_range) {
    unsigned ret = 0;

    if (rec.contig_to < local_range && rec.query_from != 0 &&
        rec.query_to == read_len - 1) {
      locks_.lock(rec.contig_id);
      fwd_mappings_[rec.contig_id].push_back(EncodeMappingRead(
          rec.contig_to, 0, rec.mismatch, rec.strand, rec.query_id));
      locks_.unlock(rec.contig_id);
      ret++;
    } else if (rec.contig_from + local_range >= contig_len &&
               rec.query_to < read_len - 1 && rec.query_from == 0) {
      locks_.lock(rec.contig_id);
      bwd_mappings_[rec.contig_id].push_back(
          EncodeMappingRead(contig_len - 1 - rec.contig_from, 0, rec.mismatch,
                            rec.strand, rec.query_id));
      locks_.unlock(rec.contig_id);
      ret++;
    }

    return ret;
  }

  unsigned AddMate(const MappingRecord &rec1, const MappingRecord &rec2,
                   int32_t contig_len, uint64_t mate_id, int32_t local_range) {
    if (rec2.valid && rec2.contig_id == rec1.contig_id) return 0;
    unsigned ret = 0;

    if (rec1.contig_to < local_range && rec1.strand == 1) {
      locks_.lock(rec1.contig_id);
      fwd_mappings_[rec1.contig_id].push_back(EncodeMappingRead(
          rec1.contig_to, 1, rec1.mismatch, rec1.strand, mate_id));
      locks_.unlock(rec1.contig_id);
      ret++;
    } else if (rec1.contig_from + local_range >= contig_len &&
               rec1.strand == 0) {
      locks_.lock(rec1.contig_id);
      bwd_mappings_[rec1.contig_id].push_back(
          EncodeMappingRead(contig_len - 1 - rec1.contig_from, 1, rec1.mismatch,
                            rec1.strand, mate_id));
      locks_.unlock(rec1.contig_id);
      ret++;
    }

    return ret;
  }

  const std::deque<uint64_t> &GetMappingResults(int64_t contig_id,
                                                uint8_t strand) {
    auto &results = strand == 0 ? fwd_mappings_ : bwd_mappings_;
    auto &ret = results[contig_id];
    std::sort(ret.begin(), ret.end());
    return ret;
  }

  static uint64_t GetContigAbsPos(uint64_t encoded) {
    return encoded >> (44u + 1u + 4u);
  }

  static uint64_t GetReadId(uint64_t encoded) {
    return encoded & ((1ull << 44u) - 1);
  }

 private:
  static uint64_t EncodeMappingRead(uint32_t contig_offset, uint8_t is_mate,
                                    uint32_t mismatch, uint8_t strand,
                                    uint64_t read_id) {
    assert(contig_offset <= (1u << 14u));
    assert(strand <= 1);
    uint64_t ret = contig_offset;
    ret = (ret << 1u) | is_mate;
    ret = (ret << 4u) | (mismatch < 15 ? mismatch : 15u);
    ret = (ret << 1u) | strand;
    ret = (ret << 44u) | read_id;  // 44 bits for read id
    return ret;
  }

  std::vector<std::deque<uint64_t>> fwd_mappings_;
  std::vector<std::deque<uint64_t>> bwd_mappings_;
  AtomicBitVector locks_;
};

#endif  // MEGAHIT_LOCALASM_MAPPING_RESULT_H
