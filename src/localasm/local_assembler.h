/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong & L3 Bioinformatics Limited
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* contact: Dinghua Li <dhli@cs.hku.hk> */

#ifndef LOCAL_ASSEMBLER_H
#define LOCAL_ASSEMBLER_H

#include <stdint.h>
#include <deque>
#include <string>
#include <vector>

#include "kmlib/kmbitvector.h"
#include "parallel_hashmap/phmap.h"
#include "sequence/kmer_plus.h"
#include "sequence/io/sequence_lib.h"
#include "sequence/sequence_package.h"
#include "utils/mutex.h"

struct LocalAssembler {
  static const int kMaxLocalRange = 650;

  struct MappingRecord {
    uint32_t contig_id;
    int32_t contig_from : 28;
    int32_t contig_to : 28;
    int16_t query_from : 15;
    int16_t query_to : 15;
    bool strand : 1;
    int mismatch : 9;
  };

  struct MappedReadRecord {
    uint64_t read_id : 44;
    uint8_t strand : 1;
    uint8_t mismatch : 4;
    bool is_mate : 1;
    uint16_t contig_offset : 14;
    MappedReadRecord() = default;
    MappedReadRecord(uint16_t contig_offset, bool is_mate, uint8_t mismatch, uint8_t strand, uint64_t read_id)
        : read_id(read_id),
          strand(strand),
          mismatch(mismatch < 15 ? mismatch : 15),
          is_mate(is_mate),
          contig_offset(contig_offset) {
      assert(contig_offset < (1 << 14));
      assert(strand < 2);
    }
    bool operator<(const MappedReadRecord &rhs) const {
      return *reinterpret_cast<const uint64_t *>(this) < *reinterpret_cast<const uint64_t *>(&rhs);
    }
  };

  typedef std::pair<double, double> tlen_t;
  typedef Kmer<2, uint32_t> kmer_t;
  typedef phmap::flat_hash_map<kmer_t, uint64_t, KmerHash> mapper_t;

  unsigned min_contig_len_;  // only align reads to these contigs
  int seed_kmer_;            // kmer size for seeding
  double similarity_;        // similarity threshold for alignment
  int sparsity_;             // sparsity of hash mapper
  int min_mapped_len_;
  std::string local_filename_;
  unsigned local_kmin_, local_kmax_, local_step_;

  // auto calculated after calling EstimateInsertSize() and MapToContigs()
  int local_range_{};
  int max_read_len_{};

  SeqPackage contigs_;
  SeqPackage reads_;
  mapper_t mapper_;
  SequenceLibCollection library_collection_;
  std::vector<tlen_t> insert_sizes_;
  SpinLock lock_;
  AtomicBitVector locks_;

  std::vector<std::deque<MappedReadRecord> > mapped_f_, mapped_r_;

  LocalAssembler(int min_contig_len, int seed_kmer, int sparsity)
      : min_contig_len_(min_contig_len), seed_kmer_(seed_kmer), sparsity_(sparsity) {
    similarity_ = 0.95;
    min_mapped_len_ = 100;
    local_kmin_ = 21;
    local_kmax_ = 41;
    local_step_ = 6;
  }

  ~LocalAssembler() = default;

  void SetKmerSize(int kmin, int kmax, int step) {
    local_kmin_ = kmin;
    local_kmax_ = kmax;
    local_step_ = step;
  }

  void SetMappingThreshold(double similarity, int mapping_len) {
    similarity_ = similarity;
    min_mapped_len_ = mapping_len;
  }

  void SetLocalContigFile(const std::string &local_filename) { local_filename_ = local_filename; }

  void ReadContigs(const std::string &file_name);
  void BuildHashMapper();
  void AddReadLib(const std::string &file_prefix);
  void EstimateInsertSize();
  void MapToContigs();
  void LocalAssemble();

  void AddToHashMapper(mapper_t &mapper, unsigned contig_id, int sparsity);
  int Match(const SeqPackage::SeqView &seq_view,
            int query_from,
            int query_to,
            size_t contig_id,
            int ref_from,
            int ref_to,
            bool strand);
  int LocalRange(int lib_id);
  int AddToMappingDeque(size_t read_id, const MappingRecord &rec, int local_range);
  int AddMateToMappingDeque(size_t read_id, size_t mate_id, const MappingRecord &rec1, const MappingRecord &rec2,
                            bool mapped2, int local_range);
  bool MapToHashMapper(const mapper_t &mapper, const SeqPackage::SeqView &seq_view, MappingRecord &rec);
};

#endif