/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong & L3 Bioinformatics
 * Limited
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

#include "local_assembler.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include <omp.h>
#include "idba/contig_graph.h"
#include "idba/hash_graph.h"
#include "idba/sequence.h"
#include "kmlib/kmbit.h"

#include "sequence/io/contig/contig_reader.h"
#include "sequence/io/contig/contig_writer.h"
#include "utils/histgram.h"
#include "utils/utils.h"

void LocalAssembler::ReadContigs(const std::string &contig_file_name) {
  ContigReader reader(contig_file_name);
  reader.SetMinLen(min_contig_len_)->SetDiscardFlag(contig_flag::kLoop);
  contigs_.Clear();
  bool contig_reverse = false;
  reader.ReadAll(&contigs_, contig_reverse);
}

void LocalAssembler::BuildHashMapper() {
  size_t sz = contigs_.seq_count();
  size_t estimate_num_kmer = 0;

#pragma omp parallel for reduction(+ : estimate_num_kmer)
  for (size_t i = 0; i < sz; ++i) {
    estimate_num_kmer += (contigs_.GetSeqView(i).length() - seed_kmer_ + sparsity_) / sparsity_;
  }

  mapper_.reserve(estimate_num_kmer);

#pragma omp parallel for
  for (size_t i = 0; i < sz; ++i) {
    AddToHashMapper(mapper_, i, sparsity_);
  }
  xinfo("Number of contigs: {}, Mapper size: {}\n", contigs_.seq_count(), mapper_.size());
}

void LocalAssembler::AddReadLib(const std::string &file_prefix) {
  library_collection_.SetPath(file_prefix);
  library_collection_.Read(&reads_);
  insert_sizes_.resize(library_collection_.size(), tlen_t(-1, -1));
}

inline uint64_t EncodeContigOffset(unsigned contig_id, unsigned contig_offset, bool strand) {
  return (uint64_t(contig_id) << 32) | (contig_offset << 1) | strand;
}

inline void DecodeContigOffset(uint64_t code, uint32_t &contig_id, uint32_t &contig_offset, bool &strand) {
  contig_id = code >> 32;
  contig_offset = (code & 0xFFFFFFFFULL) >> 1;
  strand = code & 1ULL;
}

void LocalAssembler::AddToHashMapper(mapper_t &mapper, unsigned contig_id, int sparcity) {
  kmer_t key;
  auto contig_view = contigs_.GetSeqView(contig_id);
  for (int i = 0, len = contig_view.length(); i + seed_kmer_ <= len; i += sparcity) {
    auto ptr_and_offset = contig_view.raw_address(i);
    key.InitFromPtr(ptr_and_offset.first, ptr_and_offset.second, seed_kmer_);
    auto kmer = key.unique_format(seed_kmer_);
    auto offset = EncodeContigOffset(contig_id, i, key != kmer);
    std::lock_guard<SpinLock> lk(lock_);
    auto res = mapper.emplace(kmer, offset);
    if (!res.second) {
      res.first->second |= 1ULL << 63;
    }
  }
}

inline uint32_t GetWord(const uint32_t *first_word, uint32_t first_shift, int from, int len, bool strand) {
  int from_word_idx = (first_shift + from) / 16;
  int from_word_shift = (first_shift + from) % 16;
  uint32_t ret = *(first_word + from_word_idx) << from_word_shift * 2;
  assert(len <= 16);

  if (16 - from_word_shift < len) {
    ret |= *(first_word + from_word_idx + 1) >> (16 - from_word_shift) * 2;
  }

  if (len < 16) {
    ret >>= (16 - len) * 2;
    ret <<= (16 - len) * 2;
  }

  if (strand == 1) {
    ret = kmlib::bit::ReverseComplement<2>(ret);
    ret <<= (16 - len) * 2;
  }

  return ret;
}

inline int Mismatch(uint32_t x, uint32_t y) {
  x ^= y;
  x |= x >> 1;
  x &= 0x55555555U;
  return __builtin_popcount(x);
}

int LocalAssembler::Match(const SeqPackage::SeqView &seq_view, int query_from, int query_to, size_t contig_id,
                          int ref_from, int ref_to, bool strand) {
  auto query_ptr_and_offset = seq_view.raw_address();
  const uint32_t *query_first_word = query_ptr_and_offset.first;
  int query_shift = query_ptr_and_offset.second;

  auto contig_view = contigs_.GetSeqView(contig_id);
  auto ref_ptr_and_offset = contig_view.raw_address();
  const uint32_t *ref_first_word = ref_ptr_and_offset.first;
  int ref_shift = ref_ptr_and_offset.second;

  int match_len = query_to - query_from + 1;
  int threshold = similarity_ * match_len + 0.5;

  for (int i = query_from; i <= query_to; i += 16) {
    int len = std::min(16, query_to - i + 1);
    uint32_t qw = GetWord(query_first_word, query_shift, i, len, 0);
    int ref_i = strand == 0 ? ref_from + i - query_from : ref_to - (i + len - 1 - query_from);
    uint32_t rw = GetWord(ref_first_word, ref_shift, ref_i, len, strand);

    match_len -= Mismatch(qw, rw);

    if (match_len < threshold) {
      return 0;
    }
  }

  return match_len;
}

bool LocalAssembler::MapToHashMapper(const mapper_t &mapper, const SeqPackage::SeqView &seq_view, MappingRecord &rec) {
  int len = seq_view.length();

  if (len < seed_kmer_ || len < 50) return false;  // too short reads not reliable

  int tested = 0;
  MappingRecord tested_rec[3];

  auto ptr_and_offset = seq_view.raw_address();
  kmer_t kmer_f(ptr_and_offset.first, ptr_and_offset.second, seed_kmer_);
  kmer_t kmer_r = kmer_f;
  kmer_r.ReverseComplement(seed_kmer_);
  int num_mapped = 0;
  uint32_t contig_id, contig_offset;
  bool contig_strand;

  for (int i = seed_kmer_ - 1; i < len; ++i) {
    if (i >= seed_kmer_) {
      uint8_t ch = seq_view.base_at(i);
      kmer_f.ShiftAppend(ch, seed_kmer_);
      kmer_r.ShiftPreappend(3 - ch, seed_kmer_);
    }

    bool query_strand = kmer_f.cmp(kmer_r, seed_kmer_) <= 0 ? 0 : 1;

    auto iter = mapper.find(query_strand == 0 ? kmer_f : kmer_r);

    if (iter == mapper.end() || (iter->second >> 63) != 0) {
      continue;
    }

    DecodeContigOffset(iter->second, contig_id, contig_offset, contig_strand);

    auto contig_view = contigs_.GetSeqView(contig_id);

    assert(contig_id < contigs_.seq_count());
    assert(contig_offset < contig_view.length());

    bool mapping_strand = contig_strand ^ query_strand;
    int contig_from = mapping_strand == 0 ? contig_offset - (i - seed_kmer_ + 1) : contig_offset - (len - 1 - i);
    int contig_to = mapping_strand == 0 ? contig_offset + seed_kmer_ - 1 + len - 1 - i : contig_offset + i;
    contig_from = std::max(contig_from, 0);
    contig_to = std::min(static_cast<int>(contig_view.length() - 1), contig_to);

    if (contig_to - contig_from + 1 < len && contig_to - contig_from + 1 < min_mapped_len_) {
      continue;  // clipped alignment is considered iff its length >=
                 // min_mapped_len_
    }

    int query_from =
        mapping_strand == 0 ? i - (seed_kmer_ - 1) - (contig_offset - contig_from) : i - (contig_to - contig_offset);
    int query_to =
        mapping_strand == 0 ? i - (seed_kmer_ - 1) + (contig_to - contig_offset) : i + (contig_offset - contig_from);

    bool has_tested = false;

    for (int j = 0; j < tested; ++j) {
      if (contig_id == tested_rec[j].contig_id && contig_from == tested_rec[j].contig_from &&
          contig_to == tested_rec[j].contig_to && query_from == tested_rec[j].query_from &&
          query_to == tested_rec[j].query_to && mapping_strand == tested_rec[j].strand) {
        has_tested = true;
        break;
      }
    }

    if (has_tested) {
      continue;
    } else {
      if (tested >= 3) {
        tested--;
      }

      tested_rec[tested].contig_id = contig_id;
      tested_rec[tested].query_from = query_from;
      tested_rec[tested].query_to = query_to;
      tested_rec[tested].contig_from = contig_from;
      tested_rec[tested].contig_to = contig_to;
      tested_rec[tested].strand = mapping_strand;
      ++tested;
    }

    int match_bases = Match(seq_view, query_from, query_to, contig_id, contig_from, contig_to, mapping_strand);

    if (match_bases > 0) {
      if (num_mapped > 0) {
        return false;
      } else {
        rec = tested_rec[tested - 1];
        rec.mismatch = query_to - query_from + 1 - match_bases;
        num_mapped = 1;
      }
    }
  }

  return num_mapped == 1;
}

void LocalAssembler::EstimateInsertSize() {
  for (unsigned lib_id = 0; lib_id < library_collection_.size(); ++lib_id) {
    auto lib = library_collection_.GetLib(lib_id);
    if (!lib.IsPaired()) {
      continue;
    }

    MappingRecord rec1, rec2;
    Histgram<int> insert_hist;
    const size_t min_hist_size_for_estimation = 1u << 18;
    size_t processed_reads = 0;

    while (insert_hist.size() < min_hist_size_for_estimation && processed_reads < lib.seq_count()) {
      size_t start_read_id = processed_reads;
      processed_reads = std::min(min_hist_size_for_estimation + start_read_id, lib.seq_count());

#pragma omp parallel for private(rec1, rec2)
      for (size_t i = start_read_id; i < processed_reads; i += 2) {
        auto seq1 = lib.GetSequenceView(i);
        auto seq2 = lib.GetSequenceView(i + 1);
        if (MapToHashMapper(mapper_, seq1, rec1) && MapToHashMapper(mapper_, seq2, rec2)) {
          if (rec1.contig_id == rec2.contig_id && rec1.strand != rec2.strand) {
            int insert_size;

            if (rec1.strand == 0) {
              insert_size = rec2.contig_to + seq2.length() - rec2.query_to - (rec1.contig_from - rec1.query_from);
            } else {
              insert_size = rec1.contig_to + seq1.length() - rec1.query_to - (rec2.contig_from - rec2.query_from);
            }

            if (insert_size >= (int)seq1.length() && insert_size >= (int)seq2.length()) {
              insert_hist.insert(insert_size);
            }
          }
        }
      }
    }

    insert_hist.Trim(0.01);
    insert_sizes_[lib_id] = tlen_t(insert_hist.mean(), insert_hist.sd());

    xinfo("Lib {}, insert size: {.2} sd: {.2}\n", lib_id, insert_hist.mean(), insert_hist.sd());
  }
}

int LocalAssembler::LocalRange(int lib_id) {
  auto &lib = library_collection_.GetLib(lib_id);
  int local_range = lib.GetMaxLength() - 1;

  if (insert_sizes_[lib_id].first >= lib.GetMaxLength()) {
    local_range =
        std::min(2 * insert_sizes_[lib_id].first, insert_sizes_[lib_id].first + 3 * insert_sizes_[lib_id].second);
  }

  if (local_range > kMaxLocalRange) {
    local_range = kMaxLocalRange;
  }

  return local_range;
}

int LocalAssembler::AddToMappingDeque(size_t read_id, const MappingRecord &rec, int local_range) {
  assert(read_id < reads_.seq_count());
  assert(rec.contig_id < contigs_.seq_count());

  int contig_len = contigs_.GetSeqView(rec.contig_id).length();
  int read_len = reads_.GetSeqView(read_id).length();
  unsigned ret = 0;

  if (rec.contig_to < local_range && rec.query_from != 0 && rec.query_to == read_len - 1) {
    locks_.lock(rec.contig_id);
    mapped_f_[rec.contig_id].emplace_back(rec.contig_to, 0, rec.mismatch, rec.strand, read_id);
    locks_.unlock(rec.contig_id);
    ret++;
  } else if (rec.contig_from + local_range >= contig_len && rec.query_to < read_len - 1 && rec.query_from == 0) {
    locks_.lock(rec.contig_id);
    mapped_r_[rec.contig_id].emplace_back(contig_len - 1 - rec.contig_from, 0, rec.mismatch, rec.strand, read_id);
    locks_.unlock(rec.contig_id);
    ret++;
  }

  return ret;
}

int LocalAssembler::AddMateToMappingDeque(size_t read_id, size_t mate_id, const MappingRecord &rec1,
                                          const MappingRecord &rec2, bool mapped2, int local_range) {
  assert(read_id < reads_.seq_count());
  assert(mate_id < reads_.seq_count());
  assert(rec1.contig_id < contigs_.seq_count());
  assert(!mapped2 || rec2.contig_id < contigs_.seq_count());

  if (mapped2 && rec2.contig_id == rec1.contig_id) return 0;

  int contig_len = contigs_.GetSeqView(rec1.contig_id).length();
  unsigned ret = 0;

  if (rec1.contig_to < local_range && rec1.strand == 1) {
    locks_.lock(rec1.contig_id);
    mapped_f_[rec1.contig_id].emplace_back(rec1.contig_to, 1, rec1.mismatch, rec1.strand, mate_id);
    locks_.unlock(rec1.contig_id);
    ret++;
  } else if (rec1.contig_from + local_range >= contig_len && rec1.strand == 0) {
    locks_.lock(rec1.contig_id);
    mapped_r_[rec1.contig_id].emplace_back(contig_len - 1 - rec1.contig_from, 1, rec1.mismatch, rec1.strand, mate_id);
    locks_.unlock(rec1.contig_id);
    ret++;
  }

  return ret;
}

void LocalAssembler::MapToContigs() {
  mapped_f_.resize(contigs_.seq_count());
  mapped_r_.resize(contigs_.seq_count());
  locks_.reset(contigs_.seq_count());

  max_read_len_ = 1;
  local_range_ = 0;

  for (unsigned lib_id = 0; lib_id < library_collection_.size(); ++lib_id) {
    auto lib = library_collection_.GetLib(lib_id);
    int local_range = LocalRange(lib_id);
    bool is_paired = lib.IsPaired();

    local_range_ = std::max(local_range, local_range_);
    max_read_len_ = std::max(max_read_len_, static_cast<int>(lib.GetMaxLength()));

    MappingRecord rec1, rec2;
    size_t num_added = 0, num_mapped = 0;

    if (is_paired) {
#pragma omp parallel for private(rec1, rec2) reduction(+ : num_added, num_mapped)
      for (size_t i = 0; i < lib.seq_count(); i += 2) {
        auto seq1 = lib.GetSequenceView(i);
        auto seq2 = lib.GetSequenceView(i + 1);
        bool map1 = MapToHashMapper(mapper_, seq1, rec1);
        bool map2 = MapToHashMapper(mapper_, seq2, rec2);

        if (map1) {
          ++num_mapped;
          num_added += AddToMappingDeque(seq1.id(), rec1, local_range);
          num_added += AddMateToMappingDeque(seq1.id(), seq2.id(), rec1, rec2, map2, local_range);
        }

        if (map2) {
          ++num_mapped;
          num_added += AddToMappingDeque(seq2.id(), rec2, local_range);
          num_added += AddMateToMappingDeque(seq2.id(), seq1.id(), rec2, rec1, map1, local_range);
        }
      }
    } else {
#pragma omp parallel for private(rec1, rec2) reduction(+ : num_added, num_mapped)
      for (size_t i = 0; i < lib.seq_count(); ++i) {
        auto seq1 = lib.GetSequenceView(i);
        bool map1 = MapToHashMapper(mapper_, seq1, rec1);

        if (map1) {
          num_added += AddToMappingDeque(seq1.id(), rec1, local_range);
          ++num_mapped;
        }
      }
    }

    xinfo(
        "Lib {}: total {} reads, aligned {}, added {} reads for local "
        "assembly\n",
        lib_id, lib.seq_count(), num_mapped, num_added);
  }
  locks_.reset(0);
}

inline void LaunchIDBA(HashGraph &hash_graph, ContigGraph &contig_graph, std::deque<Sequence> &reads,
                       Sequence &contig_end, std::deque<Sequence> &out_contigs,
                       std::deque<ContigInfo> &out_contig_infos, int mink, int maxk, int step) {
  int local_range = contig_end.size();
  hash_graph.clear();
  contig_graph.clear();
  out_contigs.clear();
  out_contig_infos.clear();

  int max_read_len = 0;

  for (auto &read : reads) {
    max_read_len = std::max(max_read_len, (int)read.size());
  }

  for (int kmer_size = mink; kmer_size <= std::min(maxk, max_read_len); kmer_size += step) {
    int64_t sum = 0;
    hash_graph.set_kmer_size(kmer_size);
    hash_graph.clear();

    for (auto &read : reads) {
      if ((int)read.size() < kmer_size) continue;

      const Sequence &seq(read);
      hash_graph.InsertKmers(seq);
      sum += seq.size() - kmer_size + 1;
    }

    auto histgram = hash_graph.coverage_histgram();
    double mean = histgram.percentile(1 - 1.0 * local_range / hash_graph.num_vertices());
    double threshold = mean;

    hash_graph.InsertKmers(contig_end);

    for (const auto &out_contig : out_contigs) hash_graph.InsertUncountKmers(out_contig);

    hash_graph.Assemble(out_contigs, out_contig_infos);

    contig_graph.set_kmer_size(kmer_size);
    contig_graph.Initialize(out_contigs, out_contig_infos);
    contig_graph.RemoveDeadEnd(kmer_size * 2);

    contig_graph.RemoveBubble();
    contig_graph.IterateCoverage(kmer_size * 2, 1, threshold);

    contig_graph.Assemble(out_contigs, out_contig_infos);

    if (out_contigs.size() == 1) {
      break;
    }
  }
}

void LocalAssembler::LocalAssemble() {
  int min_num_reads = local_range_ / max_read_len_;

  Sequence seq, contig_end;
  HashGraph hash_graph;
  ContigGraph contig_graph;
  std::deque<Sequence> reads;
  std::deque<Sequence> out_contigs;
  std::deque<ContigInfo> out_contig_infos;

  ContigWriter local_contig_writer(local_filename_);

#pragma omp parallel for private(hash_graph, contig_graph, seq, contig_end, reads, out_contigs, out_contig_infos) \
    schedule(dynamic)
  for (uint64_t cid = 0; cid < contigs_.seq_count(); ++cid) {
    auto contig_view = contigs_.GetSeqView(cid);
    int cl = contig_view.length();

    for (int strand = 0; strand < 2; ++strand) {
      auto &mapped_reads = strand == 0 ? mapped_f_[cid] : mapped_r_[cid];
      if (static_cast<int>(mapped_reads.size()) <= min_num_reads) {
        continue;
      }

      // collect local reads, convert them into Sequence
      reads.clear();

      std::sort(mapped_reads.begin(), mapped_reads.end());
      uint64_t last_mapping_pos = -1;
      int pos_count = 0;

      for (auto &mapped_read : mapped_reads) {
        uint64_t pos = (mapped_read.contig_offset << 1) | mapped_read.is_mate;
        pos_count = pos == last_mapping_pos ? pos_count + 1 : 1;
        last_mapping_pos = pos;

        if (pos_count <= 3) {
          seq.clear();
          auto read_view = reads_.GetSeqView(mapped_read.read_id);

          for (unsigned ri = 0, rsz = read_view.length(); ri < rsz; ++ri) {
            seq.Append(read_view.base_at(ri));
          }
          reads.push_back(seq);
        }
      }

      contig_end.clear();

      if (strand == 0) {
        for (int j = 0, e = std::min(local_range_, cl); j < e; ++j) {
          contig_end.Append(contig_view.base_at(j));
        }
      } else {
        for (int j = std::max(0, cl - local_range_); j < cl; ++j) {
          contig_end.Append(contig_view.base_at(j));
        }
      }

      out_contigs.clear();
      LaunchIDBA(hash_graph, contig_graph, reads, contig_end, out_contigs, out_contig_infos, local_kmin_, local_kmax_,
                 local_step_);

      for (uint64_t j = 0; j < out_contigs.size(); ++j) {
        if (out_contigs[j].size() > min_contig_len_ && out_contigs[j].size() > local_kmax_) {
          auto str = out_contigs[j].str();
          local_contig_writer.WriteLocalContig(str, cid, strand, j);
        }
      }
    }
  }
}