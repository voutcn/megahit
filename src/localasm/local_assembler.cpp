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

#include "sequence/lib_io.h"
#include "sequence/readers/contig_reader.h"
#include "utils/histgram.h"
#include "utils/safe_open.h"
#include "utils/utils.h"

void LocalAssembler::ReadContigs(const std::string &contig_file_name) {
  ContigReader reader(contig_file_name);
  reader.SetMinLen(min_contig_len_)->SetDiscardFlag(contig_flag::kLoop);
  contigs_.Clear();
  bool contig_reverse = false;
  reader.ReadAll(&contigs_, contig_reverse);
}

void LocalAssembler::BuildHashMapper(bool show_stat) {
  size_t sz = contigs_.SeqCount();
  size_t estimate_num_kmer = 0;

#pragma omp parallel for reduction(+ : estimate_num_kmer)
  for (size_t i = 0; i < sz; ++i) {
    estimate_num_kmer += (contigs_.SequenceLength(i) - seed_kmer_ + sparsity_) / sparsity_;
  }

  mapper_.reserve(estimate_num_kmer);

#pragma omp parallel for
  for (size_t i = 0; i < sz; ++i) {
    AddToHashMapper(mapper_, i, sparsity_);
  }

  if (show_stat) {
    xinfo("Number of contigs: %lu, Mapper size: %lu\n", contigs_.SeqCount(), mapper_.size());
  }
}

void LocalAssembler::AddReadLib(const std::string &file_prefix) {
  ReadBinaryLibs(file_prefix, reads_, lib_info_);
  insert_sizes_.resize(lib_info_.size(), tlen_t(-1, -1));
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

  for (int i = 0, len = contigs_.SequenceLength(contig_id); i + seed_kmer_ <= len; i += sparcity) {
    auto ptr_and_offset = contigs_.WordPtrAndOffset(contig_id, i);
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

int LocalAssembler::Match(size_t read_id, int query_from, int query_to, size_t contig_id, int ref_from, int ref_to,
                          bool strand) {
  auto query_ptr_and_offset = reads_.WordPtrAndOffset(read_id);
  const uint32_t *query_first_word = query_ptr_and_offset.first;
  int query_shift = query_ptr_and_offset.second;
  auto ref_ptr_and_offset = contigs_.WordPtrAndOffset(contig_id);
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

bool LocalAssembler::MapToHashMapper(const mapper_t &mapper, size_t read_id, MappingRecord &rec) {
  int len = reads_.SequenceLength(read_id);

  if (len < seed_kmer_ || len < 50) return false;  // too short reads not reliable

  int tested = 0;
  MappingRecord tested_rec[3];

  auto ptr_and_offset = reads_.WordPtrAndOffset(read_id);
  kmer_t kmer_f(ptr_and_offset.first, ptr_and_offset.second, seed_kmer_);
  kmer_t kmer_r = kmer_f;
  kmer_r.ReverseComplement(seed_kmer_);
  int num_mapped = 0;
  uint32_t contig_id, contig_offset;
  bool contig_strand;

  for (int i = seed_kmer_ - 1; i < len; ++i) {
    if (i >= seed_kmer_) {
      uint8_t ch = reads_.GetBase(read_id, i);
      kmer_f.ShiftAppend(ch, seed_kmer_);
      kmer_r.ShiftPreappend(3 - ch, seed_kmer_);
    }

    bool query_strand = kmer_f.cmp(kmer_r, seed_kmer_) <= 0 ? 0 : 1;

    auto iter = mapper.find(query_strand == 0 ? kmer_f : kmer_r);

    if (iter == mapper.end() || (iter->second >> 63) != 0) {
      continue;
    }

    DecodeContigOffset(iter->second, contig_id, contig_offset, contig_strand);
    assert(contig_id < contigs_.SeqCount());
    assert(contig_offset < contigs_.SequenceLength(contig_id));

    bool mapping_strand = contig_strand ^ query_strand;
    int contig_from = mapping_strand == 0 ? contig_offset - (i - seed_kmer_ + 1) : contig_offset - (len - 1 - i);
    int contig_to = mapping_strand == 0 ? contig_offset + seed_kmer_ - 1 + len - 1 - i : contig_offset + i;
    contig_from = std::max(contig_from, 0);
    contig_to = std::min((int)contigs_.SequenceLength(contig_id) - 1, contig_to);

    if (contig_to - contig_from + 1 < len && contig_to - contig_from + 1 < min_mapped_len_) {
      continue;  // clipped alignment is considered iff its length >= min_mapped_len_
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

    int match_bases = Match(read_id, query_from, query_to, contig_id, contig_from, contig_to, mapping_strand);

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

void LocalAssembler::EstimateInsertSize(bool show_stat) {
  for (unsigned lib_id = 0; lib_id < lib_info_.size(); ++lib_id) {
    if (!lib_info_[lib_id].is_pe) {
      continue;
    }

    MappingRecord rec1, rec2;
    Histgram<int> insert_hist;
    int64_t start_read_id = lib_info_[lib_id].from;
    int64_t end_read_id = start_read_id;

    while (insert_hist.size() < (1 << 18) && end_read_id <= lib_info_[lib_id].to) {
      start_read_id = end_read_id;
      end_read_id = std::min(lib_info_[lib_id].to + 1, start_read_id + (2 << 18));

#pragma omp parallel for private(rec1, rec2)
      for (int64_t i = start_read_id; i < end_read_id; i += 2) {
        if (MapToHashMapper(mapper_, i, rec1) && MapToHashMapper(mapper_, i + 1, rec2)) {
          if (rec1.contig_id == rec2.contig_id && rec1.strand != rec2.strand) {
            int insert_size;

            if (rec1.strand == 0) {
              insert_size =
                  rec2.contig_to + reads_.SequenceLength(i + 1) - rec2.query_to - (rec1.contig_from - rec1.query_from);
            } else {
              insert_size =
                  rec1.contig_to + reads_.SequenceLength(i) - rec1.query_to - (rec2.contig_from - rec2.query_from);
            }

            if (insert_size >= (int)reads_.SequenceLength(i) && insert_size >= (int)reads_.SequenceLength(i + 1)) {
              insert_hist.insert(insert_size);
            }
          }
        }
      }
    }

    insert_hist.Trim(0.01);
    insert_sizes_[lib_id] = tlen_t(insert_hist.mean(), insert_hist.sd());

    if (show_stat) {
      xinfo("Lib %d, insert size: %.2lf sd: %.2lf\n", lib_id, insert_hist.mean(), insert_hist.sd());
    }
  }
}

int LocalAssembler::LocalRange(int lib_id) {
  int local_range = lib_info_[lib_id].max_read_len - 1;

  if (insert_sizes_[lib_id].first >= lib_info_[lib_id].max_read_len) {
    local_range =
        std::min(2 * insert_sizes_[lib_id].first, insert_sizes_[lib_id].first + 3 * insert_sizes_[lib_id].second);
  }

  if (local_range > kMaxLocalRange) {
    local_range = kMaxLocalRange;
  }

  return local_range;
}

inline uint64_t PackMappingResult(uint64_t contig_offset, uint64_t is_mate, uint64_t mismatch, uint64_t strand,
                                  uint64_t read_id) {
  assert(contig_offset < (1ULL << 14));
  return (contig_offset << 50) | (is_mate << 49) | (std::min(uint64_t(15), mismatch) << 45) | (strand << 44) | read_id;
}

int LocalAssembler::AddToMappingDeque(size_t read_id, const MappingRecord &rec, int local_range) {
  assert(read_id < reads_.SeqCount());
  assert(rec.contig_id < contigs_.SeqCount());

  int contig_len = contigs_.SequenceLength(rec.contig_id);
  int read_len = reads_.SequenceLength(read_id);
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
  assert(read_id < reads_.SeqCount());
  assert(mate_id < reads_.SeqCount());
  assert(rec1.contig_id < contigs_.SeqCount());
  assert(!mapped2 || rec2.contig_id < contigs_.SeqCount());

  if (mapped2 && rec2.contig_id == rec1.contig_id) return 0;

  int contig_len = contigs_.SequenceLength(rec1.contig_id);
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
  mapped_f_.resize(contigs_.SeqCount());
  mapped_r_.resize(contigs_.SeqCount());
  locks_.reset(contigs_.SeqCount());

  max_read_len_ = 1;
  local_range_ = 0;

  for (unsigned lib_id = 0; lib_id < lib_info_.size(); ++lib_id) {
    int local_range = LocalRange(lib_id);
    bool is_paired = lib_info_[lib_id].is_pe;

    local_range_ = std::max(local_range, local_range_);
    max_read_len_ = std::max(max_read_len_, lib_info_[lib_id].max_read_len);

    MappingRecord rec1, rec2;
    size_t num_added = 0, num_mapped = 0;

#pragma omp parallel for private(rec1, rec2) reduction(+ : num_added, num_mapped) schedule(static, 1)
    for (int64_t i = lib_info_[lib_id].from; i <= lib_info_[lib_id].to; i += 2) {
      bool map1 = MapToHashMapper(mapper_, i, rec1);
      bool map2 = (i + 1) <= lib_info_[lib_id].to ? MapToHashMapper(mapper_, i + 1, rec2) : false;

      if (map1) {
        num_added += AddToMappingDeque(i, rec1, local_range);
        ++num_mapped;

        if (is_paired) {
          num_added += AddMateToMappingDeque(i, i + 1, rec1, rec2, map2, local_range);
        }
      }

      if (map2) {
        ++num_mapped;
        num_added += AddToMappingDeque(i + 1, rec2, local_range);

        if (is_paired) {
          num_added += AddMateToMappingDeque(i + 1, i, rec2, rec1, map1, local_range);
        }
      }
    }

    xinfo("Lib %d: total %ld reads, aligned %lu, added %lu reads for local assembly\n", lib_id,
          lib_info_[lib_id].to - lib_info_[lib_id].from + 1, num_mapped, num_added);
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

  auto local_file = xfopen(local_filename_.c_str(), "w");
  auto local_info = xfopen((local_filename_ + ".info").c_str(), "w");

  long long num_bases = 0;
  long long num_contigs = 0;

#pragma omp parallel for private(hash_graph, contig_graph, seq, contig_end, reads, out_contigs, out_contig_infos) schedule(dynamic) reduction(+ : num_contigs, num_bases)
  for (uint64_t cid = 0; cid < contigs_.SeqCount(); ++cid) {
    int cl = contigs_.SequenceLength(cid);

    for (int strand = 0; strand < 2; ++strand) {
      auto &mapped_reads = strand == 0 ? mapped_f_[cid] : mapped_r_[cid];
      if ((int)mapped_reads.size() <= min_num_reads) {
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
          uint64_t read_id = mapped_read.read_id;

          for (unsigned ri = 0, rsz = reads_.SequenceLength(read_id); ri < rsz; ++ri) {
            seq.Append(reads_.GetBase(read_id, ri));
          }
          reads.push_back(seq);
        }
      }

      contig_end.clear();

      if (strand == 0) {
        for (int j = 0, e = std::min(local_range_, cl); j < e; ++j) {
          contig_end.Append(contigs_.GetBase(cid, j));
        }
      } else {
        for (int j = std::max(0, cl - local_range_); j < cl; ++j) {
          contig_end.Append(contigs_.GetBase(cid, j));
        }
      }

      out_contigs.clear();
      LaunchIDBA(hash_graph, contig_graph, reads, contig_end, out_contigs, out_contig_infos, local_kmin_, local_kmax_,
                 local_step_);

      for (uint64_t j = 0; j < out_contigs.size(); ++j) {
        if (out_contigs[j].size() > min_contig_len_ && out_contigs[j].size() > local_kmax_) {
          auto str = out_contigs[j].str();
          fprintf(local_file, ">lc_%llu_strand_%d_id_%llu flag=0 multi=1\n%s\n", static_cast<unsigned long long>(cid),
                  strand, static_cast<unsigned long long>(j), str.c_str());
          num_contigs++;
          num_bases += out_contigs[j].size();
        }
      }
    }
  }

  fprintf(local_info, "%lld %lld\n", num_contigs, num_bases);

  fclose(local_file);
  fclose(local_info);
}