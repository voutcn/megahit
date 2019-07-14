//
// Created by vout on 7/10/19.
//

#include "hash_mapper.h"
#include "sequence/io/contig/contig_reader.h"
#include "utils/utils.h"

namespace {

inline uint64_t EncodeContigOffset(uint32_t contig_id, uint32_t contig_offset,
                                   uint8_t strand) {
  return (uint64_t(contig_id) << 32) | (contig_offset << 1) | strand;
}

inline void DecodeContigOffset(uint64_t code, uint32_t &contig_id,
                               uint32_t &contig_offset, uint8_t &strand) {
  contig_id = code >> 32;
  contig_offset = (code & 0xFFFFFFFFULL) >> 1;
  strand = code & 1ULL;
}

inline uint32_t GetWord(const uint32_t *first_word, uint32_t first_shift,
                        int from, int len, bool strand) {
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
  return kmlib::bit::Popcount(x);
}

}  // namespace

void HashMapper::LoadAndBuild(const std::string &contig_file, int32_t min_len,
                              int32_t seed_kmer_size, int32_t sparsity) {
  seed_kmer_size_ = seed_kmer_size;
  ContigReader reader(contig_file);
  reader.SetMinLen(min_len)->SetDiscardFlag(contig_flag::kLoop);
  auto sizes = reader.GetNumContigsAndBases();
  refseq_.Clear();
  refseq_.ReserveSequences(sizes.first);
  refseq_.ReserveBases(sizes.second);
  bool contig_reverse = false;
  reader.ReadAll(&refseq_, contig_reverse);

  size_t sz = refseq_.seq_count();
  size_t n_kmers = 0;

#pragma omp parallel for reduction(+ : n_kmers)
  for (size_t i = 0; i < sz; ++i) {
    n_kmers +=
        (refseq_.GetSeqView(i).length() - seed_kmer_size + sparsity) / sparsity;
  }

  index_.reserve(n_kmers);
  SpinLock mapper_lock;

#pragma omp parallel for
  for (size_t i = 0; i < sz; ++i) {
    TKmer key;
    auto contig_view = refseq_.GetSeqView(i);
    for (int j = 0, len = contig_view.length(); j + seed_kmer_size <= len;
         j += sparsity) {
      auto ptr_and_offset = contig_view.raw_address(j);
      key.InitFromPtr(ptr_and_offset.first, ptr_and_offset.second,
                      seed_kmer_size);
      auto kmer = key.unique_format(seed_kmer_size);
      auto offset = EncodeContigOffset(contig_view.id(), j, key != kmer);
      std::lock_guard<SpinLock> lk(mapper_lock);
      auto res = index_.emplace(kmer, offset);
      if (!res.second) {
        res.first->second |= 1ULL << 63;
      }
    }
  }

  xinfo("Number of contigs: {}, index size: {}\n", refseq_.seq_count(),
        index_.size());
}

int32_t HashMapper::Match(const SeqPackage::SeqView &seq_view, int query_from,
                          int query_to, size_t contig_id, int ref_from,
                          int ref_to, bool strand) const {
  auto query_ptr_and_offset = seq_view.raw_address();
  const uint32_t *query_first_word = query_ptr_and_offset.first;
  int query_shift = query_ptr_and_offset.second;

  auto contig_view = refseq_.GetSeqView(contig_id);
  auto ref_ptr_and_offset = contig_view.raw_address();
  const uint32_t *ref_first_word = ref_ptr_and_offset.first;
  int ref_shift = ref_ptr_and_offset.second;

  int match_len = query_to - query_from + 1;
  int threshold = lround(similarity_ * match_len);

  for (int i = query_from; i <= query_to; i += 16) {
    int len = std::min(16, query_to - i + 1);
    uint32_t qw = GetWord(query_first_word, query_shift, i, len, 0);
    int ref_i = strand == 0 ? ref_from + i - query_from
                            : ref_to - (i + len - 1 - query_from);
    uint32_t rw = GetWord(ref_first_word, ref_shift, ref_i, len, strand);

    match_len -= Mismatch(qw, rw);

    if (match_len < threshold) {
      return 0;
    }
  }

  return match_len;
}

MappingRecord HashMapper::TryMap(const SeqPackage::SeqView &seq_view) const {
  MappingRecord bad_record;
  bad_record.valid = false;

  int len = seq_view.length();
  if (len < seed_kmer_size_ || len < 50)
    return bad_record;  // too short reads not reliable

  // small vector optimization
  static const int kArraySize = 3;
  std::array<MappingRecord, kArraySize> mapping_records;
  int n_mapping_records = 0;
  std::unique_ptr<std::vector<MappingRecord>> v_mapping_records;

  auto ptr_and_offset = seq_view.raw_address();
  TKmer kmer_f(ptr_and_offset.first, ptr_and_offset.second, seed_kmer_size_);
  TKmer kmer_r = kmer_f;
  kmer_r.ReverseComplement(seed_kmer_size_);

  for (int i = seed_kmer_size_ - 1; i < len; ++i) {
    if (i >= seed_kmer_size_) {
      uint8_t ch = seq_view.base_at(i);
      kmer_f.ShiftAppend(ch, seed_kmer_size_);
      kmer_r.ShiftPreappend(3 - ch, seed_kmer_size_);
    }

    uint8_t query_strand = kmer_f.cmp(kmer_r, seed_kmer_size_) <= 0 ? 0 : 1;

    auto iter = index_.find(query_strand == 0 ? kmer_f : kmer_r);

    if (iter == index_.end() || (iter->second >> 63) != 0) {
      continue;
    }

    uint32_t contig_id, contig_offset;
    uint8_t contig_strand;
    DecodeContigOffset(iter->second, contig_id, contig_offset, contig_strand);

    auto contig_view = refseq_.GetSeqView(contig_id);
    assert(contig_offset < contig_view.length());

    uint8_t mapping_strand = contig_strand ^ query_strand;
    int32_t contig_from = mapping_strand == 0
                              ? contig_offset - (i - seed_kmer_size_ + 1)
                              : contig_offset - (len - 1 - i);
    int32_t contig_to = mapping_strand == 0
                            ? contig_offset + seed_kmer_size_ - 1 + len - 1 - i
                            : contig_offset + i;
    contig_from = std::max(contig_from, 0);
    contig_to =
        std::min(static_cast<int32_t>(contig_view.length() - 1), contig_to);

    if (contig_to - contig_from + 1 < len &&
        contig_to - contig_from + 1 < min_mapped_len_) {
      continue;  // clipped alignment is considered iff its length >=
      // min_mapped_len_
    }

    int32_t query_from = mapping_strand == 0 ? i - (seed_kmer_size_ - 1) -
                                                   (contig_offset - contig_from)
                                             : i - (contig_to - contig_offset);
    int32_t query_to = mapping_strand == 0 ? i - (seed_kmer_size_ - 1) +
                                                 (contig_to - contig_offset)
                                           : i + (contig_offset - contig_from);

    assert(query_from >= 0 &&
           static_cast<uint32_t>(query_from) < seq_view.length());
    assert(query_to >= 0 &&
           static_cast<uint32_t>(query_to) < seq_view.length());

    auto rec = MappingRecord{contig_id,  contig_from,
                             contig_to,  static_cast<uint64_t>(seq_view.id()),
                             query_from, query_to,
                             0,          mapping_strand,
                             true};
    auto end = mapping_records.begin() + n_mapping_records;
    if (std::find(mapping_records.begin(), end, rec) == end) {
      if (n_mapping_records < kArraySize) {
        mapping_records[n_mapping_records] = rec;
        n_mapping_records++;
      } else {
        if (v_mapping_records.get() == nullptr) {
          v_mapping_records.reset(new std::vector<MappingRecord>(1, rec));
        } else {
          v_mapping_records->push_back(rec);
        }
      }
    }
  }

  if (n_mapping_records == 0) {
    return bad_record;
  }

  MappingRecord *best = &bad_record;
  int32_t max_match = 0;

#define CHECK_BEST_UNIQ(rec)                                              \
  do {                                                                    \
    int32_t match_bases =                                                 \
        Match(seq_view, rec.query_from, rec.query_to, rec.contig_id,      \
              rec.contig_from, rec.contig_to, rec.strand);                \
    if (match_bases == max_match) {                                       \
      best = &bad_record;                                                 \
    } else if (match_bases > max_match) {                                 \
      max_match = match_bases;                                            \
      int32_t mismatch = rec.query_to - rec.query_from + 1 - match_bases; \
      rec.mismatch = mismatch;                                            \
      best = &rec;                                                        \
    }                                                                     \
  } while (0)

  if (v_mapping_records.get() != nullptr) {
    if (v_mapping_records->size() > 1) {
      std::sort(v_mapping_records->begin(), v_mapping_records->end());
      v_mapping_records->resize(
          std::unique(v_mapping_records->begin(), v_mapping_records->end()) -
          v_mapping_records->begin());
    }

    for (auto &rec : *v_mapping_records) {
      CHECK_BEST_UNIQ(rec);
    }
  }

  for (int i = 0; i < n_mapping_records; ++i) {
    auto &rec = mapping_records[i];
    CHECK_BEST_UNIQ(rec);
  }

#undef CHECK_BEST_UNIQ

  return *best;
}
