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

#include "seq_to_sdbg.h"

#include <omp.h>
#include <deque>
#include <string>
#include <vector>

#include "sequence/copy_substr.h"
#include "sequence/io/async_sequence_reader.h"
#include "sequence/io/edge/edge_reader.h"
#include "sequence/kmer.h"
#include "utils/mutex.h"
#include "utils/utils.h"

namespace {

/**
 * @brief encode seq_id and its offset in one int64_t
 */
inline int64_t EncodeEdgeOffset(int64_t seq_id, int offset, int strand,
                                const SeqPackage &p) {
  return ((p.GetSeqView(seq_id).full_offset_in_pkg() + offset) << 1) | strand;
}

inline bool IsDiffKMinusOneMer(uint32_t *item1, uint32_t *item2,
                               int64_t spacing, int k) {
  // mask extra bits
  int chars_in_last_word = (k - 1) % kCharsPerEdgeWord;
  int num_full_words = (k - 1) / kCharsPerEdgeWord;

  if (chars_in_last_word > 0) {
    uint32_t w1 = item1[num_full_words * spacing];
    uint32_t w2 = item2[num_full_words * spacing];

    if ((w1 >> (kCharsPerEdgeWord - chars_in_last_word) * kBitsPerEdgeChar) !=
        (w2 >> (kCharsPerEdgeWord - chars_in_last_word) * kBitsPerEdgeChar)) {
      return true;
    }
  }

  for (int i = num_full_words - 1; i >= 0; --i) {
    if (item1[i * spacing] != item2[i * spacing]) {
      return true;
    }
  }

  return false;
}

inline int Extract_a(uint32_t *item, int num_words, int64_t spacing,
                     unsigned k) {
  int non_dollar = (item[(num_words - 1) * spacing] >>
                    (SeqToSdbg::kBWTCharNumBits + kBitsPerMul)) &
                   1;

  if (non_dollar) {
    unsigned which_word = (k - 1) / kCharsPerEdgeWord;
    unsigned word_index = (k - 1) % kCharsPerEdgeWord;
    return (item[which_word * spacing] >>
            (kCharsPerEdgeWord - 1 - word_index) * kBitsPerEdgeChar) &
           kEdgeCharMask;
  } else {
    return SeqToSdbg::kSentinelValue;
  }
}

inline int Extract_b(uint32_t *item, int num_words, int64_t spacing) {
  return (item[(num_words - 1) * spacing] >> kBitsPerMul) &
         ((1 << SeqToSdbg::kBWTCharNumBits) - 1);
}

inline int ExtractCounting(uint32_t *item, int num_words, int64_t spacing) {
  return item[(num_words - 1) * spacing] & kMaxMul;
}

/**
 * @brief build lkt for faster binary search for mercy
 */
void InitLookupTable(int64_t *lookup_table, SeqPackage &p) {
  memset(lookup_table, 0xFF, sizeof(int64_t) * SeqToSdbg::kLookUpSize * 2);

  if (p.seq_count() == 0) {
    return;
  }

  Kmer<1, uint32_t> kmer;
  kmer.InitFromPtr(p.GetSeqView(0).raw_address().first, 0, 16);

  uint32_t cur_prefix = kmer.data()[0] >> SeqToSdbg::kLookUpShift;
  lookup_table[cur_prefix * 2] = 0;

  for (int64_t i = 1, num_edges = p.seq_count(); i < num_edges; ++i) {
    auto ptr_and_offset = p.GetSeqView(i).raw_address();
    kmer.InitFromPtr(ptr_and_offset.first, ptr_and_offset.second, 16);

    if ((kmer.data()[0] >> SeqToSdbg::kLookUpShift) > cur_prefix) {
      lookup_table[cur_prefix * 2 + 1] = i - 1;
      cur_prefix = kmer.data()[0] >> SeqToSdbg::kLookUpShift;
      lookup_table[cur_prefix * 2] = i;
    } else {
      assert(cur_prefix == (kmer.data()[0] >> SeqToSdbg::kLookUpShift));
    }
  }

  lookup_table[cur_prefix * 2 + 1] = p.seq_count() - 1;
}

/**
 * @brief search mercy kmer
 */
int64_t BinarySearchKmer(GenericKmer &kmer, int64_t *lookup_table,
                         SeqPackage &p, int kmer_size) {
  // --- first look up ---
  int64_t l = lookup_table[(kmer.data()[0] >> SeqToSdbg::kLookUpShift) * 2];

  if (l == -1) {
    return -1;
  }

  int64_t r = lookup_table[(kmer.data()[0] >> SeqToSdbg::kLookUpShift) * 2 + 1];
  GenericKmer mid_kmer;

  while (l <= r) {
    int64_t mid = (l + r) / 2;
    auto ptr_and_offset = p.GetSeqView(mid).raw_address();
    mid_kmer.InitFromPtr(ptr_and_offset.first, ptr_and_offset.second,
                         kmer_size);
    int cmp = kmer.cmp(mid_kmer, kmer_size);

    if (cmp > 0) {
      l = mid + 1;
    } else if (cmp < 0) {
      r = mid - 1;
    } else {
      return mid;
    }
  }

  return -1;
}

}  // namespace

// sorting core functions
int64_t SeqToSdbg::Lv0EncodeDiffBase(int64_t read_id) {
  assert(read_id < (int64_t)seq_pkg_.seq_count());
  return EncodeEdgeOffset(read_id, 0, 0, seq_pkg_);
}

void SeqToSdbg::GenMercyEdges() {
  std::vector<int64_t> edge_lookup(kLookUpSize * 2);
  InitLookupTable(edge_lookup.data(), seq_pkg_);

  BinaryReader binary_reader(opt_.input_prefix + ".cand");
  AsyncSequenceReader reader(&binary_reader, false, 1 << 16, 1 << 23);
  std::deque<GenericKmer> mercy_edges;
  SpinLock mercy_lock;

  int num_threads = std::max(1, opt_.n_threads - 1);
  omp_set_num_threads(num_threads);

  int64_t num_mercy_edges = 0;
  int64_t num_mercy_reads = 0;

  while (true) {
    const auto &batch_reads = reader.Next();
    if (batch_reads.seq_count() == 0) {
      break;
    }
    xinfo("Read {} reads to search for mercy k-mers\n",
          batch_reads.seq_count());

    num_mercy_reads += batch_reads.seq_count();
    mercy_edges.clear();

#pragma omp parallel for reduction(+ : num_mercy_edges)
    for (unsigned read_id = 0; read_id < batch_reads.seq_count(); ++read_id) {
      auto seq_view = batch_reads.GetSeqView(read_id);
      unsigned read_len = seq_view.length();

      if (read_len < opt_.k + 2) {
        continue;
      }

      std::vector<bool> has_in, has_out;
      GenericKmer kmer, rev_kmer;

      has_in.resize(read_len);
      has_out.resize(read_len);
      std::fill(has_in.begin(), has_in.end(), false);
      std::fill(has_out.begin(), has_out.end(), false);

      auto ptr_and_offset = seq_view.raw_address();
      kmer.InitFromPtr(ptr_and_offset.first, ptr_and_offset.second, opt_.k);
      rev_kmer = kmer;
      rev_kmer.ReverseComplement(opt_.k);

      // mark those positions with in/out
      for (int i = 0; i + opt_.k <= read_len; ++i) {
        if (!has_in[i]) {
          // search rc
          if (BinarySearchKmer(rev_kmer, edge_lookup.data(), seq_pkg_,
                               opt_.k) != -1) {
            has_in[i] = true;
          } else {
            // left append ACGT to kmer, if the (k+1)-mer exist, the kmer has in
            rev_kmer.SetBase(opt_.k, 3);  // rev kmer is used to compare to
                                          // kmer, if it's smaller, kmer
            // would not exist in the table
            kmer.ShiftPreappend(0, opt_.k + 1);

            for (int c = 0; c < 4; ++c) {
              kmer.SetBase(0, c);

              if (kmer.cmp(rev_kmer, opt_.k + 1) > 0) {
                break;
              }

              if (BinarySearchKmer(kmer, edge_lookup.data(), seq_pkg_,
                                   opt_.k + 1) != -1) {
                has_in[i] = true;
                break;
              }
            }

            rev_kmer.SetBase(opt_.k, 0);
            kmer.ShiftAppend(0, opt_.k + 1);  // clean the k+1-th char
          }
        }

        // check whether has out
        int64_t edge_id =
            BinarySearchKmer(kmer, edge_lookup.data(), seq_pkg_, opt_.k);

        if (edge_id != -1) {
          has_out[i] = true;

          // BWT see whether the next has in too
          if (i + opt_.k < read_len &&
              seq_pkg_.GetSeqView(edge_id).base_at(opt_.k) ==
                  seq_view.base_at(i + opt_.k)) {
            has_in[i + 1] = true;
          }
        } else {
          // search the rc
          kmer.SetBase(opt_.k, 3);
          int next_char =
              i + opt_.k < read_len ? 3 - seq_view.base_at(i + opt_.k) : 0;
          rev_kmer.ShiftPreappend(next_char, opt_.k + 1);

          if (rev_kmer.cmp(kmer, opt_.k + 1) <= 0 &&
              BinarySearchKmer(rev_kmer, edge_lookup.data(), seq_pkg_,
                               opt_.k + 1) != -1) {
            has_out[i] = true;
            has_in[i + 1] = true;
          } else {
            for (int c = 0; c < 4; ++c) {
              if (c == next_char) {
                continue;
              }

              rev_kmer.SetBase(0, c);

              if (rev_kmer.cmp(kmer, opt_.k + 1) > 0) {
                break;
              }

              if (BinarySearchKmer(rev_kmer, edge_lookup.data(), seq_pkg_,
                                   opt_.k + 1) != -1) {
                has_out[i] = true;
                break;
              }
            }
          }

          kmer.SetBase(opt_.k, 0);
          rev_kmer.ShiftAppend(0, opt_.k + 1);
        }

        // shift kmer and rev_kmer
        if (i + opt_.k < read_len) {
          int next_char = seq_view.base_at(i + opt_.k);
          kmer.ShiftAppend(next_char, opt_.k);
          rev_kmer.ShiftPreappend(3 - next_char, opt_.k);
        }
      }

      // adding mercy edges
      int last_no_out = -1;

      for (int i = 0; i + opt_.k <= read_len; ++i) {
        switch (has_in[i] | (int(has_out[i]) << 1)) {
          case 1: {  // has incoming only
            last_no_out = i;
            break;
          }

          case 2: {  // has outgoing only
            if (last_no_out >= 0) {
              for (int j = last_no_out; j < i; ++j) {
                auto raw_address = seq_view.raw_address();
                GenericKmer mercy_edge(raw_address.first,
                                       raw_address.second + j, opt_.k + 1);
                std::lock_guard<SpinLock> lk(mercy_lock);
                mercy_edges.emplace_back(mercy_edge);
              }

              num_mercy_edges += i - last_no_out;
            }

            last_no_out = -1;
            break;
          }

          case 3: {  // has in and out
            last_no_out = -1;
            break;
          }

          default: {
            // do nothing
            break;
          }
        }
      }
    }

    for (const auto &mercy_edge : mercy_edges) {
      seq_pkg_.AppendCompactSequence(mercy_edge.data(), opt_.k + 1);
    }
  }

  multiplicity.insert(multiplicity.end(), num_mercy_edges, 1);
  xinfo("Number of reads: {}, Number of mercy edges: {}\n", num_mercy_reads,
        num_mercy_edges);
}

SeqToSdbg::MemoryStat SeqToSdbg::Initialize() {
  // reserve space
  {
    int64_t bases_to_reserve = 0;
    int64_t num_contigs_to_reserve = 0;
    int64_t num_multiplicities_to_reserve = 0;

    if (!opt_.input_prefix.empty()) {
      EdgeReader edge_reader(opt_.input_prefix);
      int64_t num_edges = edge_reader.GetMetadata().num_edges;
      xinfo("Number edges: {}\n", num_edges);

      if (opt_.need_mercy) {
        auto mercy_factor = std::getenv("MEGAHIT_NUM_MERCY_FACTOR");
        if (mercy_factor) {
          char *_;
          num_edges *= 1 + strtod(mercy_factor, &_);
        } else {
          num_edges += num_edges >> 2;  // it is rare that # mercy > 25%
        }
      }

      bases_to_reserve += num_edges * (edge_reader.GetMetadata().kmer_size + 1);
      num_multiplicities_to_reserve += num_edges;
    }

    if (!opt_.contig.empty()) {
      auto sizes = ContigReader(opt_.contig).GetNumContigsAndBases();
      bases_to_reserve += sizes.second;
      num_contigs_to_reserve += sizes.first;
      num_multiplicities_to_reserve += sizes.first;
    }

    if (!opt_.bubble_seq.empty()) {
      auto sizes = ContigReader(opt_.bubble_seq).GetNumContigsAndBases();
      bases_to_reserve += sizes.second;
      num_contigs_to_reserve += sizes.first;
      num_multiplicities_to_reserve += sizes.first;
    }

    if (!opt_.addi_contig.empty()) {
      auto sizes = ContigReader(opt_.addi_contig).GetNumContigsAndBases();
      bases_to_reserve += sizes.second;
      num_contigs_to_reserve += sizes.first;
      num_multiplicities_to_reserve += sizes.first;
    }

    if (!opt_.local_contig.empty()) {
      auto sizes = ContigReader(opt_.local_contig).GetNumContigsAndBases();
      bases_to_reserve += sizes.second;
      num_contigs_to_reserve += sizes.first;
      num_multiplicities_to_reserve += sizes.first;
    }

    xinfo("Bases to reserve: {}, number contigs: {}, number multiplicity: {}\n",
          bases_to_reserve, num_contigs_to_reserve,
          num_multiplicities_to_reserve);
    seq_pkg_.ReserveSequences(num_contigs_to_reserve);
    seq_pkg_.ReserveBases(bases_to_reserve);
    multiplicity.reserve(num_multiplicities_to_reserve);
  }

  xinfo("Before reading, sizeof seq_package: {}, multiplicity vector: {}\n",
        seq_pkg_.size_in_byte(), multiplicity.capacity());

  if (!opt_.input_prefix.empty()) {
    EdgeReader reader(opt_.input_prefix);
    reader.SetMultiplicityVec(&multiplicity);
    auto n_read = reader.ReadAll(&seq_pkg_, false);
    xinfo("Read {} edges.\n", n_read);
    xinfo(
        "After reading, sizeof seq_package: {}/{}/{}, multiplicity vector: "
        "{}/{}\n",
        seq_pkg_.size_in_byte(), seq_pkg_.seq_count(), seq_pkg_.base_count(),
        multiplicity.size(), multiplicity.capacity());
  }

  if (opt_.need_mercy) {
    SimpleTimer timer;
    timer.reset();
    timer.start();
    xinfo("Adding mercy edges...\n");

    GenMercyEdges();
    timer.stop();
    xinfo("Done. Time elapsed: {.4}\n", timer.elapsed());
    xinfo(
        "After adding mercy, sizeof seq_package: {}/{}/{}, multiplicity "
        "vector: {}/{}\n",
        seq_pkg_.size_in_byte(), seq_pkg_.seq_count(), seq_pkg_.base_count(),
        multiplicity.size(), multiplicity.capacity());
  }

  if (!opt_.contig.empty()) {
    ContigReader reader(opt_.contig);
    reader.SetExtendLoop(opt_.k_from, opt_.k)->SetMinLen(opt_.k + 1);
    bool contig_reverse = true;
    auto n_read = reader.ReadAllWithMultiplicity(&seq_pkg_, &multiplicity,
                                                 contig_reverse);
    xinfo("Read {} contigs from {}.\n", n_read, opt_.contig.c_str());
    xinfo(
        "After reading contigs, sizeof seq_package: {}/{}/{}, multiplicity "
        "vector: {}/{}\n",
        seq_pkg_.size_in_byte(), seq_pkg_.seq_count(), seq_pkg_.base_count(),
        multiplicity.size(), multiplicity.capacity());

    // read bubble
    ContigReader bubble_reader(opt_.bubble_seq);
    bubble_reader.SetMinLen(opt_.k + 1);
    n_read = bubble_reader.ReadAllWithMultiplicity(&seq_pkg_, &multiplicity,
                                                   contig_reverse);
    xinfo("Read {} contigs from {}.\n", n_read, opt_.bubble_seq.c_str());
    xinfo(
        "After reading contigs, sizeof seq_package: {}/{}/{}, multiplicity "
        "vector: {}/{}\n",
        seq_pkg_.size_in_byte(), seq_pkg_.seq_count(), seq_pkg_.base_count(),
        multiplicity.size(), multiplicity.capacity());
  }

  if (!opt_.addi_contig.empty()) {
    ContigReader reader(opt_.addi_contig);
    reader.SetMinLen(opt_.k + 1);
    bool contig_reverse = true;
    auto n_read = reader.ReadAllWithMultiplicity(&seq_pkg_, &multiplicity,
                                                 contig_reverse);
    xinfo("Read {} contigs from {}.\n", n_read, opt_.addi_contig.c_str());
    xinfo(
        "After reading contigs, sizeof seq_package: {}/{}/{}, multiplicity "
        "vector: {}/{}\n",
        seq_pkg_.size_in_byte(), seq_pkg_.seq_count(), seq_pkg_.base_count(),
        multiplicity.size(), multiplicity.capacity());
  }

  if (!opt_.local_contig.empty()) {
    ContigReader reader(opt_.local_contig);
    reader.SetMinLen(opt_.k + 1);
    bool contig_reverse = true;
    auto n_read = reader.ReadAllWithMultiplicity(&seq_pkg_, &multiplicity,
                                                 contig_reverse);
    xinfo("Read {} contigs from {}.\n", n_read, opt_.local_contig.c_str());
    xinfo(
        "After reading contigs, sizeof seq_package: {}/{}/{}, multiplicity "
        "vector: {}/{}\n",
        seq_pkg_.size_in_byte(), seq_pkg_.seq_count(), seq_pkg_.base_count(),
        multiplicity.size(), multiplicity.capacity());
  }

  xinfo("Finally, sizeof seq_package: {}/{}/{}, multiplicity vector: {}/{}\n",
        seq_pkg_.size_in_byte(), seq_pkg_.seq_count(), seq_pkg_.base_count(),
        multiplicity.size(), multiplicity.capacity());

  seq_pkg_.BuildIndex();
  words_per_substr_ =
      DivCeiling(opt_.k * kBitsPerEdgeChar + kBWTCharNumBits + 1 + kBitsPerMul,
                 kBitsPerEdgeWord);

  sdbg_writer_.set_num_threads(opt_.n_threads);
  sdbg_writer_.set_kmer_size(opt_.k);
  sdbg_writer_.set_num_buckets(kNumBuckets);
  sdbg_writer_.set_file_prefix(opt_.output_prefix);
  sdbg_writer_.InitFiles();

  return {
      static_cast<int64_t>(seq_pkg_.seq_count()),
      static_cast<int64_t>(seq_pkg_.size_in_byte() +
                           multiplicity.capacity() * sizeof(mul_t)),
      words_per_substr_,
      0,
  };
}

void SeqToSdbg::Lv0CalcBucketSize(int64_t seq_from, int64_t seq_to,
                                  std::array<int64_t, kNumBuckets> *out) {
  auto &bucket_sizes = *out;
  std::fill(bucket_sizes.begin(), bucket_sizes.end(), 0);

  for (int64_t seq_id = seq_from; seq_id < seq_to; ++seq_id) {
    auto seq_view = seq_pkg_.GetSeqView(seq_id);
    unsigned seq_len = seq_view.length();

    if (seq_len < opt_.k + 1) {
      continue;
    }

    uint32_t key = 0;  // $$$$$$$$

    // build initial partial key
    for (int i = 0; i < static_cast<int>(kBucketPrefixLength) - 1; ++i) {
      key = key * kBucketBase + seq_view.base_at(i);
    }

    // sequence = xxxxxxxxx
    // edges = $xxxx, xxxxx, ..., xxxx$
    for (int i = kBucketPrefixLength - 1;
         i - (static_cast<int>(kBucketPrefixLength) - 1) + opt_.k - 1 <=
         seq_len;
         ++i) {
      key = (key * kBucketBase + seq_view.base_at(i)) % kNumBuckets;
      bucket_sizes[key]++;
    }

    // reverse complement
    key = 0;

    for (int i = 0; i < static_cast<int>(kBucketPrefixLength) - 1; ++i) {
      key = key * kBucketBase +
            (3 - seq_view.base_at(seq_len - 1 - i));  // complement
    }

    for (int i = kBucketPrefixLength - 1;
         i - (static_cast<int>(kBucketPrefixLength) - 1) + opt_.k - 1 <=
         seq_len;
         ++i) {
      key = key * kBucketBase + (3 - seq_view.base_at(seq_len - 1 - i));
      key %= kNumBuckets;
      bucket_sizes[key]++;
    }
  }
}

void SeqToSdbg::Lv1FillOffsets(OffsetFiller &filler, int64_t seq_from,
                               int64_t seq_to) {
// ===== this is a macro to save some copy&paste ================
#define CHECK_AND_SAVE_OFFSET(key, offset, strand)                  \
  do {                                                              \
    if (filler.IsHandling(key)) {                                   \
      filler.WriteNextOffset(                                       \
          key, EncodeEdgeOffset(seq_id, offset, strand, seq_pkg_)); \
    }                                                               \
  } while (0)
  // =========== end macro ==========================

  for (int64_t seq_id = seq_from; seq_id < seq_to; ++seq_id) {
    auto seq_view = seq_pkg_.GetSeqView(seq_id);
    unsigned seq_len = seq_view.length();
    if (seq_len < opt_.k + 1) {
      continue;
    }

    // build initial partial key
    Kmer<1, uint32_t> kmer, rev_kmer;
    auto ptr_and_offset = seq_view.raw_address();
    kmer.InitFromPtr(ptr_and_offset.first, ptr_and_offset.second,
                     kBucketPrefixLength);
    auto rev_ptr_and_offset =
        seq_view.raw_address(seq_len - kBucketPrefixLength);
    rev_kmer.InitFromPtr(rev_ptr_and_offset.first, rev_ptr_and_offset.second,
                         kBucketPrefixLength);
    rev_kmer.ReverseComplement(kBucketPrefixLength);

    int key = kmer.data()[0] >> (32 - kBucketPrefixLength * 2);
    int rev_key = rev_kmer.data()[0] >> (32 - kBucketPrefixLength * 2);
    CHECK_AND_SAVE_OFFSET(key, 0, 0);
    CHECK_AND_SAVE_OFFSET(rev_key, 0, 1);

    // sequence = xxxxxxxxx
    // edges = $xxxx, xxxxx, ..., xxxx$
    for (int i = kBucketPrefixLength;
         i - (static_cast<int>(kBucketPrefixLength) - 1) + opt_.k - 1 <=
         seq_len;
         ++i) {
      key = (key * kBucketBase + seq_view.base_at(i)) % kNumBuckets;
      rev_key = rev_key * kBucketBase + (3 - seq_view.base_at(seq_len - 1 - i));
      rev_key %= kNumBuckets;
      CHECK_AND_SAVE_OFFSET(key, i - kBucketPrefixLength + 1, 0);
      CHECK_AND_SAVE_OFFSET(rev_key, i - kBucketPrefixLength + 1, 1);
    }
  }
#undef CHECK_AND_SAVE_OFFSET
}

void SeqToSdbg::Lv2ExtractSubString(OffsetFetcher &fetcher, SubstrPtr substr) {
  while (fetcher.HasNext()) {
    int64_t full_offset = fetcher.Next();
    auto seq_view = seq_pkg_.GetSeqViewByOffset(full_offset >> 1);
    int offset = (full_offset >> 1) - seq_view.full_offset_in_pkg();
    unsigned strand = full_offset & 1;

    unsigned seq_len = seq_view.length();
    unsigned num_chars_to_copy = opt_.k - (offset + opt_.k > seq_len);
    int counting = 0;

    if (offset > 0 && offset + opt_.k <= seq_len) {
      counting = multiplicity[seq_view.id()];
    }

    auto ptr_and_offset = seq_view.raw_address();
    unsigned start_offset = ptr_and_offset.second;
    unsigned words_this_seq = DivCeiling(start_offset + seq_len, 16);
    const uint32_t *edge_p = ptr_and_offset.first;

    if (strand == 0) {
      // copy counting and W char
      unsigned prev_char;

      if (offset == 0) {
        assert(num_chars_to_copy == opt_.k);
        prev_char = kSentinelValue;
      } else {
        prev_char = seq_view.base_at(offset - 1);
      }

      CopySubstring(substr, edge_p, offset + start_offset, num_chars_to_copy, 1,
                    words_this_seq, words_per_substr_);

      auto last_word = substr + words_per_substr_ - 1;
      *last_word |= unsigned(num_chars_to_copy == opt_.k)
                    << (kBWTCharNumBits + kBitsPerMul);
      *last_word |= prev_char << kBitsPerMul;
      *last_word |= std::max(
          0,
          kMaxMul - counting);  // then larger counting come first after sorting
    } else {
      unsigned prev_char;

      if (offset == 0) {
        assert(num_chars_to_copy == opt_.k);
        prev_char = kSentinelValue;
      } else {
        prev_char = 3 - seq_view.base_at(seq_len - 1 - offset + 1);
      }

      offset = seq_len - 1 - offset - (opt_.k - 1);  // switch to normal strand

      if (offset < 0) {
        assert(num_chars_to_copy == opt_.k - 1);
        offset = 0;
      }

      CopySubstringRC(substr, edge_p, offset + start_offset, num_chars_to_copy,
                      1, words_this_seq, words_per_substr_);

      auto last_word = substr + words_per_substr_ - 1;
      *last_word |= unsigned(num_chars_to_copy == opt_.k)
                    << (kBWTCharNumBits + kBitsPerMul);
      *last_word |= prev_char << kBitsPerMul;
      *last_word |= std::max(0, kMaxMul - counting);
    }

    substr += words_per_substr_;
  }
}

void SeqToSdbg::Lv2Postprocess(int64_t from, int64_t to, int tid,
                               uint32_t *substr) {
  int64_t start_idx, end_idx;
  int has_solid_a = 0;  // has solid (k+1)-mer aSb
  int has_solid_b = 0;  // has solid aSb
  int64_t last_a[4], outputed_b;
  SdbgWriter::Snapshot snapshot;

  for (start_idx = from; start_idx < to; start_idx = end_idx) {
    end_idx = start_idx + 1;
    uint32_t *item = substr + start_idx * words_per_substr_;

    while (end_idx < to &&
           !IsDiffKMinusOneMer(item, substr + end_idx * words_per_substr_, 1,
                               opt_.k)) {
      ++end_idx;
    }

    // clean marking
    has_solid_a = has_solid_b = 0;
    outputed_b = 0;

    for (int64_t i = start_idx; i < end_idx; ++i) {
      uint32_t *cur_item = substr + i * words_per_substr_;
      int a = Extract_a(cur_item, words_per_substr_, 1, opt_.k);
      int b = Extract_b(cur_item, words_per_substr_, 1);

      if (a != kSentinelValue && b != kSentinelValue) {
        has_solid_a |= 1 << a;
        has_solid_b |= 1 << b;
      }

      if (a != kSentinelValue &&
          (b != kSentinelValue || !(has_solid_a & (1 << a)))) {
        last_a[a] = i;
      }
    }

    for (int64_t i = start_idx, j; i < end_idx; i = j) {
      uint32_t *cur_item = substr + i * words_per_substr_;
      int a = Extract_a(cur_item, words_per_substr_, 1, opt_.k);
      int b = Extract_b(cur_item, words_per_substr_, 1);

      j = i + 1;

      while (j < end_idx) {
        uint32_t *next_item = substr + j * words_per_substr_;

        if (Extract_a(next_item, words_per_substr_, 1, opt_.k) != a ||
            Extract_b(next_item, words_per_substr_, 1) != b) {
          break;
        } else {
          ++j;
        }
      }

      int w, last, is_dollar = 0;

      if (a == kSentinelValue) {
        assert(b != kSentinelValue);

        if (has_solid_b & (1 << b)) {
          continue;
        }

        is_dollar = 1;
      }

      if (b == kSentinelValue) {
        assert(a != kSentinelValue);

        if (has_solid_a & (1 << a)) {
          continue;
        }
      }

      w = (b == kSentinelValue) ? 0 : ((outputed_b & (1 << b)) ? b + 5 : b + 1);
      last = (a == kSentinelValue) ? 0 : ((last_a[a] == j - 1) ? 1 : 0);
      outputed_b |= 1 << b;

      sdbg_writer_.Write(
          tid, cur_item[0] >> (32 - kBucketPrefixLength * 2), w, last,
          is_dollar, kMaxMul - ExtractCounting(cur_item, words_per_substr_, 1),
          cur_item, &snapshot);
    }
  }
  sdbg_writer_.SaveSnapshot(snapshot);
}

void SeqToSdbg::Lv0Postprocess() {
  sdbg_writer_.Finalize();
  xinfo("Number of $ A C G T A- C- G- T-:\n");
  xinfo("");
  for (int i = 0; i < 9; ++i) {
    xinfoc("{} ", sdbg_writer_.final_meta().w_count(i));
  }

  xinfoc("{s}", "\n");
  xinfo("Total number of edges: {}\n", sdbg_writer_.final_meta().item_count());
  xinfo("Total number of ONEs: {}\n", sdbg_writer_.final_meta().ones_in_last());
  xinfo("Total number of $v edges: {}\n",
        sdbg_writer_.final_meta().tip_count());

  assert(sdbg_writer_.final_meta().w_count(0) ==
         sdbg_writer_.final_meta().tip_count());
}
