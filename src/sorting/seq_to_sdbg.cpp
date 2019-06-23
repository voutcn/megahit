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

#include "seq_to_sdbg.h"

#include <omp.h>
#include <sequence/readers/edge_io.h>
#include <mutex>
#include <string>
#include <vector>

#include "sequence/kmer.h"
#include "sequence/packed_reads.h"
#include "sequence/readers/async_sequence_reader.h"
#include "sequence/readers/edge_reader.h"
#include "utils/utils.h"

namespace {

/**
 * @brief encode seq_id and its offset in one int64_t
 */
inline int64_t EncodeEdgeOffset(int64_t seq_id, int offset, int strand, SeqPackage &p) {
  return ((p.StartPos(seq_id) + offset) << 1) | strand;
}

inline bool IsDiffKMinusOneMer(uint32_t *item1, uint32_t *item2, int64_t spacing, int k) {
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

inline int Extract_a(uint32_t *item, int num_words, int64_t spacing, unsigned k) {
  int non_dollar = (item[(num_words - 1) * spacing] >> (SeqToSdbg::kBWTCharNumBits + kBitsPerMul)) & 1;

  if (non_dollar) {
    unsigned which_word = (k - 1) / kCharsPerEdgeWord;
    unsigned word_index = (k - 1) % kCharsPerEdgeWord;
    return (item[which_word * spacing] >> (kCharsPerEdgeWord - 1 - word_index) * kBitsPerEdgeChar) & kEdgeCharMask;
  } else {
    return SeqToSdbg::kSentinelValue;
  }
}

inline int Extract_b(uint32_t *item, int num_words, int64_t spacing) {
  return (item[(num_words - 1) * spacing] >> kBitsPerMul) & ((1 << SeqToSdbg::kBWTCharNumBits) - 1);
}

inline int ExtractCounting(uint32_t *item, int num_words, int64_t spacing) {
  return item[(num_words - 1) * spacing] & kMaxMul;
}

/**
 * @brief build lkt for faster binary search for mercy
 */
void InitLookupTable(int64_t *lookup_table, SeqPackage &p) {
  memset(lookup_table, 0xFF, sizeof(int64_t) * SeqToSdbg::kLookUpSize * 2);

  if (p.Size() == 0) {
    return;
  }

  Kmer<1, uint32_t> kmer;
  kmer.InitFromPtr(p.WordPtrAndOffset(0).first, 0, 16);

  uint32_t cur_prefix = kmer.data()[0] >> SeqToSdbg::kLookUpShift;
  lookup_table[cur_prefix * 2] = 0;

  for (int64_t i = 1, num_edges = p.Size(); i < num_edges; ++i) {
    auto ptr_and_offset = p.WordPtrAndOffset(i);
    kmer.InitFromPtr(ptr_and_offset.first, ptr_and_offset.second, 16);

    if ((kmer.data()[0] >> SeqToSdbg::kLookUpShift) > cur_prefix) {
      lookup_table[cur_prefix * 2 + 1] = i - 1;
      cur_prefix = kmer.data()[0] >> SeqToSdbg::kLookUpShift;
      lookup_table[cur_prefix * 2] = i;
    } else {
      assert(cur_prefix == (kmer.data()[0] >> SeqToSdbg::kLookUpShift));
    }
  }

  lookup_table[cur_prefix * 2 + 1] = p.Size() - 1;
}

/**
 * @brief search mercy kmer
 */
int64_t BinarySearchKmer(GenericKmer &kmer, int64_t *lookup_table, SeqPackage &p, int kmer_size) {
  // --- first look up ---
  int64_t l = lookup_table[(kmer.data()[0] >> SeqToSdbg::kLookUpShift) * 2];

  if (l == -1) {
    return -1;
  }

  int64_t r = lookup_table[(kmer.data()[0] >> SeqToSdbg::kLookUpShift) * 2 + 1];
  GenericKmer mid_kmer;

  while (l <= r) {
    int64_t mid = (l + r) / 2;
    auto ptr_and_offset = p.WordPtrAndOffset(mid);
    mid_kmer.InitFromPtr(ptr_and_offset.first, ptr_and_offset.second, kmer_size);
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

} // namespace

// sorting core functions
int64_t SeqToSdbg::Lv0EncodeDiffBase(int64_t read_id) {
  assert(read_id < (int64_t) seq_pkg_.Size());
  return EncodeEdgeOffset(read_id, 0, 0, seq_pkg_);
}

void SeqToSdbg::GenMercyEdges() {
  std::vector<int64_t> edge_lookup(kLookUpSize * 2);
  InitLookupTable(edge_lookup.data(), seq_pkg_);

  std::vector<GenericKmer> mercy_edges;
  std::mutex mercy_lock;
  AsyncReadReader reader(opt_.input_prefix + ".cand");

  int num_threads = std::max(1, opt_.n_threads - 1);
  omp_set_num_threads(num_threads);

  int64_t num_mercy_edges = 0;
  int64_t num_mercy_reads = 0;

  while (true) {
    SeqPackage &rp = reader.Next();
    if (rp.Size() == 0) {
      break;
    }
    xinfo("Read %u reads to search for mercy k-mers\n", rp.Size());

    num_mercy_reads += rp.Size();
    mercy_edges.clear();
#pragma omp parallel for reduction(+ : num_mercy_edges)
    for (unsigned read_id = 0; read_id < rp.Size(); ++read_id) {
      int read_len = rp.SequenceLength(read_id);

      if (read_len < opt_.k + 2) {
        continue;
      }

      std::vector<bool> has_in, has_out;
      GenericKmer kmer, rev_kmer;

      has_in.resize(read_len);
      has_out.resize(read_len);
      std::fill(has_in.begin(), has_in.end(), false);
      std::fill(has_out.begin(), has_out.end(), false);

      auto ptr_and_offset = rp.WordPtrAndOffset(read_id);
      kmer.InitFromPtr(ptr_and_offset.first, ptr_and_offset.second, opt_.k);
      rev_kmer = kmer;
      rev_kmer.ReverseComplement(opt_.k);

      // mark those positions with in/out
      for (int i = 0; i + opt_.k <= read_len; ++i) {
        if (!has_in[i]) {
          // search rc
          if (BinarySearchKmer(rev_kmer, edge_lookup.data(), seq_pkg_, opt_.k) != -1) {
            has_in[i] = true;
          } else {
            // left append ACGT to kmer, if the (k+1)-mer exist, the kmer has in
            rev_kmer.SetBase(opt_.k,
                             3);  // rev kmer is used to compare to kmer, if it's smaller, kmer
            // would not exist in the table
            kmer.ShiftPreappend(0, opt_.k + 1);

            for (int c = 0; c < 4; ++c) {
              kmer.SetBase(0, c);

              if (kmer.cmp(rev_kmer, opt_.k + 1) > 0) {
                break;
              }

              if (BinarySearchKmer(kmer, edge_lookup.data(), seq_pkg_, opt_.k + 1) != -1) {
                has_in[i] = true;
                break;
              }
            }

            rev_kmer.SetBase(opt_.k, 0);
            kmer.ShiftAppend(0, opt_.k + 1);  // clean the k+1-th char
          }
        }

        // check whether has out
        int64_t edge_id = BinarySearchKmer(kmer, edge_lookup.data(), seq_pkg_, opt_.k);

        if (edge_id != -1) {
          has_out[i] = true;

          // BWT see whether the next has in too
          if (i + opt_.k < read_len &&
              seq_pkg_.GetBase(edge_id, opt_.k) == rp.GetBase(read_id, i + opt_.k)) {
            has_in[i + 1] = true;
          }
        } else {
          // search the rc
          kmer.SetBase(opt_.k, 3);
          int next_char = i + opt_.k < read_len ? 3 - rp.GetBase(read_id, i + opt_.k) : 0;
          rev_kmer.ShiftPreappend(next_char, opt_.k + 1);

          if (rev_kmer.cmp(kmer, opt_.k + 1) <= 0 &&
              BinarySearchKmer(rev_kmer, edge_lookup.data(), seq_pkg_, opt_.k + 1) != -1) {
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

              if (BinarySearchKmer(rev_kmer, edge_lookup.data(), seq_pkg_, opt_.k + 1) != -1) {
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
          int next_char = rp.GetBase(read_id, i + opt_.k);
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
                std::lock_guard<std::mutex> lk(mercy_lock);
                auto ptr_and_offset = rp.WordPtrAndOffset(read_id);
                mercy_edges.emplace_back(ptr_and_offset.first, ptr_and_offset.second + j, opt_.k + 1);
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

    for (auto &mercy_edge : mercy_edges) {
      seq_pkg_.AppendCompactSequence(mercy_edge.data(), opt_.k + 1);
    }
  }

  multiplicity.insert(multiplicity.end(), num_mercy_edges, 1);
  xinfo("Number of reads: %ld, Number of mercy edges: %ld\n", num_mercy_reads, num_mercy_edges);
}

SeqToSdbg::Meta SeqToSdbg::Initialize() {
  // reserve space
  {
    long long bases_to_reserve = 0;
    long long num_contigs_to_reserve = 0;
    long long num_multiplicities_to_reserve = 0;

    if (!opt_.input_prefix.empty()) {
      MegahitEdgeReader edge_reader;
      edge_reader.SetFilePrefix(opt_.input_prefix);
      edge_reader.ReadInfo();
      int64_t num_edges = edge_reader.num_edges();
      xinfo("Number edges: %lld\n", (long long) num_edges);

      if (opt_.need_mercy) {
        num_edges += num_edges >> 2;  // it is rare that # mercy > 25%
      }

      bases_to_reserve += num_edges * (edge_reader.k() + 1);
      num_multiplicities_to_reserve += num_edges;
    }

    if (!opt_.contig.empty()) {
      long long num_contigs, num_bases;
      FILE *contig_info = xfopen((opt_.contig + ".info").c_str(), "r");
      if (fscanf(contig_info, "%lld%lld", &num_contigs, &num_bases) != 2) {
        xfatal("Invalid format\n");
      }
      bases_to_reserve += num_bases;
      num_contigs_to_reserve += num_contigs;
      num_multiplicities_to_reserve += num_contigs;
      fclose(contig_info);
    }

    if (!opt_.addi_contig.empty()) {
      long long num_contigs, num_bases;
      FILE *contig_info = xfopen((opt_.addi_contig + ".info").c_str(), "r");
      if (fscanf(contig_info, "%lld%lld", &num_contigs, &num_bases) != 2) {
        xfatal("Invalid format\n");
      }
      bases_to_reserve += num_bases;
      num_contigs_to_reserve += num_contigs;
      num_multiplicities_to_reserve += num_contigs;
      fclose(contig_info);
    }

    if (!opt_.local_contig.empty()) {
      long long num_contigs, num_bases;
      FILE *contig_info = xfopen((opt_.local_contig + ".info").c_str(), "r");
      if (fscanf(contig_info, "%lld%lld", &num_contigs, &num_bases) != 2) {
        xfatal("Invalid format\n");
      }
      bases_to_reserve += num_bases;
      num_contigs_to_reserve += num_contigs;
      num_multiplicities_to_reserve += num_contigs;
      fclose(contig_info);
    }

    xinfo("Bases to reserve: %lld, number contigs: %lld, number multiplicity: %lld\n", bases_to_reserve,
          num_contigs_to_reserve, num_multiplicities_to_reserve);
    seq_pkg_.ReserveSequences(num_contigs_to_reserve);
    seq_pkg_.ReserveBases(bases_to_reserve);
    multiplicity.reserve(num_multiplicities_to_reserve);
  }

  xinfo("Before reading, sizeof seq_package: %lld, multiplicity vector: %lld\n", seq_pkg_.SizeInByte(),
        multiplicity.capacity());

  if (!opt_.input_prefix.empty()) {
    EdgeReader reader(opt_.input_prefix);
    if (opt_.need_mercy) {
      reader.ReadSorted(&seq_pkg_, &multiplicity, 1LL << 60);
    } else {
      reader.ReadUnsorted(&seq_pkg_, &multiplicity, 1LL << 60);
    }
  }

  if (opt_.need_mercy) {
    SimpleTimer timer;
    timer.reset();
    timer.start();
    xinfo("Adding mercy edges...\n");

    GenMercyEdges();
    timer.stop();
    xinfo("Done. Time elapsed: %.4lf\n", timer.elapsed());
  }

  if (!opt_.contig.empty()) {
    ContigReader reader(opt_.contig);
    reader.SetExtendLoop(opt_.k_from, opt_.k)->SetMinLen(opt_.k + 1);
    bool contig_reverse = true;
    reader.ReadAllWithMultiplicity(&seq_pkg_, &multiplicity, contig_reverse);

    // read bubble
    ContigReader bubble_reader({opt_.bubble_seq});
    bubble_reader.SetExtendLoop(opt_.k_from, opt_.k)->SetMinLen(opt_.k + 1);
    bubble_reader.ReadAllWithMultiplicity(&seq_pkg_, &multiplicity, contig_reverse);
  }

  if (!opt_.addi_contig.empty()) {
    ContigReader reader({opt_.addi_contig});
    reader.SetExtendLoop(opt_.k_from, opt_.k)->SetMinLen(opt_.k + 1);
    bool contig_reverse = true;
    reader.ReadAllWithMultiplicity(&seq_pkg_, &multiplicity, contig_reverse);
  }

  if (!opt_.local_contig.empty()) {
    ContigReader reader({opt_.local_contig});
    reader.SetExtendLoop(opt_.k_from, opt_.k)->SetMinLen(opt_.k + 1);
    bool contig_reverse = true;
    reader.ReadAllWithMultiplicity(&seq_pkg_, &multiplicity, contig_reverse);
  }

  xinfo("After reading, sizeof seq_package: %lld, multiplicity vector: %lld\n", seq_pkg_.SizeInByte(),
        multiplicity.capacity());

  seq_pkg_.BuildIndex();
  words_per_substr_ =
      DivCeiling(opt_.k * kBitsPerEdgeChar + kBWTCharNumBits + 1 + kBitsPerMul, kBitsPerEdgeWord);
  words_per_dummy_node_ = DivCeiling(opt_.k * kBitsPerEdgeChar, kBitsPerEdgeWord);

  sdbg_writer_.set_num_threads(opt_.n_threads);
  sdbg_writer_.set_kmer_size(opt_.k);
  sdbg_writer_.set_num_buckets(kNumBuckets);
  sdbg_writer_.set_file_prefix(opt_.output_prefix);
  sdbg_writer_.InitFiles();

  return {
      static_cast<int64_t>(seq_pkg_.Size()),
      static_cast<int64_t>(seq_pkg_.SizeInByte() + multiplicity.size() * sizeof(mul_t)),
      words_per_substr_,
      0,
  };
}

void SeqToSdbg::Lv0CalcBucketSize(ReadPartition *_data) {
  auto &rp = *_data;
  auto &bucket_sizes = rp.rp_bucket_sizes;
  std::fill(bucket_sizes.begin(), bucket_sizes.end(), 0);

  for (int64_t seq_id = rp.rp_start_id; seq_id < rp.rp_end_id; ++seq_id) {
    int seq_len = seq_pkg_.SequenceLength(seq_id);

    if (seq_len < opt_.k + 1) {
      continue;
    }

    uint32_t key = 0;  // $$$$$$$$

    // build initial partial key
    for (int i = 0; i < static_cast<int>(kBucketPrefixLength) - 1; ++i) {
      key = key * kBucketBase + seq_pkg_.GetBase(seq_id, i);
    }

    // sequence = xxxxxxxxx
    // edges = $xxxx, xxxxx, ..., xxxx$
    for (int i = kBucketPrefixLength - 1;
         i - (static_cast<int>(kBucketPrefixLength) - 1) + opt_.k - 1 <= seq_len; ++i) {
      key = (key * kBucketBase + seq_pkg_.GetBase(seq_id, i)) % kNumBuckets;
      bucket_sizes[key]++;
    }

    // reverse complement
    key = 0;

    for (int i = 0; i < static_cast<int>(kBucketPrefixLength) - 1; ++i) {
      key = key * kBucketBase + (3 - seq_pkg_.GetBase(seq_id, seq_len - 1 - i));  // complement
    }

    for (int i = kBucketPrefixLength - 1;
         i - (static_cast<int>(kBucketPrefixLength) - 1) + opt_.k - 1 <= seq_len; ++i) {
      key = key * kBucketBase + (3 - seq_pkg_.GetBase(seq_id, seq_len - 1 - i));
      key %= kNumBuckets;
      bucket_sizes[key]++;
    }
  }
}

void SeqToSdbg::Lv1FillOffsets(ReadPartition *_data) {
  auto &rp = *_data;
  std::array<int64_t, kNumBuckets> prev_full_offsets;

  for (auto b = GetLv1StartBucket(); b < GetLv1EndBucket(); ++b)
    prev_full_offsets[b] = rp.rp_lv1_differential_base;

  // this loop is VERY similar to that in PreprocessScanToFillBucketSizesThread

  // ===== this is a macro to save some copy&paste ================
#define CHECK_AND_SAVE_OFFSET(key, offset, strand)                                                    \
  do {                                                                                                \
    if (HandlingBucket(key)) {                                                           \
      int key_ = GetBucketRank(key);                                                     \
      int64_t full_offset = EncodeEdgeOffset(seq_id, offset, strand, seq_pkg_);                \
      int64_t differential = full_offset - prev_full_offsets[key_];                                   \
      int64_t index = rp.rp_bucket_offsets[key_]++;                                                   \
      WriteOffset(index, differential, full_offset);                                     \
      assert(rp.rp_bucket_offsets[key_] <= GetLv1NumItems());                            \
      prev_full_offsets[key_] = full_offset;                                                          \
    }                                                                                                 \
  } while (0)
  // ^^^^^ why is the macro surrounded by a do-while? please ask Google
  // =========== end macro ==========================

  for (int64_t seq_id = rp.rp_start_id; seq_id < rp.rp_end_id; ++seq_id) {
    int seq_len = seq_pkg_.SequenceLength(seq_id);

    if (seq_len < opt_.k + 1) {
      continue;
    }

    // build initial partial key
    Kmer<1, uint32_t> kmer, rev_kmer;
    auto ptr_and_offset = seq_pkg_.WordPtrAndOffset(seq_id);
    kmer.InitFromPtr(ptr_and_offset.first, ptr_and_offset.second, kBucketPrefixLength);
    auto rev_ptr_and_offset = seq_pkg_.WordPtrAndOffset(seq_id, seq_len - kBucketPrefixLength);
    rev_kmer.InitFromPtr(rev_ptr_and_offset.first, rev_ptr_and_offset.second, kBucketPrefixLength);
    rev_kmer.ReverseComplement(kBucketPrefixLength);

    int key = kmer.data()[0] >> (32 - kBucketPrefixLength * 2);
    int rev_key = rev_kmer.data()[0] >> (32 - kBucketPrefixLength * 2);
    CHECK_AND_SAVE_OFFSET(key, 0, 0);
    CHECK_AND_SAVE_OFFSET(rev_key, 0, 1);

    // sequence = xxxxxxxxx
    // edges = $xxxx, xxxxx, ..., xxxx$
    for (int i = kBucketPrefixLength;
         i - (static_cast<int>(kBucketPrefixLength) - 1) + opt_.k - 1 <= seq_len; ++i) {
      key = (key * kBucketBase + seq_pkg_.GetBase(seq_id, i)) % kNumBuckets;
      rev_key = rev_key * kBucketBase + (3 - seq_pkg_.GetBase(seq_id, seq_len - 1 - i));
      rev_key %= kNumBuckets;
      CHECK_AND_SAVE_OFFSET(key, i - kBucketPrefixLength + 1, 0);
      CHECK_AND_SAVE_OFFSET(rev_key, i - kBucketPrefixLength + 1, 1);
    }
  }
#undef CHECK_AND_SAVE_OFFSET
}

void SeqToSdbg::Lv2ExtractSubString(unsigned from_bucket,
                                    unsigned to_bucket, uint32_t *substr) {
  auto lv1_p = GetLv1Iterator(from_bucket);

  for (auto bucket = from_bucket; bucket < to_bucket; ++bucket) {
    for (int t = 0; t < opt_.n_threads; ++t) {
      int64_t full_offset = GetReadPartition(t).rp_lv1_differential_base;
      int64_t num = GetReadPartition(t).rp_bucket_sizes[bucket];

      for (int64_t i = 0; i < num; ++i) {
        if (*lv1_p >= 0) {
          full_offset += *(lv1_p++);
        } else {
          full_offset = GetSpecialOffset(-1 - *(lv1_p++));
        }

        int64_t seq_id = seq_pkg_.GetSeqID(full_offset >> 1);
        int offset = (full_offset >> 1) - seq_pkg_.StartPos(seq_id);
        int strand = full_offset & 1;

        int seq_len = seq_pkg_.SequenceLength(seq_id);
        int num_chars_to_copy = opt_.k - (offset + opt_.k > seq_len);
        int counting = 0;

        if (offset > 0 && offset + opt_.k <= seq_len) {
          counting = multiplicity[seq_id];
        }

        auto ptr_and_offset = seq_pkg_.WordPtrAndOffset(seq_id);
        int start_offset = ptr_and_offset.second;
        int words_this_seq = DivCeiling(start_offset + seq_len, 16);
        const uint32_t *edge_p = ptr_and_offset.first;

        if (strand == 0) {
          // copy counting and W char
          int prev_char;

          if (offset == 0) {
            assert(num_chars_to_copy == opt_.k);
            prev_char = kSentinelValue;
          } else {
            prev_char = seq_pkg_.GetBase(seq_id, offset - 1);
          }

          CopySubstring(substr, edge_p, offset + start_offset, num_chars_to_copy, 1, words_this_seq,
                        words_per_substr_);

          uint32_t *last_word = substr + words_per_substr_ - 1;
          *last_word |= int(num_chars_to_copy == opt_.k) << (kBWTCharNumBits + kBitsPerMul);
          *last_word |= prev_char << kBitsPerMul;
          *last_word |= std::max(0, kMaxMul - counting);  // then larger counting come first after sorting
        } else {
          int prev_char;

          if (offset == 0) {
            assert(num_chars_to_copy == opt_.k);
            prev_char = kSentinelValue;
          } else {
            prev_char = 3 - seq_pkg_.GetBase(seq_id, seq_len - 1 - offset + 1);
          }

          offset = seq_len - 1 - offset - (opt_.k - 1);  // switch to normal strand

          if (offset < 0) {
            assert(num_chars_to_copy == opt_.k - 1);
            offset = 0;
          }

          CopySubstringRC(substr, edge_p, offset + start_offset, num_chars_to_copy, 1, words_this_seq,
                          words_per_substr_);

          uint32_t *last_word = substr + words_per_substr_ - 1;
          *last_word |= int(num_chars_to_copy == opt_.k) << (kBWTCharNumBits + kBitsPerMul);
          *last_word |= prev_char << kBitsPerMul;
          *last_word |= std::max(0, kMaxMul - counting);
        }

        substr += words_per_substr_;
      }
    }
  }
}

void SeqToSdbg::Lv2Postprocess(int64_t from, int64_t to, int tid, uint32_t *substr) {
  int64_t start_idx, end_idx;
  int has_solid_a = 0;  // has solid (k+1)-mer aSb
  int has_solid_b = 0;  // has solid aSb
  int64_t last_a[4], outputed_b;
  uint32_t tip_label[32];
  SdbgWriter::Snapshot snapshot;

  for (start_idx = from; start_idx < to; start_idx = end_idx) {
    end_idx = start_idx + 1;
    uint32_t *item = substr + start_idx * words_per_substr_;

    while (end_idx < to &&
        !IsDiffKMinusOneMer(item, substr + end_idx * words_per_substr_, 1, opt_.k)) {
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

      if (a != kSentinelValue && (b != kSentinelValue || !(has_solid_a & (1 << a)))) {
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

      if (is_dollar) {
        for (int64_t i = 0; i < words_per_dummy_node_; ++i) {
          tip_label[i] = cur_item[i];
        }
      }

      sdbg_writer_.Write(tid, cur_item[0] >> (32 - kBucketPrefixLength * 2), w, last, is_dollar,
                         kMaxMul - ExtractCounting(cur_item, words_per_substr_, 1), tip_label,
                         &snapshot);
    }
  }
  sdbg_writer_.SaveSnapshot(snapshot);
}

void SeqToSdbg::Lv0Postprocess() {
  sdbg_writer_.Finalize();
  xinfo("Number of $ A C G T A- C- G- T-:\n");
  xinfo("");
  for (int i = 0; i < 9; ++i) {
    xinfoc("%lld ", (long long) sdbg_writer_.final_meta().w_count(i));
  }

  xinfoc("\n");
  xinfo("Total number of edges: %lld\n", (long long) sdbg_writer_.final_meta().item_count());
  xinfo("Total number of ONEs: %lld\n", (long long) sdbg_writer_.final_meta().ones_in_last());
  xinfo("Total number of $v edges: %lld\n", (long long) sdbg_writer_.final_meta().tip_count());

  assert(sdbg_writer_.final_meta().w_count(0) == sdbg_writer_.final_meta().tip_count());
}
