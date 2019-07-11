#include "kmer_counter.h"

#include <algorithm>

#include "sequence/copy_substr.h"
#include "sequence/io/sequence_lib.h"
#include "sequence/kmer.h"
#include "utils/utils.h"

/**
 * @brief encode read_id and its offset in one int64_t
 */
inline int64_t EncodeOffset(int64_t read_id, int offset, int strand,
                            const SeqPackage &p) {
  return ((p.GetSeqView(read_id).full_offset_in_pkg() + offset) << 1) | strand;
}

inline bool IsDifferentEdges(uint32_t *item1, uint32_t *item2, int num_words,
                             int64_t spacing) {
  for (int i = num_words - 1; i >= 0; --i) {
    if (*(item1 + i * spacing) != *(item2 + i * spacing)) {
      return true;
    }
  }

  return false;
}

/**
 * @brief pack an edge and its multiplicity to word-aligned spaces
 */
void KmerCounter::PackEdge(uint32_t *dest, uint32_t *item, int64_t counting) {
  for (int i = 0; i < words_per_edge_ && i < words_per_substr_; ++i) {
    dest[i] = *(item + i);
  }

  int chars_in_last_word = (opt_.k + 1) % kCharsPerEdgeWord;
  int which_word = (opt_.k + 1) / kCharsPerEdgeWord;

  if (chars_in_last_word > 0) {
    dest[which_word] >>= (kCharsPerEdgeWord - chars_in_last_word) * 2;
    dest[which_word] <<= (kCharsPerEdgeWord - chars_in_last_word) * 2;
  } else {
    dest[which_word] = 0;
  }

  while (++which_word < words_per_edge_) {
    dest[which_word] = 0;
  }

  dest[words_per_edge_ - 1] |= std::min(int64_t(kMaxMul), counting);
}

// function pass to BaseSequenceSortingEngine

int64_t KmerCounter::Lv0EncodeDiffBase(int64_t read_id) {
  return EncodeOffset(read_id, 0, 0, seq_pkg_);
}

KmerCounter::MemoryStat KmerCounter::Initialize() {
  bool is_reverse = true;

  int64_t num_bases, num_reads;
  SequenceLibCollection seq_collection(opt_.read_lib_file);
  auto collection_size = seq_collection.GetSize();
  num_bases = collection_size.first;
  num_reads = collection_size.second;

  seq_pkg_.ReserveSequences(num_reads);
  seq_pkg_.ReserveBases(num_bases);

  seq_collection.Read(&seq_pkg_, is_reverse);

  seq_pkg_.BuildIndex();
  num_reads = seq_pkg_.seq_count();
  xinfo("{} reads, {} max read length\n", num_reads, seq_pkg_.max_length());

  words_per_substr_ =
      DivCeiling((opt_.k + 1) * kBitsPerEdgeChar, kBitsPerEdgeWord);
  words_per_edge_ = DivCeiling((opt_.k + 1) * kBitsPerEdgeChar + kBitsPerMul,
                               kBitsPerEdgeWord);
  xinfo("{} words per substring, {} words per edge\n", words_per_substr_,
        words_per_edge_);

  // --- malloc read first_in / last_out ---
  first_0_out_ =
      std::vector<AtomicWrapper<uint32_t>>(seq_pkg_.seq_count(), 0xFFFFFFFFU);
  last_0_in_ =
      std::vector<AtomicWrapper<uint32_t>>(seq_pkg_.seq_count(), 0xFFFFFFFFU);

  // --- initialize stat ---
  edge_counter_.SetNumThreads(opt_.n_threads);

  // --- initialize writer ---
  edge_writer_.SetFilePrefix(opt_.output_prefix);
  edge_writer_.SetNumThreads(opt_.n_threads);
  edge_writer_.SetKmerSize(opt_.k);
  edge_writer_.SetNumBuckets(kNumBuckets);
  edge_writer_.InitFiles();

  int64_t memory_for_data = seq_pkg_.size_in_byte() +
                            +seq_pkg_.seq_count() * sizeof(first_0_out_[0]) *
                                2  // first_in0 & last_out0
                            + edge_counter_.size_in_byte();  // edge_counting

  return {
      num_reads,
      memory_for_data,
      words_per_substr_ + 2,
      2,
  };
}

void KmerCounter::Lv0CalcBucketSize(int64_t seq_from, int64_t seq_to,
                                    std::array<int64_t, kNumBuckets> *out) {
  std::array<int64_t, kNumBuckets> &bucket_sizes = *out;
  std::fill(bucket_sizes.begin(), bucket_sizes.end(), 0);
  GenericKmer edge, rev_edge;  // (k+1)-mer and its rc

  for (int64_t read_id = seq_from; read_id < seq_to; ++read_id) {
    auto seq_view = seq_pkg_.GetSeqView(read_id);
    auto read_length = seq_view.length();

    if (read_length < opt_.k + 1) {
      continue;
    }

    auto start_ptr_and_offset = seq_view.raw_address();
    edge.InitFromPtr(start_ptr_and_offset.first, start_ptr_and_offset.second,
                     opt_.k + 1);
    rev_edge = edge;
    rev_edge.ReverseComplement(opt_.k + 1);

    unsigned last_char_offset = opt_.k;

    while (true) {
      if (rev_edge.cmp(edge, opt_.k + 1) < 0) {
        bucket_sizes[rev_edge.data()[0] >>
                     (kCharsPerEdgeWord - kBucketPrefixLength) *
                         kBitsPerEdgeChar]++;
      } else {
        bucket_sizes[edge.data()[0] >>
                     (kCharsPerEdgeWord - kBucketPrefixLength) *
                         kBitsPerEdgeChar]++;
      }

      if (++last_char_offset >= read_length) {
        break;
      } else {
        int c = seq_view.base_at(last_char_offset);
        edge.ShiftAppend(c, opt_.k + 1);
        rev_edge.ShiftPreappend(3 - c, opt_.k + 1);
      }
    }
  }
}

void KmerCounter::Lv1FillOffsets(OffsetFiller &filler, int64_t seq_from,
                                 int64_t seq_to) {
  GenericKmer edge, rev_edge;  // (k+1)-mer and its rc
  unsigned key;
  for (int64_t read_id = seq_from; read_id < seq_to; ++read_id) {
    auto seq_view = seq_pkg_.GetSeqView(read_id);
    auto read_length = seq_view.length();

    if (read_length < opt_.k + 1) {
      continue;
    }

    auto ptr_and_offset = seq_view.raw_address();
    edge.InitFromPtr(ptr_and_offset.first, ptr_and_offset.second, opt_.k + 1);
    rev_edge = edge;
    rev_edge.ReverseComplement(opt_.k + 1);

    // shift the key char by char
    unsigned last_char_offset = opt_.k;

    while (true) {
      if (rev_edge.cmp(edge, opt_.k + 1) < 0) {
        key = rev_edge.data()[0] >>
              (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
        if (filler.IsHandling(key)) {
          filler.WriteNextOffset(
              key,
              EncodeOffset(read_id, last_char_offset - opt_.k, 1, seq_pkg_));
        }
      } else {
        key = edge.data()[0] >>
              (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
        if (filler.IsHandling(key)) {
          filler.WriteNextOffset(
              key,
              EncodeOffset(read_id, last_char_offset - opt_.k, 0, seq_pkg_));
        }
      }

      if (++last_char_offset >= read_length) {
        break;
      } else {
        int c = seq_view.base_at(last_char_offset);
        edge.ShiftAppend(c, opt_.k + 1);
        rev_edge.ShiftPreappend(3 - c, opt_.k + 1);
      }
    }
  }
}

void KmerCounter::Lv2ExtractSubString(OffsetFetcher &fetcher,
                                      SubstrPtr substr_ptr) {
  while (fetcher.HasNext()) {
    int64_t full_offset = fetcher.Next();
    auto seq_view = seq_pkg_.GetSeqViewByOffset(full_offset >> 1u);
    unsigned strand = full_offset & 1;
    unsigned offset = (full_offset >> 1) - seq_view.full_offset_in_pkg();
    unsigned num_chars_to_copy = opt_.k + 1;

    unsigned read_length = seq_view.length();
    auto ptr_and_offset = seq_view.raw_address();
    unsigned start_offset = ptr_and_offset.second;
    unsigned words_this_seq = DivCeiling(start_offset + read_length, 16);
    const uint32_t *read_p = ptr_and_offset.first;

    unsigned char prev, next;

    if (offset > 0) {
      prev = seq_view.base_at(offset - 1);
    } else {
      prev = kSentinelValue;
    }

    if (offset + opt_.k + 1 < read_length) {
      next = seq_view.base_at(offset + opt_.k + 1);
    } else {
      next = kSentinelValue;
    }
    uint64_t read_info;
    if (strand == 0) {
      CopySubstring(substr_ptr, read_p, offset + start_offset,
                    num_chars_to_copy, 1, words_this_seq, words_per_substr_);
      read_info = (full_offset << 6) | (prev << 3) | next;
    } else {
      CopySubstringRC(substr_ptr, read_p, offset + start_offset,
                      num_chars_to_copy, 1, words_this_seq, words_per_substr_);
      read_info =
          (full_offset << 6) |
          ((next == kSentinelValue ? kSentinelValue : (3 - next)) << 3) |
          (prev == kSentinelValue ? kSentinelValue : (3 - prev));
    }
    DecomposeUint64(substr_ptr + words_per_substr_, read_info);
    substr_ptr += words_per_substr_ + 2;
  }
}

void KmerCounter::Lv2Postprocess(int64_t start_index, int64_t end_index,
                                 int thread_id, uint32_t *substr_ptr) {
  uint32_t packed_edge[32];
  int64_t count_prev[5], count_next[5];

  int64_t from_;
  int64_t to_;

  EdgeWriter::Snapshot snapshot;

  for (int64_t i = start_index; i < end_index; i = to_) {
    from_ = i;
    to_ = i + 1;
    uint32_t *first_item = substr_ptr + i * (words_per_substr_ + 2);

    while (to_ < end_index) {
      if (IsDifferentEdges(first_item,
                           substr_ptr + to_ * (words_per_substr_ + 2),
                           words_per_substr_, 1)) {
        break;
      }

      ++to_;
    }

    int64_t count = to_ - from_;

    // update read's first and last

    memset(count_prev, 0, sizeof(count_prev[0]) * 4);
    memset(count_next, 0, sizeof(count_next[0]) * 4);
    bool has_in = false;
    bool has_out = false;

    for (int64_t j = from_; j < to_; ++j) {
      auto read_info_ptr =
          substr_ptr + j * (words_per_substr_ + 2) + words_per_substr_;
      uint64_t read_info = ComposeUint64(read_info_ptr);
      int prev_and_next = read_info & ((1 << 6) - 1);
      count_prev[prev_and_next >> 3]++;
      count_next[prev_and_next & 7]++;
    }

    for (int j = 0; j < 4; ++j) {
      if (count_prev[j] >= opt_.solid_threshold) {
        has_in = true;
      }

      if (count_next[j] >= opt_.solid_threshold) {
        has_out = true;
      }
    }

    if (!has_in && count >= opt_.solid_threshold) {
      for (int64_t j = from_; j < to_; ++j) {
        auto *read_info_ptr =
            substr_ptr + j * (words_per_substr_ + 2) + words_per_substr_;
        int64_t read_info = ComposeUint64(read_info_ptr) >> 6;
        auto seq_view = seq_pkg_.GetSeqViewByOffset(read_info >> 1);
        int strand = read_info & 1;
        uint32_t offset = (read_info >> 1) - seq_view.full_offset_in_pkg();

        if (strand == 0) {
          // update last
          uint32_t old_value = last_0_in_[seq_view.id()].v.load(
              std::memory_order::memory_order_acquire);
          while ((old_value == kSentinelOffset || old_value < offset) &&
                 !last_0_in_[seq_view.id()].v.compare_exchange_weak(
                     old_value, offset, std::memory_order::memory_order_release,
                     std::memory_order::memory_order_relaxed)) {
          }
        } else {
          // update first
          offset++;
          uint32_t old_value = first_0_out_[seq_view.id()].v.load(
              std::memory_order::memory_order_acquire);
          while (old_value > offset &&
                 !first_0_out_[seq_view.id()].v.compare_exchange_weak(
                     old_value, offset, std::memory_order::memory_order_release,
                     std::memory_order::memory_order_relaxed)) {
          }
        }
      }
    }

    if (!has_out && count >= opt_.solid_threshold) {
      for (int64_t j = from_; j < to_; ++j) {
        auto *read_info_ptr =
            substr_ptr + j * (words_per_substr_ + 2) + words_per_substr_;
        int64_t read_info = ComposeUint64(read_info_ptr) >> 6;
        auto seq_view = seq_pkg_.GetSeqViewByOffset(read_info >> 1);
        int strand = read_info & 1;
        uint32_t offset = (read_info >> 1) - seq_view.full_offset_in_pkg();

        if (strand == 0) {
          // update first
          offset++;
          uint32_t old_value = first_0_out_[seq_view.id()].v.load(
              std::memory_order::memory_order_acquire);
          while (old_value > offset &&
                 !first_0_out_[seq_view.id()].v.compare_exchange_weak(
                     old_value, offset, std::memory_order::memory_order_release,
                     std::memory_order::memory_order_relaxed)) {
          }
        } else {
          // update last
          uint32_t old_value = last_0_in_[seq_view.id()].v.load(
              std::memory_order::memory_order_acquire);
          while ((old_value == kSentinelOffset || old_value < offset) &&
                 !last_0_in_[seq_view.id()].v.compare_exchange_weak(
                     old_value, offset, std::memory_order::memory_order_release,
                     std::memory_order::memory_order_relaxed)) {
          }
        }
      }
    }
    edge_counter_.Add(count, thread_id);

    if (count >= opt_.solid_threshold) {
      PackEdge(packed_edge, first_item, count);
      edge_writer_.Write(packed_edge,
                         packed_edge[0] >> (32 - 2 * kBucketPrefixLength),
                         thread_id, &snapshot);
    }
  }

  edge_writer_.SaveSnapshot(snapshot);
}

void KmerCounter::Lv0Postprocess() {
  // --- output reads for mercy ---
  int64_t num_candidate_reads = 0;
  int64_t num_has_tips = 0;
  std::ofstream candidate_file(opt_.output_prefix + ".cand",
                               std::ofstream::binary | std::ofstream::out);

  for (size_t i = 0; i < seq_pkg_.seq_count(); ++i) {
    auto first =
        first_0_out_[i].v.load(std::memory_order::memory_order_relaxed);
    auto last = last_0_in_[i].v.load(std::memory_order::memory_order_relaxed);

    if (first != kSentinelOffset && last != kSentinelOffset) {
      ++num_has_tips;

      if (last > first) {
        ++num_candidate_reads;
        seq_pkg_.WriteSequences(candidate_file, i, i);
      }
    }
  }

  xinfo("Total number of candidate reads: {} ({})\n", num_candidate_reads,
        num_has_tips);
  xinfo("Total number of solid edges: {}\n",
        edge_counter_.GetNumSolidEdges(opt_.solid_threshold));
  std::ofstream counting_file(opt_.output_prefix + ".counting");
  edge_counter_.DumpStat(counting_file);

  // --- cleaning ---
  edge_writer_.Finalize();
}
