#include "kmer_counter.h"

#include <omp.h>
#include <algorithm>
#include <mutex>

#include "sequence/kmer.h"
#include "sequence/packed_reads.h"
#include "sequence/lib_io.h"
#include "sequence/readers/kseq.h"
#include "utils/safe_open.h"
#include "utils/utils.h"

/**
 * @brief encode read_id and its offset in one int64_t
 */
inline int64_t EncodeOffset(int64_t read_id, int offset, int strand, SeqPackage &p) {
  return ((p.StartPos(read_id) + offset) << 1) | strand;
}

inline bool IsDifferentEdges(uint32_t *item1, uint32_t *item2, int num_words, int64_t spacing) {
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

KmerCounter::Meta KmerCounter::Initialize() {
  bool is_reverse = true;

  int64_t num_bases, num_reads;
  GetBinaryLibSize(opt_.read_lib_file, num_bases, num_reads);

  seq_pkg_.ReserveSequences(num_reads);
  seq_pkg_.ReserveBases(num_bases);

  std::vector<lib_info_t> lib_info;
  ReadBinaryLibs(opt_.read_lib_file, seq_pkg_, lib_info, is_reverse);

  seq_pkg_.BuildIndex();
  num_reads = seq_pkg_.SeqCount();
  xinfo("{} reads, {} max read length\n", num_reads, seq_pkg_.MaxSequenceLength());

  words_per_substr_ = DivCeiling((opt_.k + 1) * kBitsPerEdgeChar, kBitsPerEdgeWord);
  words_per_edge_ = DivCeiling((opt_.k + 1) * kBitsPerEdgeChar + kBitsPerMul, kBitsPerEdgeWord);
  xinfo("{} words per substring, {} words per edge\n", words_per_substr_, words_per_edge_);

  // --- malloc read first_in / last_out ---
  first_0_out_ = std::vector<AtomicWrapper<uint32_t>>(seq_pkg_.SeqCount(), 0xFFFFFFFFU);
  last_0_in_ = std::vector<AtomicWrapper<uint32_t>>(seq_pkg_.SeqCount(), 0xFFFFFFFFU);

  // --- initialize stat ---
  thread_edge_counting_.resize(opt_.n_threads);
  for (auto &c : thread_edge_counting_) {
    c.resize(kMaxMul + 1);
    std::fill(c.begin(), c.end(), 0);
  }

  // --- initialize writer ---
  edge_writer_.SetFilePrefix(opt_.output_prefix);
  edge_writer_.SetNumThreads(opt_.n_threads);
  edge_writer_.SetKmerSize(opt_.k);
  edge_writer_.SetNumBuckets(kNumBuckets);
  edge_writer_.InitFiles();

  int64_t memory_for_data = seq_pkg_.SizeInByte() +
      +seq_pkg_.SeqCount() * sizeof(first_0_out_[0]) * 2     // first_in0 & last_out0
      + (kMaxMul + 1) * (opt_.n_threads + 1) * sizeof(int64_t);  // edge_counting

  return {
      num_reads,
      memory_for_data,
      words_per_substr_ + 2,
      2,
  };
}

void KmerCounter::Lv0CalcBucketSize(int64_t seq_from, int64_t seq_to, std::array<int64_t, kNumBuckets> *out) {
  std::array<int64_t, kNumBuckets> &bucket_sizes = *out;
  std::fill(bucket_sizes.begin(), bucket_sizes.end(), 0);
  GenericKmer edge, rev_edge;  // (k+1)-mer and its rc

  for (int64_t read_id = seq_from; read_id < seq_to; ++read_id) {
    auto read_length = seq_pkg_.SequenceLength(read_id);

    if (read_length < opt_.k + 1) {
      continue;
    }

    auto start_ptr_and_offset = seq_pkg_.WordPtrAndOffset(read_id);
    edge.InitFromPtr(start_ptr_and_offset.first, start_ptr_and_offset.second, opt_.k + 1);
    rev_edge = edge;
    rev_edge.ReverseComplement(opt_.k + 1);

    unsigned last_char_offset = opt_.k;

    while (true) {
      if (rev_edge.cmp(edge, opt_.k + 1) < 0) {
        bucket_sizes[rev_edge.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
      } else {
        bucket_sizes[edge.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
      }

      if (++last_char_offset >= read_length) {
        break;
      } else {
        int c = seq_pkg_.GetBase(read_id, last_char_offset);
        edge.ShiftAppend(c, opt_.k + 1);
        rev_edge.ShiftPreappend(3 - c, opt_.k + 1);
      }
    }
  }
}

void KmerCounter::Lv1FillOffsets(OffsetFiller &filler, int64_t seq_from, int64_t seq_to) {
  GenericKmer edge, rev_edge;  // (k+1)-mer and its rc
  unsigned key;
  for (int64_t read_id = seq_from; read_id < seq_to; ++read_id) {
    auto read_length = seq_pkg_.SequenceLength(read_id);

    if (read_length < opt_.k + 1) {
      continue;
    }

    auto ptr_and_offset = seq_pkg_.WordPtrAndOffset(read_id);
    edge.InitFromPtr(ptr_and_offset.first, ptr_and_offset.second, opt_.k + 1);
    rev_edge = edge;
    rev_edge.ReverseComplement(opt_.k + 1);

    // shift the key char by char
    unsigned last_char_offset = opt_.k;

    while (true) {
      if (rev_edge.cmp(edge, opt_.k + 1) < 0) {
        key = rev_edge.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
        if (filler.IsHandling(key)) {
          filler.WriteNextOffset(key, EncodeOffset(read_id, last_char_offset - opt_.k, 1, seq_pkg_));
        }
      } else {
        key = edge.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
        if (filler.IsHandling(key)) {
          filler.WriteNextOffset(key, EncodeOffset(read_id, last_char_offset - opt_.k, 0, seq_pkg_));
        }
      }

      if (++last_char_offset >= read_length) {
        break;
      } else {
        int c = seq_pkg_.GetBase(read_id, last_char_offset);
        edge.ShiftAppend(c, opt_.k + 1);
        rev_edge.ShiftPreappend(3 - c, opt_.k + 1);
      }
    }
  }
}

void KmerCounter::Lv2ExtractSubString(unsigned start_bucket, unsigned end_bucket, uint32_t *substr_ptr) {
  auto offset_iterator = GetOffsetFetcher(start_bucket, end_bucket);

  while (offset_iterator.HasNext()) {
    int64_t full_offset = offset_iterator.Next();
    int64_t read_id = seq_pkg_.GetSeqID(full_offset >> 1);
    unsigned strand = full_offset & 1;
    unsigned offset = (full_offset >> 1) - seq_pkg_.StartPos(read_id);
    unsigned num_chars_to_copy = opt_.k + 1;

    unsigned read_length = seq_pkg_.SequenceLength(read_id);
    auto ptr_and_offset = seq_pkg_.WordPtrAndOffset(read_id);
    unsigned start_offset = ptr_and_offset.second;
    unsigned words_this_seq = DivCeiling(start_offset + read_length, 16);
    const uint32_t *read_p = ptr_and_offset.first;

    unsigned char prev, next;

    if (offset > 0) {
      prev = seq_pkg_.GetBase(read_id, offset - 1);
    } else {
      prev = kSentinelValue;
    }

    if (offset + opt_.k + 1 < read_length) {
      next = seq_pkg_.GetBase(read_id, offset + opt_.k + 1);
    } else {
      next = kSentinelValue;
    }
    uint64_t read_info;
    if (strand == 0) {
      CopySubstring(substr_ptr, read_p, offset + start_offset, num_chars_to_copy, 1, words_this_seq,
                    words_per_substr_);
      read_info = (full_offset << 6) | (prev << 3) | next;
    } else {
      CopySubstringRC(substr_ptr, read_p, offset + start_offset, num_chars_to_copy, 1, words_this_seq,
                      words_per_substr_);
      read_info = (full_offset << 6) | ((next == kSentinelValue ? kSentinelValue : (3 - next)) << 3) |
          (prev == kSentinelValue ? kSentinelValue : (3 - prev));
    }
    DecomposeUint64(substr_ptr + words_per_substr_, read_info);
    substr_ptr += words_per_substr_ + 2;
  }
}

void KmerCounter::Lv2Postprocess(int64_t start_index, int64_t end_index, int thread_id, uint32_t *substrings) {
  uint32_t packed_edge[32];
  int64_t count_prev[5], count_next[5];
  auto &thread_edge_count = thread_edge_counting_[thread_id];

  int64_t from_;
  int64_t to_;

  EdgeWriter::Snapshot snapshot;

  for (int64_t i = start_index; i < end_index; i = to_) {
    from_ = i;
    to_ = i + 1;
    uint32_t *first_item = substrings + i * (words_per_substr_ + 2);

    while (to_ < end_index) {
      if (IsDifferentEdges(first_item, substrings + to_ * (words_per_substr_ + 2),
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
      auto read_info_ptr = substrings + j * (words_per_substr_ + 2) + words_per_substr_;
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
        auto *read_info_ptr = substrings + j * (words_per_substr_ + 2) + words_per_substr_;
        uint64_t read_info_context = ComposeUint64(read_info_ptr);
        int64_t read_info = read_info_context >> 6;
        int64_t read_id = seq_pkg_.GetSeqID(read_info >> 1);
        int strand = read_info & 1;
        uint32_t offset = (read_info >> 1) - seq_pkg_.StartPos(read_id);

        if (strand == 0) {
          // update last
          uint32_t old_value = last_0_in_[read_id].v.load(std::memory_order::memory_order_acquire);
          while ((old_value == kSentinelOffset || old_value < offset) &&
              !last_0_in_[read_id].v.compare_exchange_weak(old_value, offset,
                                                           std::memory_order::memory_order_release,
                                                           std::memory_order::memory_order_relaxed)) {
          }
        } else {
          // update first
          offset++;
          uint32_t old_value = first_0_out_[read_id].v.load(std::memory_order::memory_order_acquire);
          while (old_value > offset && !first_0_out_[read_id].v.compare_exchange_weak(
              old_value, offset, std::memory_order::memory_order_release,
              std::memory_order::memory_order_relaxed)) {
          }
        }
      }
    }

    if (!has_out && count >= opt_.solid_threshold) {
      for (int64_t j = from_; j < to_; ++j) {
        auto *read_info_ptr = substrings + j * (words_per_substr_ + 2) + words_per_substr_;
        uint64_t read_info_context = ComposeUint64(read_info_ptr);
        int64_t read_info = read_info_context >> 6;
        int64_t read_id = seq_pkg_.GetSeqID(read_info >> 1);
        int strand = read_info & 1;
        uint32_t offset = (read_info >> 1) - seq_pkg_.StartPos(read_id);

        if (strand == 0) {
          // update first
          offset++;
          uint32_t old_value = first_0_out_[read_id].v.load(std::memory_order::memory_order_acquire);
          while (old_value > offset && !first_0_out_[read_id].v.compare_exchange_weak(
              old_value, offset, std::memory_order::memory_order_release,
              std::memory_order::memory_order_relaxed)) {
          }
        } else {
          // update last
          uint32_t old_value = last_0_in_[read_id].v.load(std::memory_order::memory_order_acquire);
          while ((old_value == kSentinelOffset || old_value < offset) &&
              !last_0_in_[read_id].v.compare_exchange_weak(old_value, offset,
                                                           std::memory_order::memory_order_release,
                                                           std::memory_order::memory_order_relaxed)) {
          }
        }
      }
    }

    ++thread_edge_count[std::min(count, int64_t(kMaxMul))];

    if (count >= opt_.solid_threshold) {
      PackEdge(packed_edge, first_item, count);
      edge_writer_.Write(packed_edge, packed_edge[0] >> (32 - 2 * kBucketPrefixLength), thread_id, &snapshot);
    }
  }

  edge_writer_.SaveSnapshot(snapshot);
}

void KmerCounter::Lv0Postprocess() {
  std::vector<int64_t> edge_counting(kMaxMul + 1, 0);
  for (int t = 0; t < opt_.n_threads; ++t) {
    for (int i = 1; i <= kMaxMul; ++i) {
      edge_counting[i] += thread_edge_counting_[t][i];
    }
  }

  // --- output reads for mercy ---
  int64_t num_candidate_reads = 0;
  int64_t num_has_tips = 0;
  FILE *candidate_file = xfopen((opt_.output_prefix + ".cand").c_str(), "wb");

  for (size_t i = 0; i < seq_pkg_.SeqCount(); ++i) {
    auto first = first_0_out_[i].v.load(std::memory_order::memory_order_relaxed);
    auto last = last_0_in_[i].v.load(std::memory_order::memory_order_relaxed);

    if (first != kSentinelOffset && last != kSentinelOffset) {
      ++num_has_tips;

      if (last > first) {
        ++num_candidate_reads;
        WriteBinarySequences(seq_pkg_, candidate_file, i, i);
      }
    }
  }

  fclose(candidate_file);

  xinfo("Total number of candidate reads: {} ({})\n", num_candidate_reads, num_has_tips);

  // --- stat ---
  int64_t num_solid_edges = 0;

  for (int i = opt_.solid_threshold; i <= kMaxMul; ++i) {
    num_solid_edges += edge_counting[i];
  }
  xinfo("Total number of solid edges: {}\n", num_solid_edges);

  FILE *counting_file = xfopen((opt_.output_prefix + ".counting").c_str(), "w");

  for (int64_t i = 1, acc = 0; i <= kMaxMul; ++i) {
    acc += edge_counting[i];
    pfprintf(counting_file, "{} {}\n", i, acc);
  }

  fclose(counting_file);

  // --- cleaning ---
  edge_writer_.Finalize();
}
