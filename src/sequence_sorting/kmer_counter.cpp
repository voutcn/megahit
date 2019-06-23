#include "kmer_counter.h"

#include <omp.h>
#include <string.h>
#include <zlib.h>
#include <algorithm>
#include <mutex>

#include "sequence/kmer.h"
#include "sequence/packed_reads.h"
#include "sequence/lib_io.h"
#include "sequence/readers/kseq.h"
#include "utils/safe_open.h"
#include "utils/utils.h"

#include "sorting.h"

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
  for (int i = 0; i < words_per_edge && i < words_per_substring; ++i) {
    dest[i] = *(item + i);
  }

  int chars_in_last_word = (kmer_k + 1) % kCharsPerEdgeWord;
  int which_word = (kmer_k + 1) / kCharsPerEdgeWord;

  if (chars_in_last_word > 0) {
    dest[which_word] >>= (kCharsPerEdgeWord - chars_in_last_word) * 2;
    dest[which_word] <<= (kCharsPerEdgeWord - chars_in_last_word) * 2;
  } else {
    dest[which_word] = 0;
  }

  while (++which_word < words_per_edge) {
    dest[which_word] = 0;
  }

  dest[words_per_edge - 1] |= std::min(int64_t(kMaxMul), counting);
}

// function pass to BaseSequenceSortingEngine

int64_t KmerCounter::encode_lv1_diff_base_func_(int64_t read_id) {
  return EncodeOffset(read_id, 0, 0, package);
}

void KmerCounter::prepare_func_() {
  bool is_reverse = true;

  int64_t num_bases, num_reads;
  GetBinaryLibSize(read_lib_file, num_bases, num_reads);

  if (assist_seq_file != "") {
    FILE *assist_seq_info = xfopen((assist_seq_file + ".info").c_str(), "r");
    long long num_ass_bases, num_ass_seq;
    if (fscanf(assist_seq_info, "%lld%lld", &num_ass_seq, &num_ass_bases) != 2) {
      xfatal("Invalid format\n");
    }
    fclose(assist_seq_info);

    num_bases += num_ass_bases;
    num_reads += num_ass_seq;
  }

  package.ReserveSequences(num_reads);
  package.ReserveBases(num_bases);

  std::vector<lib_info_t> lib_info;
  ReadBinaryLibs(read_lib_file, package, lib_info, is_reverse);

  if (assist_seq_file != "") {
    FastxReader reader({assist_seq_file});
    reader.ReadAll(&package, is_reverse);
  }

  package.BuildIndex();
  num_reads = package.Size();

  xinfo("%ld reads, %d max read length\n", num_reads, package.MaxSequenceLength());

  int64_t mem_packed_reads = package.SizeInByte();

  int64_t mem_low_bound = mem_packed_reads +
      num_reads * sizeof(unsigned char) * 2  // first_in0 & last_out0
      + kNumBuckets * sizeof(int64_t) * (num_cpu_threads * 3 + 1) +
      (kMaxMul + 1) * (num_cpu_threads + 1) * sizeof(int64_t);
  mem_low_bound *= 1.05;

  if (mem_low_bound > host_mem) {
    xfatal("%lld bytes is not enough for BaseSequenceSortingEngine sorting, please set -m parameter to at least %lld\n", host_mem,
           mem_low_bound);
  }

  // set sequence_sorting param
  SetNumCpuThreads(num_cpu_threads);;
  SetNumItems(num_reads);
}

void KmerCounter::lv0_calc_bucket_size_func_(ReadPartition *_data) {
  auto &rp = *(_data);
  std::array<int64_t, kNumBuckets> &bucket_sizes = rp.rp_bucket_sizes;
  std::fill(bucket_sizes.begin(), bucket_sizes.end(), 0);
  GenericKmer edge, rev_edge;  // (k+1)-mer and its rc

  for (int64_t read_id = rp.rp_start_id; read_id < rp.rp_end_id; ++read_id) {
    int read_length = package.SequenceLength(read_id);

    if (read_length < kmer_k + 1) {
      continue;
    }

    auto start_ptr_and_offset = package.WordPtrAndOffset(read_id);
    edge.InitFromPtr(start_ptr_and_offset.first, start_ptr_and_offset.second, kmer_k + 1);
    rev_edge = edge;
    rev_edge.ReverseComplement(kmer_k + 1);

    int last_char_offset = kmer_k;

    while (true) {
      if (rev_edge.cmp(edge, kmer_k + 1) < 0) {
        bucket_sizes[rev_edge.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
      } else {
        bucket_sizes[edge.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
      }

      if (++last_char_offset >= read_length) {
        break;
      } else {
        int c = package.GetBase(read_id, last_char_offset);
        edge.ShiftAppend(c, kmer_k + 1);
        rev_edge.ShiftPreappend(3 - c, kmer_k + 1);
      }
    }
  }
}

void KmerCounter::init_global_and_set_cx1_func_() {
  // --- calculate lv2 memory ---
  int64_t max_bucket_size =
      *std::max_element(GetBucketSizes().begin(), GetBucketSizes().end());
  int64_t tot_bucket_size = 0;
  int num_non_empty = 0;

  for (unsigned i = 0; i < kNumBuckets; ++i) {
    tot_bucket_size += GetBucketSizes()[i];

    if (GetBucketSizes()[i] > 0) {
      num_non_empty++;
    }
  }

  words_per_substring = DivCeiling((kmer_k + 1) * kBitsPerEdgeChar, kBitsPerEdgeWord);
  words_per_edge = DivCeiling((kmer_k + 1) * kBitsPerEdgeChar + kBitsPerMul, kBitsPerEdgeWord);
  xinfo("%d words per substring, %d words per edge\n", words_per_substring, words_per_edge);

  num_non_empty = std::max(1, num_non_empty);

  int64_t max_sorting_items =
      std::max(3 * tot_bucket_size / num_non_empty * num_cpu_threads, max_bucket_size);

  // lv2 bytes: substring, readinfo
  int64_t lv2_bytes_per_item = words_per_substring * sizeof(uint32_t) + sizeof(int64_t);

  int64_t mem_remained = host_mem - package.SizeInByte() -
      num_cpu_threads * 65536 * sizeof(uint64_t)  // radix sort buckets
      - package.Size() * sizeof(unsigned char) * 2     // first_in0 & last_out0
      - kNumBuckets * sizeof(int64_t) * (num_cpu_threads * 3 + 1) -
      (kMaxMul + 1) * (num_cpu_threads + 1) * sizeof(int64_t);
  int64_t min_lv1_items = tot_bucket_size / (kMaxLv1ScanTime - 0.5);
  int64_t max_lv1_items = 0;

  if (mem_flag == 1) {
    // auto set memory
    max_lv1_items = int64_t(tot_bucket_size / (kDefaultLv1ScanTime - 0.5));
    max_lv1_items = std::max(max_lv1_items, max_bucket_size);
    int64_t mem_needed =
        max_lv1_items * kLv1BytePerItem + max_sorting_items * lv2_bytes_per_item;
    if (mem_needed > mem_remained) {
      adjust_mem(mem_remained, lv2_bytes_per_item, min_lv1_items, max_bucket_size,
                 max_sorting_items, max_lv1_items, max_sorting_items);
    }

  } else if (mem_flag == 0) {
    // min memory
    max_lv1_items = int64_t(tot_bucket_size / (kMaxLv1ScanTime - 0.5));
    max_lv1_items = std::max(max_lv1_items, max_bucket_size);
    int64_t mem_needed =
        max_lv1_items * kLv1BytePerItem + max_sorting_items * lv2_bytes_per_item;

    if (mem_needed > mem_remained) {
      adjust_mem(mem_remained, lv2_bytes_per_item, min_lv1_items, max_bucket_size,
                 max_sorting_items, max_lv1_items, max_sorting_items);
    } else {
      adjust_mem(mem_needed, lv2_bytes_per_item, min_lv1_items, max_bucket_size,
                 max_sorting_items, max_lv1_items, max_sorting_items);
    }

  } else {
    // use all
    adjust_mem(mem_remained, lv2_bytes_per_item, min_lv1_items, max_bucket_size,
               max_sorting_items, max_lv1_items, max_sorting_items);
  }

  if (max_lv1_items < min_lv1_items) {
    xfatal("No enough memory to process.");
  }

  SetMaxLv1Lv2Items(max_lv1_items, max_sorting_items, lv2_bytes_per_item / sizeof(uint32_t), 2);
  xinfo("Memory for reads: %lld\n", package.SizeInByte());
  xinfo("max # lv.1 items = %lld\n", max_lv1_items);

  // --- malloc read first_in / last_out ---
  first_0_out = std::vector<AtomicWrapper<uint32_t>>(package.Size(), 0xFFFFFFFFU);
  last_0_in = std::vector<AtomicWrapper<uint32_t>>(package.Size(), 0xFFFFFFFFU);

  // --- initialize stat ---
  thread_edge_counting.resize(num_cpu_threads);
  for (auto &c : thread_edge_counting) {
    c.resize(kMaxMul + 1);
    std::fill(c.begin(), c.end(), 0);
  }

  // --- initialize writer ---
  edge_writer.SetFilePrefix(output_prefix);
  edge_writer.SetNumThreads(num_cpu_threads);
  edge_writer.SetKmerSize(kmer_k);
  edge_writer.SetNumBuckets(kNumBuckets);
  edge_writer.InitFiles();
}

void KmerCounter::lv1_fill_offset_func_(ReadPartition *_data) {
  auto &rp = *_data;
  std::array<int64_t, kNumBuckets> prev_full_offsets;

  for (auto b = GetLv1StartBucket(); b < GetLv1EndBucket(); ++b)
    prev_full_offsets[b] = rp.rp_lv1_differential_base;

  // this loop is VERY similar to that in PreprocessScanToFillBucketSizesThread
  GenericKmer edge, rev_edge;  // (k+1)-mer and its rc
  int key;

  for (int64_t read_id = rp.rp_start_id; read_id < rp.rp_end_id; ++read_id) {
    int read_length = package.SequenceLength(read_id);

    if (read_length < kmer_k + 1) {
      continue;
    }

    auto ptr_and_offset = package.WordPtrAndOffset(read_id);
    edge.InitFromPtr(ptr_and_offset.first, ptr_and_offset.second, kmer_k + 1);
    rev_edge = edge;
    rev_edge.ReverseComplement(kmer_k + 1);

    // ===== this is a macro to save some copy&paste ================
#define CHECK_AND_SAVE_OFFSET(offset, strand)                                            \
  do {                                                                                   \
    if (HandlingBucket(key)) {                                                           \
      int key_ = GetBucketRank(key);                                                     \
      int64_t full_offset = EncodeOffset(read_id, offset, strand, package);              \
      int64_t differential = full_offset - prev_full_offsets[key_];                      \
      int64_t index = rp.rp_bucket_offsets[key_]++;                                      \
      WriteOffset(index, differential, full_offset);                                     \
      assert(rp.rp_bucket_offsets[key_] <= GetLv1NumItems());                            \
      prev_full_offsets[key_] = full_offset;                                             \
    }                                                                                    \
  } while (0)
    // ^^^^^ why is the macro surrounded by a do-while? please ask Google
    // =========== end macro ==========================

    // shift the key char by char
    int last_char_offset = kmer_k;

    while (true) {
      if (rev_edge.cmp(edge, kmer_k + 1) < 0) {
        key = rev_edge.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
        CHECK_AND_SAVE_OFFSET(last_char_offset - kmer_k, 1);
      } else {
        key = edge.data()[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
        CHECK_AND_SAVE_OFFSET(last_char_offset - kmer_k, 0);
      }

      if (++last_char_offset >= read_length) {
        break;
      } else {
        int c = package.GetBase(read_id, last_char_offset);
        edge.ShiftAppend(c, kmer_k + 1);
        rev_edge.ShiftPreappend(3 - c, kmer_k + 1);
      }
    }
  }

#undef CHECK_AND_SAVE_OFFSET
}

void KmerCounter::lv2_extract_substr_(unsigned start_bucket, unsigned end_bucket, uint32_t *substrings_p) {
  auto lv1_p = GetLv1Iterator(start_bucket);

  for (auto b = start_bucket; b < end_bucket; ++b) {
    for (int t = 0; t < num_cpu_threads; ++t) {
      int64_t full_offset = GetReadPartition(t).rp_lv1_differential_base;
      int64_t num = GetReadPartition(t).rp_bucket_sizes[b];

      for (int64_t i = 0; i < num; ++i) {
        if (*lv1_p >= 0) {
          full_offset += *(lv1_p++);
        } else {
          full_offset = GetSpecialOffset(-1 - *(lv1_p++));
        }

        int64_t read_id = package.GetSeqID(full_offset >> 1);
        int strand = full_offset & 1;
        int offset = (full_offset >> 1) - package.StartPos(read_id);
        int num_chars_to_copy = kmer_k + 1;

        int read_length = package.SequenceLength(read_id);
        auto ptr_and_offset = package.WordPtrAndOffset(read_id);
        int start_offset = ptr_and_offset.second;
        int words_this_seq = DivCeiling(start_offset + read_length, 16);
        const uint32_t *read_p = ptr_and_offset.first;

        unsigned char prev, next;

        if (offset > 0) {
          prev = package.GetBase(read_id, offset - 1);
        } else {
          prev = kSentinelValue;
        }

        if (offset + kmer_k + 1 < read_length) {
          next = package.GetBase(read_id, offset + kmer_k + 1);
        } else {
          next = kSentinelValue;
        }
        uint64_t *read_info_p = reinterpret_cast<uint64_t *>(substrings_p + words_per_substring);
        if (strand == 0) {
          CopySubstring(substrings_p, read_p, offset + start_offset, num_chars_to_copy, 1, words_this_seq,
                        words_per_substring);
          *read_info_p = (full_offset << 6) | (prev << 3) | next;
        } else {
          CopySubstringRC(substrings_p, read_p, offset + start_offset, num_chars_to_copy, 1, words_this_seq,
                          words_per_substring);
          *read_info_p = (full_offset << 6) | ((next == kSentinelValue ? kSentinelValue : (3 - next)) << 3) |
              (prev == kSentinelValue ? kSentinelValue : (3 - prev));
        }

        substrings_p += words_per_substring + 2;
      }
    }
  }
}

void KmerCounter::output_(int64_t start_index, int64_t end_index, int thread_id, uint32_t *substrings) {
  uint32_t packed_edge[32];
  int64_t count_prev[5], count_next[5];
  auto &thread_edge_count = thread_edge_counting[thread_id];

  int64_t from_;
  int64_t to_;

  EdgeWriter::Snapshot snapshot;

  for (int64_t i = start_index; i < end_index; i = to_) {
    from_ = i;
    to_ = i + 1;
    uint32_t *first_item = substrings + i * (words_per_substring + 2);

    while (to_ < end_index) {
      if (IsDifferentEdges(first_item, substrings + to_ * (words_per_substring + 2),
                           words_per_substring, 1)) {
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
      auto *read_info = reinterpret_cast<uint64_t *>(substrings + j * (words_per_substring + 2) +
          words_per_substring);
      int prev_and_next = *read_info & ((1 << 6) - 1);
      count_prev[prev_and_next >> 3]++;
      count_next[prev_and_next & 7]++;
    }

    for (int j = 0; j < 4; ++j) {
      if (count_prev[j] >= kmer_freq_threshold) {
        has_in = true;
      }

      if (count_next[j] >= kmer_freq_threshold) {
        has_out = true;
      }
    }

    if (!has_in && count >= kmer_freq_threshold) {
      for (int64_t j = from_; j < to_; ++j) {
        auto *read_info_ptr = reinterpret_cast<uint64_t *>(substrings + j * (words_per_substring + 2) +
            words_per_substring);
        int64_t read_info = *read_info_ptr >> 6;
        int64_t read_id = package.GetSeqID(read_info >> 1);
        int strand = read_info & 1;
        uint32_t offset = (read_info >> 1) - package.StartPos(read_id);

        if (strand == 0) {
          // update last
          uint32_t old_value = last_0_in[read_id].v.load(std::memory_order::memory_order_acquire);
          while ((old_value == kSentinelOffset || old_value < offset) &&
              !last_0_in[read_id].v.compare_exchange_weak(old_value, offset,
                                                          std::memory_order::memory_order_release,
                                                          std::memory_order::memory_order_relaxed)) {
          }
        } else {
          // update first
          offset++;
          uint32_t old_value = first_0_out[read_id].v.load(std::memory_order::memory_order_acquire);
          while (old_value > offset && !first_0_out[read_id].v.compare_exchange_weak(
              old_value, offset, std::memory_order::memory_order_release,
              std::memory_order::memory_order_relaxed)) {
          }
        }
      }
    }

    if (!has_out && count >= kmer_freq_threshold) {
      for (int64_t j = from_; j < to_; ++j) {
        auto *read_info_ptr = reinterpret_cast<uint64_t *>(substrings + j * (words_per_substring + 2) +
            words_per_substring);
        int64_t read_info = *read_info_ptr >> 6;
        int64_t read_id = package.GetSeqID(read_info >> 1);
        int strand = read_info & 1;
        uint32_t offset = (read_info >> 1) - package.StartPos(read_id);

        if (strand == 0) {
          // update first
          offset++;
          uint32_t old_value = first_0_out[read_id].v.load(std::memory_order::memory_order_acquire);
          while (old_value > offset && !first_0_out[read_id].v.compare_exchange_weak(
              old_value, offset, std::memory_order::memory_order_release,
              std::memory_order::memory_order_relaxed)) {
          }
        } else {
          // update last
          uint32_t old_value = last_0_in[read_id].v.load(std::memory_order::memory_order_acquire);
          while ((old_value == kSentinelOffset || old_value < offset) &&
              !last_0_in[read_id].v.compare_exchange_weak(old_value, offset,
                                                          std::memory_order::memory_order_release,
                                                          std::memory_order::memory_order_relaxed)) {
          }
        }
      }
    }

    ++thread_edge_count[std::min(count, int64_t(kMaxMul))];

    if (count >= kmer_freq_threshold) {
      PackEdge(packed_edge, first_item, count);
      edge_writer.Write(packed_edge, packed_edge[0] >> (32 - 2 * kBucketPrefixLength), thread_id, &snapshot);
    }
  }

  edge_writer.SaveSnapshot(snapshot);
}

void KmerCounter::post_proc_func_() {
  std::vector<int64_t> edge_counting(kMaxMul + 1, 0);
  for (int t = 0; t < num_cpu_threads; ++t) {
    for (int i = 1; i <= kMaxMul; ++i) {
      edge_counting[i] += thread_edge_counting[t][i];
    }
  }

  // --- output reads for mercy ---
  int64_t num_candidate_reads = 0;
  int64_t num_has_tips = 0;
  FILE *candidate_file = xfopen((output_prefix + ".cand").c_str(), "wb");

  for (size_t i = 0; i < package.Size(); ++i) {
    auto first = first_0_out[i].v.load(std::memory_order::memory_order_relaxed);
    auto last = last_0_in[i].v.load(std::memory_order::memory_order_relaxed);

    if (first != kSentinelOffset && last != kSentinelOffset) {
      ++num_has_tips;

      if (last > first) {
        ++num_candidate_reads;
        WriteBinarySequences(package, candidate_file, i, i);
      }
    }
  }

  fclose(candidate_file);

  xinfo("Total number of candidate reads: %lld(%lld)\n", num_candidate_reads, num_has_tips);

  // --- stat ---
  int64_t num_solid_edges = 0;

  for (int i = kmer_freq_threshold; i <= kMaxMul; ++i) {
    num_solid_edges += edge_counting[i];
  }
  xinfo("Total number of solid edges: %llu\n", num_solid_edges);

  FILE *counting_file = xfopen((output_prefix + ".counting").c_str(), "w");

  for (int64_t i = 1, acc = 0; i <= kMaxMul; ++i) {
    acc += edge_counting[i];
    fprintf(counting_file, "%lld %lld\n", (long long) i, (long long) acc);
  }

  fclose(counting_file);

  // --- cleaning ---
  edge_writer.Finalize();
}