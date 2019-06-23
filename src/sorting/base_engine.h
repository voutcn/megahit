#ifndef MEGAHIT_BASE_ENGINE_H
#define MEGAHIT_BASE_ENGINE_H

#include <cassert>
#include <cstdint>
#include <algorithm>
#include <vector>
#include <mutex>
#include <omp.h>

#include "utils/utils.h"
#include "kmsort_selector.h"

/**
 * The base class of sequence sorting engine
 */
class BaseSequenceSortingEngine {
 public:
  static const unsigned kBucketBase = 4;
  static const unsigned kBucketPrefixLength = 8;
  static const unsigned kNumBuckets = 65536;
  static const int64_t kDefaultLv1ScanTime = 8;
  static const int64_t kMaxLv1ScanTime = 128;

  static const unsigned kLv1BytePerItem = 4;  // 32-bit differential offset
  static const int64_t kDifferentialLimit = (1llu << 31u) - 1;

  //  Sequence divided into n_threads_ partitions
  //  Each partition divided into kNumBuckets blocks
  //  Each block divided into n_threads_ sub-blocks
  //
  //  <--------------p0--------------><--------------p1-------------->
  //  <------b0------><------b1------><------b0------><------b1------>
  //  <--t0--><--t1--><--t0--><--t1--><--t0--><--t1--><--t0--><--t1-->
  //
  //  thread 0 fills/fetches from t0
  //  thread 1 fills/fetches from t1

  struct SeqPartition {
    // local data for each read partition (i.e. a sub-range of input reads)
    int64_t from, to;  // start and end IDs of this read partition (end is exclusive)
    std::array<int64_t, kNumBuckets> rp_bucket_sizes;
    std::array<int64_t, kNumBuckets> rp_bucket_offsets;
    int64_t rp_lv1_differential_base;  // the initial offset globals.lv1_items
  };

  struct Meta {
    int64_t num_sequences;
    int64_t memory_for_data;
    int64_t words_per_lv2;
    int64_t aux_words_per_lv2;
  };

 public:
  BaseSequenceSortingEngine(int64_t mem, int mem_flag, int n_threads)
      : host_mem_(mem), mem_flag_(mem_flag), n_threads_(n_threads) {};
  virtual ~BaseSequenceSortingEngine() = default;

  /**
   * Data that will change at every Lv1 stage
   */
 private:
  unsigned lv1_start_bucket_{}, lv1_end_bucket_{};
  int64_t lv1_num_items_{};
  std::vector<bool> cur_lv1_buckets_;
  std::vector<int64_t> lv1_special_offsets_;
  std::vector<uint32_t> lv1_offsets_;
  std::mutex special_item_lock_;

  std::vector<SeqPartition> rp_;
  std::array<int64_t, kNumBuckets> bucket_sizes_{};
  std::array<int, kNumBuckets> ori_bucket_id_{};
  std::array<int, kNumBuckets> bucket_rank_{};

 protected:
  /**
   * Iterator to fetch full offset from given buckets
   */
  class OffsetIterator {
   public:
    OffsetIterator(BaseSequenceSortingEngine *engine, unsigned start_bucket, unsigned end_bucket) :
        engine_(engine), end_bucket_(end_bucket), cur_bucket_(start_bucket) {
      cur_thread_id_ = 0;
      cur_item_index_ = 0;
      cur_full_offset_ = engine_->rp_[cur_thread_id_].rp_lv1_differential_base;
      cur_n_items_ = engine_->rp_[cur_thread_id_].rp_bucket_sizes[cur_bucket_];
      diff_iter_ = engine_->lv1_offsets_.cbegin() + engine_->rp_[0].rp_bucket_offsets[cur_bucket_];
      TryRefill();
    }

    bool HasNext() const { return cur_bucket_ < end_bucket_; }

    int64_t Next() {
      int64_t diff = *diff_iter_;
      if (diff <= kDifferentialLimit) {
        cur_full_offset_ += diff;
      } else {
        cur_full_offset_ = engine_->GetSpecialOffset(diff);
      }

      auto ret = cur_full_offset_;
      ++diff_iter_;
      ++cur_item_index_;
      TryRefill();
      return ret;
    }

   private:
    void TryRefill() {
      while (cur_item_index_ >= cur_n_items_) {
        ++cur_thread_id_;
        if (cur_thread_id_ == engine_->n_threads_) {
          cur_thread_id_ = 0;
          ++cur_bucket_;
          if (cur_bucket_ >= end_bucket_) {
            return;
          }
        }
        cur_full_offset_ = engine_->rp_[cur_thread_id_].rp_lv1_differential_base;
        cur_item_index_ = 0;
        cur_n_items_ = engine_->rp_[cur_thread_id_].rp_bucket_sizes[cur_bucket_];
      }
    }
   private:
    BaseSequenceSortingEngine *engine_;
    unsigned end_bucket_;
    unsigned cur_bucket_;
    unsigned cur_thread_id_;
    int64_t cur_full_offset_;
    size_t cur_item_index_;
    size_t cur_n_items_;
    std::vector<uint32_t>::const_iterator diff_iter_;
  };

  unsigned GetLv1StartBucket() const { return lv1_start_bucket_; }
  unsigned GetLv1EndBucket() const { return lv1_end_bucket_; }
  int64_t GetLv1NumItems() const { return lv1_num_items_; }

  /**
   * whether the bucket is being handled in current stage
   */
  bool HandlingBucket(unsigned bucket) const {
    return cur_lv1_buckets_[bucket];
  }

  void WriteOffset(size_t index, int64_t diff, int64_t full_offset) {
    if (diff > kDifferentialLimit) {
      int64_t key = kDifferentialLimit + 1 + AddSpecialOffset(full_offset);
      if (key >= std::numeric_limits<uint32_t>::max()) {
        xfatal("Too many large difference items!");
      }
      lv1_offsets_[index] = key;
    } else {
      assert(diff >= 0);
      lv1_offsets_[index] = static_cast<int32_t>(diff);
    }
  }

  OffsetIterator GetOffsetIterator(unsigned start_bucket, unsigned end_bucket) {
    return OffsetIterator(this, start_bucket, end_bucket);
  }

  int GetBucketRank(unsigned bucket) const { return bucket_rank_[bucket]; }

 private:
  int64_t AddSpecialOffset(int64_t offset) {
    std::lock_guard<std::mutex> lk(special_item_lock_);
    int64_t ret = lv1_special_offsets_.size();
    lv1_special_offsets_.push_back(offset);
    return ret;
  }

  int64_t GetSpecialOffset(uint32_t key) const {
    assert(key > kDifferentialLimit);
    return lv1_special_offsets_[key - (kDifferentialLimit + 1)];
  }

 private:
  /**
   * Memory related data
   */
  int64_t host_mem_{};
  int mem_flag_{};
  unsigned n_threads_{};
  Meta meta_{}; // set by Initialize return
  std::function<void(uint32_t *, int64_t)> substr_sort_;  // set after Initialize

  /**
   * Interfaces used by `Run` and must be implemented in derived class
   */
 public:
  virtual Meta Initialize() = 0;

 protected:
  virtual int64_t Lv0EncodeDiffBase(int64_t) = 0;
  virtual void Lv0CalcBucketSize(SeqPartition *) = 0;
  virtual void Lv1FillOffsets(SeqPartition *) = 0;
  virtual void Lv2ExtractSubString(unsigned bucket_from, unsigned bucket_to, uint32_t *substr_ptr) = 0;
  virtual void Lv2Postprocess(int64_t start_index, int64_t end_index, int thread_id, uint32_t *substrings) = 0;
  virtual void Lv0Postprocess() = 0;

 private:
  /**
   * Adjust memory layout of the engine, i.e., how many lv1 and lv2 items
   */
  void AdjustMemory();

  void Lv0PrepareReadPartition() {
    rp_.resize(n_threads_);
    for (unsigned t = 0; t < n_threads_; ++t) {
      SeqPartition &rp = rp_[t];
      // distribute reads to partitions
      int64_t average = meta_.num_sequences / n_threads_;
      rp.from = t * average;
      rp.to = t < n_threads_ - 1 ? (t + 1) * average : meta_.num_sequences;
      rp.rp_lv1_differential_base = Lv0EncodeDiffBase(rp.from);
    }

    for (unsigned i = 0; i < kNumBuckets; ++i) {
      ori_bucket_id_[i] = i;
      bucket_rank_[i] = i;
    }
  }

  void Lv0ReorderBuckets() {
    std::vector<std::pair<int64_t, int> > tmp_v(kNumBuckets);

    for (unsigned i = 0; i < kNumBuckets; ++i) {
      tmp_v[i] = std::make_pair(bucket_sizes_[i], i);
    }

    std::sort(tmp_v.rbegin(), tmp_v.rend());

    for (unsigned i = 0; i < kNumBuckets; ++i) {
      bucket_sizes_[i] = tmp_v[i].first;
      ori_bucket_id_[i] = tmp_v[i].second;
      bucket_rank_[tmp_v[i].second] = i;
    }

    for (unsigned tid = 0; tid < n_threads_; ++tid) {
      auto old_rp_bucket_sizes = rp_[tid].rp_bucket_sizes;
      for (unsigned i = 0; i < kNumBuckets; ++i) {
        rp_[tid].rp_bucket_sizes[i] = old_rp_bucket_sizes[tmp_v[i].second];
      }
    }
  }

  unsigned Lv1FindEndBuckets(unsigned start_bucket) {
    unsigned end_bucket = start_bucket;
    unsigned used_threads = 0;
    int64_t num_lv2 = 0;
    lv1_num_items_ = 0;

    cur_lv1_buckets_.resize(kNumBuckets);
    std::fill(cur_lv1_buckets_.begin(), cur_lv1_buckets_.end(), false);

    while (end_bucket < kNumBuckets) {
      if (used_threads < n_threads_) {
        num_lv2 += bucket_sizes_[end_bucket];
        ++used_threads;
      }

      if (num_lv2 * meta_.words_per_lv2 + lv1_num_items_ + bucket_sizes_[end_bucket] >
          static_cast<int64_t>(lv1_offsets_.size())) {
        return end_bucket;
      }

      lv1_num_items_ += bucket_sizes_[end_bucket];
      cur_lv1_buckets_[ori_bucket_id_[end_bucket]] = true;
      ++end_bucket;
    }

    return kNumBuckets;
  }

  void Lv1ComputeOffset() {
    // compute "global" (thread 0) offsets first
    auto &offsets = rp_[0].rp_bucket_offsets;
    offsets[lv1_start_bucket_] = 0;

    for (unsigned b = lv1_start_bucket_ + 1; b < lv1_end_bucket_; ++b) {
      offsets[b] = offsets[b - 1] + bucket_sizes_[b - 1];  // accumulate
    }

    // then for each read partition
    for (unsigned t = 1; t < n_threads_; ++t) {
      auto &this_offsets = rp_[t].rp_bucket_offsets;
      auto &prev_offsets = rp_[t - 1].rp_bucket_offsets;
      auto &sizes = rp_[t - 1].rp_bucket_sizes;

      for (unsigned b = lv1_start_bucket_; b < lv1_end_bucket_; ++b) {
        this_offsets[b] = prev_offsets[b] + sizes[b];
      }
    }
  }

  void Lv0CalcBucketSizeLaunchMt() {
#pragma omp parallel for
    for (unsigned t = 0; t < n_threads_; ++t) {
      Lv0CalcBucketSize(&rp_[t]);
    }
    // sum up readpartitions bucketsizes to form global bucketsizes
    std::fill(bucket_sizes_.begin(), bucket_sizes_.end(), 0);

    // the array accesses in this loop are optimized by the compiler??
    for (unsigned t = 0; t < n_threads_; ++t) {
      for (unsigned b = 0; b < kNumBuckets; ++b) {
        bucket_sizes_[b] += rp_[t].rp_bucket_sizes[b];
      }
    }
  }

  inline void Lv1FillOffsetsLaunchMt() {
    lv1_special_offsets_.clear();
    Lv1ComputeOffset();

#pragma omp parallel for
    for (unsigned t = 0; t < n_threads_; ++t) {
      Lv1FillOffsets(&rp_[t]);
    }
    Lv1ComputeOffset();
  }

  struct Lv2ThreadStatus {
    std::vector<int64_t> thread_offset;
    std::vector<int64_t> rank;
    std::mutex mutex;
    int64_t acc{0};
    int seen{0};
  };

  void Lv2Sort(Lv2ThreadStatus *thread_status, unsigned b, int tid) {
    if (thread_status->thread_offset[tid] == -1) {
      std::lock_guard<std::mutex> lk(thread_status->mutex);
      thread_status->thread_offset[tid] = thread_status->acc;
      thread_status->acc += bucket_sizes_[b];
      thread_status->rank[tid] = thread_status->seen;
      thread_status->seen++;
    }

    if (bucket_sizes_[b] == 0) {
      return;
    }

    size_t offset = lv1_num_items_ + thread_status->thread_offset[tid] * meta_.words_per_lv2;
    auto substr_ptr = lv1_offsets_.data() + offset;
    Lv2ExtractSubString(b, b + 1, substr_ptr);
    substr_sort_(substr_ptr, bucket_sizes_[b]);
    Lv2Postprocess(0, bucket_sizes_[b], tid, substr_ptr);
  }

  void Lv1FetchAndSortLaunchMt() {
    Lv2ThreadStatus thread_status{};
    thread_status.thread_offset.resize(n_threads_, -1);
    thread_status.rank.resize(n_threads_, 0);
    omp_set_num_threads(n_threads_);
#pragma omp parallel for schedule(dynamic)
    for (auto i = GetLv1StartBucket(); i < GetLv1EndBucket(); ++i) {
      Lv2Sort(&thread_status, i, omp_get_thread_num());
    }
  }

 public:
  void Run();
};

#endif  //MEGAHIT_BASE_ENGINE_H