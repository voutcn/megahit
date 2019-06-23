#ifndef MEGAHIT_BASE_ENGINE_H
#define MEGAHIT_BASE_ENGINE_H

#include <cassert>
#include <cstdint>
#include <algorithm>
#include <vector>
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
  static const int64_t kDifferentialLimit = (1ULL << 31) - 1;

  struct ReadPartition {
    // local data for each read partition (i.e. a sub-range of input reads)
    int64_t rp_start_id, rp_end_id;  // start and end IDs of this read partition (end is exclusive)
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
  // change during in different Lv1 stage
  unsigned lv1_start_bucket_{}, lv1_end_bucket_{};
  int64_t lv1_num_items_{};
  std::vector<bool> cur_lv1_buckets_;
  std::vector<int64_t> lv1_special_offsets_;
  std::vector<uint32_t> lv1_offsets_;
  std::mutex special_item_lock_;

  std::array<int64_t, kNumBuckets> bucket_sizes_{};
  std::array<int, kNumBuckets> ori_bucket_id_{};
  std::array<int, kNumBuckets> bucket_rank_{};
  std::vector<ReadPartition> rp_;

  /**
   * APIs used by derived classes
   */
 protected:
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

    bool HasNext() { return cur_bucket_ < end_bucket_; }
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

  bool HandlingBucket(unsigned bucket) const {
    return cur_lv1_buckets_[bucket];
  }

  OffsetIterator GetOffsetIterator(unsigned start_bucket, unsigned end_bucket) {
    return OffsetIterator(this, start_bucket, end_bucket);
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

  /**
   * Memory related data
   */
 private:
  // set by initialization
  int64_t host_mem_{};
  int mem_flag_{};
  unsigned n_threads_{};
  Meta meta_{};

  // set according to memory
  std::function<void(uint32_t *, int64_t)> substr_sort_;

  /**
   * Interfaces used by `Run` and must be implemented in derived class
   */
 public:
  virtual Meta Initialize() = 0;

 protected:
  virtual int64_t Lv0EncodeDiffBase(int64_t) = 0;
  virtual void Lv0CalcBucketSize(ReadPartition *) = 0;
  virtual void Lv1FillOffsets(ReadPartition *) = 0;
  virtual void Lv2ExtractSubString(unsigned bucket_from, unsigned bucket_to, uint32_t *substr_ptr) = 0;
  virtual void Lv2Postprocess(int64_t start_index, int64_t end_index, int thread_id, uint32_t *substrings) = 0;
  virtual void Lv0Postprocess() = 0;

 private:
  void AdjustMemory() {
    int64_t max_bucket_size = *std::max_element(bucket_sizes_.begin(), bucket_sizes_.end());
    int64_t total_bucket_size = 0;
    int num_non_empty = 0;
    for (unsigned i = 0; i < kNumBuckets; ++i) {
      total_bucket_size += bucket_sizes_[i];
      if (bucket_sizes_[i] > 0) {
        num_non_empty++;
      }
    }

    int64_t est_lv2_items = std::max(3 * total_bucket_size / std::max(1, num_non_empty) * n_threads_,
                                     max_bucket_size);

    const int64_t mem_remained = host_mem_ - meta_.memory_for_data;
    const int64_t min_lv1_items = total_bucket_size / (kMaxLv1ScanTime - 0.5);
    const int64_t lv2_bytes_per_item = meta_.words_per_lv2 * sizeof(uint32_t);
    std::pair<int64_t, int64_t> n_items;

    if (mem_flag_ == 1) {
      // auto set memory
      int64_t est_lv1_items = total_bucket_size / (kDefaultLv1ScanTime - 0.5);
      est_lv1_items = std::max(est_lv1_items, max_bucket_size);
      int64_t mem_needed = est_lv1_items * kLv1BytePerItem + est_lv2_items * lv2_bytes_per_item;
      if (mem_needed > mem_remained) {
        n_items = AdjustItemNumbers(mem_remained, lv2_bytes_per_item, min_lv1_items, max_bucket_size,
                                    est_lv2_items);
      } else {
        n_items = {est_lv1_items, est_lv2_items};
      }
    } else if (mem_flag_ == 0) {
      // min memory
      int64_t est_lv1_items = total_bucket_size / (kMaxLv1ScanTime - 0.5);
      est_lv1_items = std::max(est_lv1_items, max_bucket_size);
      int64_t mem_needed = est_lv1_items * kLv1BytePerItem + est_lv2_items * lv2_bytes_per_item;

      if (mem_needed > mem_remained) {
        n_items = AdjustItemNumbers(mem_remained, lv2_bytes_per_item, min_lv1_items, max_bucket_size,
                                    est_lv2_items);
      } else {
        n_items = AdjustItemNumbers(mem_needed, lv2_bytes_per_item, min_lv1_items, max_bucket_size,
                                    est_lv2_items);
      }
    } else {
      // use all
      n_items = AdjustItemNumbers(mem_remained, lv2_bytes_per_item, min_lv1_items, max_bucket_size,
                                  est_lv2_items);
    }

    if (n_items.first < min_lv1_items) {
      xfatal("No enough memory");
    }

    if (n_items.first > total_bucket_size) {
      n_items.first = total_bucket_size;
    }

    auto words_required = n_items.first + n_items.second * meta_.words_per_lv2;
    lv1_offsets_.reserve(words_required);
    lv1_offsets_.resize(words_required);
    substr_sort_ = SelectSortingFunc(meta_.words_per_lv2 - meta_.aux_words_per_lv2, meta_.aux_words_per_lv2);

    xinfo("Lv1 items: %lld, Lv2 items: %lld\n", n_items.first, n_items.second);
    xinfo("Memory of derived class: %lld, Memory for Lv1+Lv2: %lld\n",
          meta_.memory_for_data, words_required * kLv1BytePerItem);
  }

  static std::pair<int64_t, int64_t> AdjustItemNumbers(
      int64_t mem_avail, int64_t bytes_per_lv2_item, int64_t min_lv1_items,
      int64_t min_lv2_items, int64_t max_lv2_items) {

    int64_t num_lv1_items = 0;
    int64_t num_lv2_items = max_lv2_items;
    while (num_lv2_items >= min_lv2_items) {
      int64_t mem_for_lv2 = bytes_per_lv2_item * num_lv2_items;

      xinfo("Adjusting memory layout: num_lv1_items=%lld, num_lv2_items=%lld, "
            "mem_sorting_items=%lld, mem_avail=%lld\n",
            num_lv1_items, num_lv2_items, mem_for_lv2, mem_avail);

      if (mem_avail < mem_for_lv2) {
        num_lv2_items *= 0.95;
        continue;
      }

      num_lv1_items = (mem_avail - mem_for_lv2) / kLv1BytePerItem;
      if (num_lv1_items < min_lv1_items || num_lv1_items < num_lv2_items) {
        num_lv2_items *= 0.95;
      } else {
        break;
      }
    }

    if (num_lv2_items < min_lv2_items) {
      xfatal("No enough memory to perform sequence sorting.\n");
    }

    // --- adjust num_lv2_items to fit more lv1 item ---
    while (num_lv2_items * 4 > num_lv1_items) {
      if (num_lv2_items * 0.95 >= min_lv2_items) {
        num_lv2_items *= 0.95;
        num_lv1_items = (mem_avail - bytes_per_lv2_item * num_lv2_items) / kLv1BytePerItem;
      } else {
        break;
      }
    }

    return {num_lv1_items, num_lv2_items};
  }

  void Lv0PrepareReadPartition() {
    rp_.resize(n_threads_);
    for (unsigned t = 0; t < n_threads_; ++t) {
      ReadPartition &rp = rp_[t];
      // distribute reads to partitions
      int64_t average = meta_.num_sequences / n_threads_;
      rp.rp_start_id = t * average;
      rp.rp_end_id = t < n_threads_ - 1 ? (t + 1) * average : meta_.num_sequences;
      rp.rp_lv1_differential_base = Lv0EncodeDiffBase(rp.rp_start_id);
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
  void Run() {
    SimpleTimer lv0_timer;
    // read input & prepare
    lv0_timer.reset();
    lv0_timer.start();
    xinfo("Preparing data...\n");

    meta_ = Initialize();

    lv0_timer.stop();
    xinfo("Preparing data... Done. Time elapsed: %.4f\n", lv0_timer.elapsed());
    lv0_timer.reset();
    lv0_timer.start();
    xinfo("Preparing partitions and calculating bucket sizes...\n");

    // prepare rp bp and op
    Lv0PrepareReadPartition();
    // calc bucket size
    Lv0CalcBucketSizeLaunchMt();
    Lv0ReorderBuckets();
    AdjustMemory();
    lv0_timer.stop();
    xinfo("Preparing partitions and calculating bucket sizes... Done. Time elapsed: %.4f\n", lv0_timer.elapsed());

    lv0_timer.reset();
    lv0_timer.start();
    xinfo("Start main loop...\n");
    int lv1_iteration = 0;
    lv1_start_bucket_ = 0;

    while (lv1_start_bucket_ < kNumBuckets) {
      SimpleTimer lv1_timer;
      lv1_iteration++;
      // --- finds the bucket range for this iteration ---
      lv1_end_bucket_ = Lv1FindEndBuckets(lv1_start_bucket_);
      assert(lv1_start_bucket_ < lv1_end_bucket_);

      lv1_timer.reset();
      lv1_timer.start();
      xinfo("Lv1 scanning from bucket %d to %d\n", lv1_start_bucket_, lv1_end_bucket_);

      // --- scan to fill offset ---
      Lv1FillOffsetsLaunchMt();

      lv1_timer.stop();
      xinfo("Lv1 scanning done. Large diff: %lu. Time elapsed: %.4f\n",
            lv1_special_offsets_.size(),
            lv1_timer.elapsed());
      lv1_timer.reset();
      lv1_timer.start();
      Lv1FetchAndSortLaunchMt();
      lv1_timer.stop();
      xinfo("Lv1 fetching & sorting done. Time elapsed: %.4f\n", lv1_timer.elapsed());
      lv1_start_bucket_ = lv1_end_bucket_;
    }

    lv0_timer.stop();
    xinfo("Main loop done. Time elapsed: %.4f\n", lv0_timer.elapsed());
    lv0_timer.reset();
    lv0_timer.start();
    xinfo("Postprocessing...\n");
    Lv0Postprocess();
    lv0_timer.stop();
    xinfo("Postprocess done. Time elapsed: %.4f\n", lv0_timer.elapsed());
  }
};

#endif  //MEGAHIT_BASE_ENGINE_H