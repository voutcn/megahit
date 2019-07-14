#ifndef MEGAHIT_BASE_ENGINE_H
#define MEGAHIT_BASE_ENGINE_H

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <mutex>
#include <vector>

#include "kmsort_selector.h"
#include "utils/utils.h"

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

  struct MemoryStat {
    int64_t num_sequences;
    int64_t memory_for_data;
    int64_t words_per_lv2;
    int64_t aux_words_per_lv2;
  };

 public:
  BaseSequenceSortingEngine(int64_t mem, int mem_flag, int n_threads)
      : host_mem_(mem), mem_flag_(mem_flag), n_threads_(n_threads){};
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

 private:
  //  Sequence divided into n_threads_ partitions
  //  Each thread read one partition when filling Lv1 offsets, but it will fill
  //  multiple buckets
  //
  //  SeqView
  //  <--------------p0--------------><--------------p1-------------->
  //  |t0---------------------------->|t1---------------------------->
  //
  //  BucketView:
  //  <------b0------><------b1------><------b2------><------b3------>
  //  |t0---->|t1---->|t0---->|t1---->|t0---->|t1---->|t0---->|t1---->
  //  |                                       |
  //  \ bucket_begin[0] in thread_meta_[0]    \ bucket_begin[2] in
  //  thread_meta_[1]

  struct ThreadMeta {
    int64_t seq_from, seq_to;  // start and end IDs of this sequence partition
                               // (end is exclusive)
    std::array<int64_t, kNumBuckets> bucket_sizes;
    std::array<int64_t, kNumBuckets> bucket_begin;
    int64_t offset_base;  // the initial offset globals.lv1_items
  };

  std::vector<ThreadMeta> thread_meta_;
  std::array<int64_t, kNumBuckets> bucket_sizes_{};
  std::array<unsigned, kNumBuckets> bucket_real_id{};
  std::array<unsigned, kNumBuckets> bucket_rank_{};

 protected:
  using SubstrPtr = std::vector<uint32_t>::iterator;
  /**
   * For derived class to fill offsets
   */
  class OffsetFiller {
   public:
    OffsetFiller(BaseSequenceSortingEngine *engine, unsigned start_bucket,
                 unsigned end_bucket, const ThreadMeta &sp)
        : engine_(engine) {
      std::fill(prev_full_offsets_.begin() + start_bucket,
                prev_full_offsets_.begin() + end_bucket, sp.offset_base);
      std::copy(sp.bucket_begin.begin() + start_bucket,
                sp.bucket_begin.begin() + end_bucket,
                bucket_index_.begin() + start_bucket);
    }

    void WriteNextOffset(unsigned bucket, int64_t full_offset) {
      assert(IsHandling(bucket));
      bucket = engine_->bucket_rank_[bucket];
      int64_t diff = full_offset - prev_full_offsets_[bucket];
      int64_t index = bucket_index_[bucket]++;
      engine_->WriteOffset(index, diff, full_offset);
      prev_full_offsets_[bucket] = full_offset;
    }

    bool IsHandling(unsigned bucket) const {
      return engine_->cur_lv1_buckets_[bucket];
    }

   private:
    BaseSequenceSortingEngine *engine_;
    std::array<int64_t, kNumBuckets> prev_full_offsets_;
    std::array<int64_t, kNumBuckets> bucket_index_;
  };

  /**
   * Iterator to fetch full offset from given buckets
   */
  class OffsetFetcher {
   public:
    OffsetFetcher(BaseSequenceSortingEngine *engine, unsigned start_bucket,
                  unsigned end_bucket)
        : engine_(engine), end_bucket_(end_bucket), cur_bucket_(start_bucket) {
      cur_thread_id_ = 0;
      cur_item_index_ = 0;
      cur_full_offset_ = engine_->thread_meta_[cur_thread_id_].offset_base;
      cur_n_items_ =
          engine_->thread_meta_[cur_thread_id_].bucket_sizes[cur_bucket_];
      diff_iter_ = engine_->lv1_offsets_.cbegin() +
                   engine_->thread_meta_[0].bucket_begin[cur_bucket_];
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
        cur_full_offset_ = engine_->thread_meta_[cur_thread_id_].offset_base;
        cur_item_index_ = 0;
        cur_n_items_ =
            engine_->thread_meta_[cur_thread_id_].bucket_sizes[cur_bucket_];
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

 private:
  void WriteOffset(size_t index, int64_t diff, int64_t full_offset) {
    if (diff > kDifferentialLimit) {
      int64_t bucket = kDifferentialLimit + 1 + AddSpecialOffset(full_offset);
      if (bucket >= std::numeric_limits<uint32_t>::max()) {
        xfatal("Too many large difference items!");
      }
      lv1_offsets_[index] = bucket;
    } else {
      assert(diff >= 0);
      lv1_offsets_[index] = static_cast<int32_t>(diff);
    }
  }

  int64_t AddSpecialOffset(int64_t offset) {
    std::lock_guard<std::mutex> lk(special_item_lock_);
    int64_t ret = lv1_special_offsets_.size();
    lv1_special_offsets_.push_back(offset);
    return ret;
  }

  int64_t GetSpecialOffset(uint32_t bucket) const {
    assert(bucket > kDifferentialLimit);
    return lv1_special_offsets_[bucket - (kDifferentialLimit + 1)];
  }

 private:
  /**
   * Memory related data
   */
  int64_t host_mem_{};
  int mem_flag_{};
  unsigned n_threads_{};
  MemoryStat meta_{};  // set by Initialize return
  std::function<void(uint32_t *, int64_t)>
      substr_sort_;  // set after Initialize

  /**
   * Interfaces used by `Run` and must be implemented in derived class
   */
 public:
  virtual MemoryStat Initialize() = 0;

 protected:
  virtual int64_t Lv0EncodeDiffBase(int64_t) = 0;
  virtual void Lv0CalcBucketSize(int64_t seq_from, int64_t seq_to,
                                 std::array<int64_t, kNumBuckets> *out) = 0;
  virtual void Lv1FillOffsets(OffsetFiller &filler, int64_t seq_from,
                              int64_t seq_to) = 0;
  virtual void Lv2ExtractSubString(OffsetFetcher &fetcher,
                                   SubstrPtr substr_ptr) = 0;
  virtual void Lv2Postprocess(int64_t start_index, int64_t end_index,
                              int thread_id, uint32_t *substr_ptr) = 0;
  virtual void Lv0Postprocess() = 0;

 private:
  void AdjustMemory();
  void Lv0PrepareThreadPartition();
  void Lv0CalcBucketSizeLaunchMt();
  void Lv0ReorderBuckets();
  unsigned Lv1FindEndBuckets(unsigned start_bucket);
  void Lv1ComputeThreadBegin();
  inline void Lv1FillOffsetsLaunchMt();
  void Lv1FetchAndSortLaunchMt();

  struct Lv2ThreadStatus {
    std::vector<int64_t> thread_offset;
    std::vector<int64_t> rank;
    std::mutex mutex;
    int64_t acc{0};
    int seen{0};
  };

  void Lv2Sort(Lv2ThreadStatus *thread_status, unsigned b, int tid);

 public:
  void Run();
};

#endif  // MEGAHIT_BASE_ENGINE_H