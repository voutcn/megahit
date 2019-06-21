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

#ifndef CX1_H__
#define CX1_H__

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include <omp.h>

#include "utils/utils.h"
#include "sorting.h"


template<typename TGlobal, unsigned NBuckets>
struct CX1 {
 public:
  typedef TGlobal global_data_t;
  // other settings, don't change
  static const int kLv1BytePerItem = 4;  // 32-bit differatial offset
  static const uint64_t kSpDiffMaxNum = (1ULL << 32) - 1;
  static const int64_t kDifferentialLimit = (1ULL << 31) - 1;

  struct ReadPartition {
    // local data for each read partition (i.e. a subrange of input reads)
    global_data_t *globals;
    int64_t rp_start_id, rp_end_id;  // start and end IDs of this read partition (end is exclusive)
    std::array<int64_t, NBuckets> rp_bucket_sizes;
    std::array<int64_t, NBuckets> rp_bucket_offsets;
    int64_t rp_lv1_differential_base;  // the initial offset globals.lv1_items
  };
  // param: must be set
  global_data_t *g_;

 private:
  // may change as cx1 goes
  int64_t lv1_num_items_{};
  unsigned lv1_start_bucket_{}, lv1_end_bucket_{};
  std::vector<bool> cur_lv1_buckets_;
  std::vector<int64_t> lv1_special_offsets_;
  std::vector<int32_t> lv1_offsets_;
  std::mutex special_item_lock_;

  std::array<int64_t, NBuckets> bucket_sizes_;
  std::array<int, NBuckets> ori_bucket_id_;
  std::array<int, NBuckets> bucket_rank_;
  std::vector<ReadPartition> rp_;

 public:
  int64_t GetLv1NumItems() const { return lv1_num_items_; }
  unsigned GetLv1StartBucket() const { return lv1_start_bucket_; }
  unsigned GetLv1EndBucket() const { return lv1_end_bucket_; }

  bool HandlingBucket(unsigned bucket) const {
    return cur_lv1_buckets_[bucket];
  }

  void WriteOffset(size_t index, int64_t diff, int64_t full_offset) {
    if (diff > kDifferentialLimit) {
      lv1_offsets_[index] = -1 - AddSpecialOffset(full_offset);
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

  int64_t GetSpecialOffset(size_t index) const {
    return lv1_special_offsets_[index];
  }

  std::vector<int32_t>::const_iterator GetLv1Iterator(unsigned bucket) const {
    return lv1_offsets_.cbegin() + GetReadPartition(0).rp_bucket_offsets[bucket];
  }

  uint32_t *Lv1DataPtr() {
    return reinterpret_cast<uint32_t *>(lv1_offsets_.data());
  }

  const std::array<int64_t, NBuckets> &GetBucketSizes() const {
    return bucket_sizes_;
  }

  int GetBucketRank(unsigned bucket) const { return bucket_rank_[bucket]; }

  const ReadPartition &GetReadPartition(unsigned tid) const {
    return rp_[tid];
  }

 private:
  // set by initialization
  int64_t num_items_{};
  unsigned num_cpu_threads_{};

  // set according to memory
  int64_t max_mem_remain_{};
  size_t words_per_sorting_item_{};
  size_t extra_words_{};
  int64_t max_lv1_items_{};
  std::function<void(uint32_t *, int64_t)> substr_sort_;

 public:
  void SetNumItems(int64_t num) { num_items_ = num; }
  void SetNumCpuThreads(unsigned n) { num_cpu_threads_ = n; }
  void SetMaxLv1Lv2Items(int64_t n_lv1, int64_t n_lv2, int64_t words_per_lv2, int64_t extra_words = 0) {
    max_lv1_items_ = n_lv1;
    words_per_sorting_item_ = words_per_lv2;
    extra_words_ = extra_words;
    auto words_required = n_lv1 + n_lv2 * words_per_lv2;
    max_mem_remain_ = words_required * kLv1BytePerItem;
    lv1_offsets_.reserve(words_required);
    lv1_offsets_.resize(words_required);
    substr_sort_ = SelectSortingFunc(words_per_sorting_item_ - extra_words_, extra_words_);
  }
  size_t GetWordsPerSortingItem() const { return words_per_sorting_item_; }

  // === functions to specify a CX1 instance ===
 public:
  virtual int64_t encode_lv1_diff_base_func_(int64_t, global_data_t &) = 0;
  virtual void prepare_func_(global_data_t &) = 0;  // num_items_, num_cpu_threads_ and num_output_threads_ must be set here
  virtual void lv0_calc_bucket_size_func_(ReadPartition *) = 0;
  virtual void init_global_and_set_cx1_func_(global_data_t &) = 0;  // xxx set here
  virtual void lv1_fill_offset_func_(ReadPartition *) = 0;
  virtual void lv2_extract_substr_(unsigned bucket_from, unsigned bucket_to, TGlobal &g, uint32_t *substr_ptr) = 0;
  virtual void output_(int64_t start_index, int64_t end_index, int thread_id, TGlobal &g, uint32_t *substrings) = 0;
  virtual void post_proc_func_(global_data_t &) = 0;

  CX1() = default;

  // === single thread functions ===
  static void adjust_mem(int64_t mem_avail, int64_t bytes_per_sorting_item, int64_t min_lv1_items,
                         int64_t min_sorting_items, int64_t max_sorting_items, int64_t &max_lv1_items,
                         int64_t &num_sorting_items) {
    num_sorting_items = max_sorting_items;

    while (num_sorting_items >= min_sorting_items) {
      int64_t mem_sorting_items = bytes_per_sorting_item * num_sorting_items;
      xinfo(
          "Adjusting memory layout: max_lv1_items=%lld, num_sorting_items=%lld, "
          "mem_sorting_items=%lld, mem_avail=%lld\n",
          max_lv1_items, num_sorting_items, mem_sorting_items, mem_avail);

      if (mem_avail < mem_sorting_items) {
        num_sorting_items = num_sorting_items * 0.95;
        continue;
      }

      max_lv1_items = (mem_avail - mem_sorting_items) / kLv1BytePerItem;
      if (max_lv1_items < min_lv1_items || max_lv1_items < num_sorting_items) {
        num_sorting_items *= 0.95;
      } else {
        break;
      }
    }

    if (num_sorting_items < min_sorting_items) {
      xfatal("No enough memory to process CX1.\n");
    }

    // --- adjust num_sorting_items to fit more lv1 item ---
    // TODO: 4 is arbitrary chosen, not fine tune
    while (num_sorting_items * 4 > max_lv1_items) {
      if (num_sorting_items * 0.95 >= min_sorting_items) {
        num_sorting_items *= 0.95;
        max_lv1_items = (mem_avail - bytes_per_sorting_item * num_sorting_items) / kLv1BytePerItem;
      } else {
        break;
      }
    }
  }

  inline void prepare_rp_and_bp_() {  // call after prepare_func_
    rp_.resize(num_cpu_threads_);

    for (unsigned t = 0; t < num_cpu_threads_; ++t) {
      struct ReadPartition &rp = rp_[t];
      rp.globals = g_;
      // distribute reads to partitions
      int64_t average = num_items_ / num_cpu_threads_;
      rp.rp_start_id = t * average;
      rp.rp_end_id = t < num_cpu_threads_ - 1 ? (t + 1) * average : num_items_;
      rp.rp_lv1_differential_base = encode_lv1_diff_base_func_(rp.rp_start_id, *g_);
    }

    for (unsigned i = 0; i < NBuckets; ++i) {
      ori_bucket_id_[i] = i;
      bucket_rank_[i] = i;
    }
  }

  inline void reorder_buckets_() {
    std::vector<std::pair<int64_t, int> > tmp_v(NBuckets);

    for (unsigned i = 0; i < NBuckets; ++i) {
      tmp_v[i] = std::make_pair(bucket_sizes_[i], i);
    }

    std::sort(tmp_v.rbegin(), tmp_v.rend());

    for (unsigned i = 0; i < NBuckets; ++i) {
      bucket_sizes_[i] = tmp_v[i].first;
      ori_bucket_id_[i] = tmp_v[i].second;
      bucket_rank_[tmp_v[i].second] = i;
    }

    for (unsigned tid = 0; tid < num_cpu_threads_; ++tid) {
      auto old_rp_bucket_sizes = rp_[tid].rp_bucket_sizes;
      for (unsigned i = 0; i < NBuckets; ++i) {
        rp_[tid].rp_bucket_sizes[i] = old_rp_bucket_sizes[tmp_v[i].second];
      }
    }
  }

  inline int find_end_buckets_with_rank_(unsigned start_bucket, unsigned end_limit, int64_t mem_limit,
                                         unsigned bytes_per_sorting_items, int64_t &num_items) {
    num_items = 0;
    unsigned end_bucket = start_bucket;
    unsigned used_threads = 0;
    int64_t mem_sorting_items = 0;

    cur_lv1_buckets_.resize(NBuckets);
    std::fill(cur_lv1_buckets_.begin(), cur_lv1_buckets_.end(), false);

    while (end_bucket < end_limit) {
      if (used_threads < num_cpu_threads_) {
        mem_sorting_items += bytes_per_sorting_items * bucket_sizes_[end_bucket];
        ++used_threads;
      }

      if (mem_sorting_items + (num_items + bucket_sizes_[end_bucket]) * kLv1BytePerItem > mem_limit) {
        return end_bucket;
      }

      num_items += bucket_sizes_[end_bucket];
      cur_lv1_buckets_[ori_bucket_id_[end_bucket]] = true;
      ++end_bucket;
    }

    return end_limit;
  }

  inline void lv1_compute_offset_() {
    // compute "global" (thread 0) offsets first
    auto &offsets = rp_[0].rp_bucket_offsets;
    offsets[lv1_start_bucket_] = 0;

    for (unsigned b = lv1_start_bucket_ + 1; b < lv1_end_bucket_; ++b) {
      offsets[b] = offsets[b - 1] + bucket_sizes_[b - 1];  // accumulate
    }

    // then for each read partition
    for (unsigned t = 1; t < num_cpu_threads_; ++t) {
      auto &this_offsets = rp_[t].rp_bucket_offsets;
      auto &prev_offsets = rp_[t - 1].rp_bucket_offsets;
      auto &sizes = rp_[t - 1].rp_bucket_sizes;

      for (unsigned b = lv1_start_bucket_; b < lv1_end_bucket_; ++b) {
        this_offsets[b] = prev_offsets[b] + sizes[b];
      }
    }
  }

  // === multi-thread wrappers ====
  inline void lv0_calc_bucket_size_mt_() {
#pragma omp parallel for
    for (unsigned t = 0; t < num_cpu_threads_; ++t) {
      lv0_calc_bucket_size_func_(&rp_[t]);
    }
    // sum up readpartitions bucketsizes to form global bucketsizes
    std::fill(bucket_sizes_.begin(), bucket_sizes_.end(), 0);

    // the array accesses in this loop are optimized by the compiler??
    for (unsigned t = 0; t < num_cpu_threads_; ++t) {
      for (unsigned b = 0; b < NBuckets; ++b) {
        bucket_sizes_[b] += rp_[t].rp_bucket_sizes[b];
      }
    }
  }

  inline void lv1_fill_offset_mt_() {
    lv1_special_offsets_.clear();
    lv1_compute_offset_();

    // create threads
#pragma omp parallel for
    for (unsigned t = 0; t < num_cpu_threads_; ++t) {
      lv1_fill_offset_func_(&rp_[t]);
    }
    lv1_compute_offset_();
  }

  struct Lv2ThreadStatus {
    std::vector<int64_t> thread_offset;
    std::vector<int64_t> rank;
    std::mutex mutex;
    int64_t acc{0};
    int seen{0};
  };

  void sort_bucket(Lv2ThreadStatus *thread_status, unsigned b, int tid) {
    if (thread_status->thread_offset[tid] == -1) {
      std::lock_guard<std::mutex> lk(thread_status->mutex);
      thread_status->thread_offset[tid] = thread_status->acc;
      thread_status->acc += GetBucketSizes()[b];
      thread_status->rank[tid] = thread_status->seen;
      thread_status->seen++;
    }

    if (GetBucketSizes()[b] == 0) {
      return;
    }

    size_t offset = GetLv1NumItems() + thread_status->thread_offset[tid] * GetWordsPerSortingItem();
    auto substr_ptr = Lv1DataPtr() + offset;
    lv2_extract_substr_(b, b + 1, *g_, substr_ptr);
    substr_sort_(substr_ptr, GetBucketSizes()[b]);
    output_(0, GetBucketSizes()[b], tid, *g_, substr_ptr);
  }

  void lv1_fetch_and_sort_mt_() {
    Lv2ThreadStatus thread_status{};
    thread_status.thread_offset.resize(num_cpu_threads_, -1);
    thread_status.rank.resize(num_cpu_threads_, 0);
    omp_set_num_threads(num_cpu_threads_);
#pragma omp parallel for schedule(dynamic)
    for (auto i = GetLv1StartBucket(); i < GetLv1EndBucket(); ++i) {
      sort_bucket(&thread_status, i, omp_get_thread_num());
    }
  }

  // entry
  void run() {
    SimpleTimer lv0_timer;
    // read input & prepare
    lv0_timer.reset();
    lv0_timer.start();
    xinfo("Preparing data...\n");

    prepare_func_(*g_);

    lv0_timer.stop();
    xinfo("Preparing data... Done. Time elapsed: %.4f\n", lv0_timer.elapsed());
    lv0_timer.reset();
    lv0_timer.start();
    xinfo("Preparing partitions and initialing global data...\n");

    // prepare rp bp and op
    prepare_rp_and_bp_();
    // calc bucket size
    lv0_calc_bucket_size_mt_();
    // init global datas
    init_global_and_set_cx1_func_(*g_);

    reorder_buckets_();
    lv0_timer.stop();
    xinfo("Preparing partitions and initialing global data... Done. Time elapsed: %.4f\n", lv0_timer.elapsed());

    lv0_timer.reset();
    lv0_timer.start();
    xinfo("Start main loop...\n");
    // === start main loop ===
    int lv1_iteration = 0;
    lv1_start_bucket_ = 0;

    while (lv1_start_bucket_ < NBuckets) {
      SimpleTimer lv1_timer;
      lv1_iteration++;
      // --- finds the bucket range for this iteration ---
      lv1_end_bucket_ = find_end_buckets_with_rank_(lv1_start_bucket_, NBuckets, max_mem_remain_,
                                                    words_per_sorting_item_ * sizeof(uint32_t), lv1_num_items_);

      if (lv1_num_items_ == 0 && lv1_end_bucket_ < NBuckets) {
        fprintf(stderr, "Bucket %d too large for lv1: %lld > %lld\n", lv1_end_bucket_,
                static_cast<long long>(bucket_sizes_[lv1_end_bucket_]), static_cast<long long>(max_lv1_items_));
        exit(1);
      }

      lv1_timer.reset();
      lv1_timer.start();
      xinfo("Lv1 scanning from bucket %d to %d\n", lv1_start_bucket_, lv1_end_bucket_);

      // --- scan to fill offset ---
      lv1_fill_offset_mt_();

      if (lv1_special_offsets_.size() > kSpDiffMaxNum) {
        fprintf(stderr, "Too many large diff items (%lu) from in buckets [%d, %d)\n", lv1_special_offsets_.size(),
                lv1_start_bucket_, lv1_end_bucket_);
        exit(1);
      }

      lv1_timer.stop();
      xinfo("Lv1 scanning done. Large diff: %lu. Time elapsed: %.4f\n",
            lv1_special_offsets_.size(),
            lv1_timer.elapsed());
      lv1_timer.reset();
      lv1_timer.start();
      lv1_fetch_and_sort_mt_();
      lv1_timer.stop();
      xinfo("Lv1 fetching & sorting done. Time elapsed: %.4f\n", lv1_timer.elapsed());
      lv1_start_bucket_ = lv1_end_bucket_;
    }

    lv0_timer.stop();
    xinfo("Main loop done. Time elapsed: %.4f\n", lv0_timer.elapsed());
    lv0_timer.reset();
    lv0_timer.start();
    xinfo("Postprocessing...\n");
    post_proc_func_(*g_);
    lv0_timer.stop();
    xinfo("Postprocess done. Time elapsed: %.4f\n", lv0_timer.elapsed());
  }
};

#endif  // CX1_H__