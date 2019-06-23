//
// Created by vout on 6/23/19.
//

#include "base_engine.h"
#include <omp.h>


namespace {

/**
 * Try to fully use available memory to fit as much lv1 items as possible
 */
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
      num_lv2_items -= 0.95;
      continue;
    }

    num_lv1_items = (mem_avail - mem_for_lv2) / BaseSequenceSortingEngine::kLv1BytePerItem;
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
      num_lv1_items = (mem_avail - bytes_per_lv2_item * num_lv2_items) / BaseSequenceSortingEngine::kLv1BytePerItem;
    } else {
      break;
    }
  }

  return {num_lv1_items, num_lv2_items};
}

}

void BaseSequenceSortingEngine::AdjustMemory() {
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

void BaseSequenceSortingEngine::Run() {
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
  Lv0PrepareThreadPartition();
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

void BaseSequenceSortingEngine::Lv0PrepareThreadPartition() {
  thread_meta_.resize(n_threads_);
  for (unsigned t = 0; t < n_threads_; ++t) {
    ThreadMeta &rp = thread_meta_[t];
    // distribute reads to partitions
    int64_t average = meta_.num_sequences / n_threads_;
    rp.seq_from = t * average;
    rp.seq_to = t < n_threads_ - 1 ? (t + 1) * average : meta_.num_sequences;
    rp.offset_base = Lv0EncodeDiffBase(rp.seq_from);
  }

  for (unsigned i = 0; i < kNumBuckets; ++i) {
    ori_bucket_id_[i] = i;
    bucket_rank_[i] = i;
  }
}

void BaseSequenceSortingEngine::Lv0ReorderBuckets() {
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
    auto old_rp_bucket_sizes = thread_meta_[tid].bucket_sizes;
    for (unsigned i = 0; i < kNumBuckets; ++i) {
      thread_meta_[tid].bucket_sizes[i] = old_rp_bucket_sizes[tmp_v[i].second];
    }
  }
}

unsigned BaseSequenceSortingEngine::Lv1FindEndBuckets(unsigned start_bucket) {
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

void BaseSequenceSortingEngine::Lv1ComputeThreadBegin() {
  // set the bucket begin for the first thread
  thread_meta_[0].bucket_begin[lv1_start_bucket_] = 0;
  for (unsigned b = lv1_start_bucket_ + 1; b < lv1_end_bucket_; ++b) {
    thread_meta_[0].bucket_begin[b] = thread_meta_[0].bucket_begin[b - 1] + bucket_sizes_[b - 1];  // accumulate
  }

  // then for the remaining threads
  for (unsigned t = 1; t < n_threads_; ++t) {
    auto &current_thread_begin = thread_meta_[t].bucket_begin;
    auto &prev_thread_begin = thread_meta_[t - 1].bucket_begin;
    auto &prev_thread_sizes = thread_meta_[t - 1].bucket_sizes;
    for (unsigned b = lv1_start_bucket_; b < lv1_end_bucket_; ++b) {
      current_thread_begin[b] = prev_thread_begin[b] + prev_thread_sizes[b];
    }
  }
}

void BaseSequenceSortingEngine::Lv0CalcBucketSizeLaunchMt() {
#pragma omp parallel for
  for (unsigned t = 0; t < n_threads_; ++t) {
    auto &thread_meta = thread_meta_[t];
    Lv0CalcBucketSize(thread_meta.seq_from, thread_meta.seq_to, &thread_meta.bucket_sizes);
  }
  std::fill(bucket_sizes_.begin(), bucket_sizes_.end(), 0);

  for (unsigned t = 0; t < n_threads_; ++t) {
    for (unsigned b = 0; b < kNumBuckets; ++b) {
      bucket_sizes_[b] += thread_meta_[t].bucket_sizes[b];
    }
  }
}

void BaseSequenceSortingEngine::Lv1FetchAndSortLaunchMt() {
  Lv2ThreadStatus thread_status{};
  thread_status.thread_offset.resize(n_threads_, -1);
  thread_status.rank.resize(n_threads_, 0);
  omp_set_num_threads(n_threads_);
#pragma omp parallel for schedule(dynamic)
  for (auto i = GetLv1StartBucket(); i < GetLv1EndBucket(); ++i) {
    Lv2Sort(&thread_status, i, omp_get_thread_num());
  }
}

void BaseSequenceSortingEngine::Lv2Sort(BaseSequenceSortingEngine::Lv2ThreadStatus *thread_status, unsigned b, int tid) {
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

void BaseSequenceSortingEngine::Lv1FillOffsetsLaunchMt() {
  lv1_special_offsets_.clear();
  Lv1ComputeThreadBegin();

#pragma omp parallel for
  for (unsigned t = 0; t < n_threads_; ++t) {
    OffsetFiller filler(this, lv1_start_bucket_, lv1_end_bucket_, thread_meta_[t]);
    Lv1FillOffsets(filler, thread_meta_[t].seq_from, thread_meta_[t].seq_to);
  }
}