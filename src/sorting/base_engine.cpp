//
// Created by vout on 6/23/19.
//

#include "base_engine.h"

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
