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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include <stdint.h>
#include <vector>
#include <algorithm>

#include "mem_file_checker-inl.h"
#include "utils.h"

/**
 * @brief    an CX1 engine
 * @details  use CX1 algorithm to do all kinds of things related to substring sorting
 *
 * @tparam global_data_type the type of global datas used in a specified CX1 engine
 *                          must contain a member `CX1* cx1`, where cx1->g_ is itself,
 *                          that enables interactions with CX1 functions
 * @tparam kNumBuckets      number of buckets
 */
template <typename global_data_type, int kNumBuckets>
struct CX1 {
    typedef global_data_type global_data_t;
    static const int kCX1Verbose = 3; // tunable
    // other settings, don't change
    static const int kGPUBytePerItem = 16; // key & value, 4 byte each. double for radix sort internal buffer
    static const int kLv1BytePerItem = 4; // 32-bit differatial offset
    static const uint64_t kSpDiffMaxNum = (1ULL << 32) - 1;
    static const int64_t kDifferentialLimit = (1ULL << 31) - 1;

    struct readpartition_data_t {
        // local data for each read partition (i.e. a subrange of input reads)
        global_data_t *globals;
        int rp_id; // ID of this read partition, in [ 0, num_cpu_threads ).
        pthread_t thread;
        int64_t rp_start_id, rp_end_id; // start and end IDs of this read partition (end is exclusive)
        int64_t *rp_bucket_sizes; // bucket sizes for this read partition; len =
        int64_t *rp_bucket_offsets;
        int64_t rp_lv1_differential_base; // the initial offset globals.lv1_items
    };

    struct bucketpartition_data_t {
        // local data for each bucket partition (i.e. a range of buckets), used in lv.2 (extract substring)
        global_data_t *globals;
        int bp_id;
        pthread_t thread;
        int bp_start_bucket, bp_end_bucket;
    };

    struct outputpartition_data_t {
        // output data for each thread
        global_data_t *globals;
        int op_id;
        pthread_t thread;
        int64_t op_start_index, op_end_index;
    };

    // param: must be set
    global_data_t *g_;
    int64_t num_items_;
    int num_cpu_threads_;
    int num_output_threads_;
    int64_t max_lv1_items_, max_lv2_items_;

    bool lv1_just_go_;
    int64_t max_mem_remain_;
    int64_t bytes_per_sorting_item_;
    std::vector<bool> cur_lv1_buckets_;

    // other data
    int64_t *bucket_sizes_;
    int *ori_bucket_id_;
    int *bucket_rank_;
    readpartition_data_t *rp_;
    bucketpartition_data_t *bp_;
    outputpartition_data_t *op_;

    // may change as cx1 goes
    int64_t lv1_num_items_, lv2_num_items_;
    int lv1_start_bucket_, lv1_end_bucket_;
    int lv2_start_bucket_, lv2_end_bucket_;
    std::vector<int64_t> lv1_items_special_;

    // === functions to specify a CX1 instance ===
    int64_t (*encode_lv1_diff_base_func_) (int64_t, global_data_t &);
    void (*prepare_func_) (global_data_t &); // num_items_, num_cpu_threads_ and num_output_threads_ must be set here
    void *(*lv0_calc_bucket_size_func_) (void *);
    void (*init_global_and_set_cx1_func_) (global_data_t &); // xxx set here
    void *(*lv1_fill_offset_func_) (void *);
    void (*lv1_sort_and_proc) (global_data_t &);
    void *(*lv2_extract_substr_func_) (void *);
    void (*lv2_sort_func_) (global_data_t &);
    void (*lv2_pre_output_partition_func_) (global_data_t &); // op_ set here
    void *(*lv2_output_func_) (void *);
    void (*lv2_post_output_func_) (global_data_t &);
    void (*post_proc_func_) (global_data_t &);

    CX1() : lv1_just_go_(false) {}

    // === single thread functions ===
    inline void adjust_mem(int64_t mem_avail, int64_t lv2_bytes_per_item, int64_t min_lv1_items, int64_t min_lv2_items) {
        // --- adjust max_lv2_items to fit memory ---
        while (max_lv2_items_ >= min_lv2_items) {
            int64_t mem_lv2 = lv2_bytes_per_item * max_lv2_items_;

            if (mem_avail <= mem_lv2) {
                max_lv2_items_ *= 0.95;
                continue;
            }

            max_lv1_items_ = (mem_avail - mem_lv2) / kLv1BytePerItem;

            if (max_lv1_items_ < min_lv1_items ||
                    max_lv1_items_ < max_lv2_items_) {
                max_lv2_items_ *= 0.95;
            }
            else {
                break;
            }
        }

        if (max_lv2_items_ < min_lv2_items) {
            xerr_and_exit("No enough memory to process CX1.\n");
        }

        // --- adjust max_lv2_items to fit more lv1 item ---
        // TODO: 4 is arbitrary chosen, not fine tune
        while (max_lv2_items_ * 4 > max_lv1_items_) {
            if (max_lv2_items_ * 0.95 >= min_lv2_items) {
                max_lv2_items_ *= 0.95;
                max_lv1_items_ = (mem_avail - lv2_bytes_per_item * max_lv2_items_) / kLv1BytePerItem;
            }
            else {
                break;
            }
        }
    }

    inline void adjust_mem_just_go(int64_t mem_avail, int64_t bytes_per_sorting_item, int64_t min_lv1_items, int64_t min_sorting_items,
                                   int64_t max_sorting_items, int64_t &max_lv1_items, int64_t &num_sorting_items) {
        num_sorting_items = max_sorting_items;

        while (num_sorting_items >= min_sorting_items) {
            int64_t mem_sorting_items = bytes_per_sorting_item * num_sorting_items;

            if (mem_avail < mem_sorting_items) {
                mem_sorting_items = num_sorting_items * 0.95;
                continue;
            }

            max_lv1_items = (mem_avail - mem_sorting_items) / kLv1BytePerItem;

            if (max_lv1_items < min_lv1_items || max_lv1_items < num_sorting_items) {
                num_sorting_items *= 0.95;
            }
            else {
                break;
            }
        }

        if (num_sorting_items < min_sorting_items) {
            if (num_sorting_items < min_sorting_items) {
                xerr_and_exit("No enough memory to process CX1.\n");
            }
        }

        // --- adjust num_sorting_items to fit more lv1 item ---
        // TODO: 4 is arbitrary chosen, not fine tune
        while (num_sorting_items * 4 > max_lv1_items) {
            if (num_sorting_items * 0.95 >= min_sorting_items) {
                num_sorting_items *= 0.95;
                max_lv1_items = (mem_avail - bytes_per_sorting_item * num_sorting_items) / kLv1BytePerItem;
            }
            else {
                break;
            }
        }
    }

    inline void prepare_rp_and_bp_() { // call after prepare_func_
        rp_ = (readpartition_data_t *) MallocAndCheck(sizeof(readpartition_data_t) * num_cpu_threads_, __FILE__, __LINE__);
        bp_ = (bucketpartition_data_t *) MallocAndCheck(sizeof(bucketpartition_data_t) * (num_cpu_threads_ - num_output_threads_), __FILE__, __LINE__);
        op_ = (outputpartition_data_t *) MallocAndCheck(sizeof(outputpartition_data_t) * num_output_threads_, __FILE__, __LINE__);

        for (int t = 0; t < num_cpu_threads_; ++t) {
            struct readpartition_data_t &rp = rp_[t];
            rp.rp_id = t;
            rp.globals = g_;
            rp.rp_bucket_sizes = (int64_t *) MallocAndCheck(kNumBuckets * sizeof(int64_t), __FILE__, __LINE__);
            rp.rp_bucket_offsets = (int64_t *) MallocAndCheck(kNumBuckets * sizeof(int64_t), __FILE__, __LINE__);
            // distribute reads to partitions
            int64_t average = num_items_ / num_cpu_threads_;
            rp.rp_start_id = t * average;
            rp.rp_end_id = t < num_cpu_threads_ - 1 ? (t + 1) * average : num_items_;
            rp.rp_lv1_differential_base = encode_lv1_diff_base_func_(rp.rp_start_id, *g_);
        }

        // init bucket partitions
        for (int t = 0; t < num_cpu_threads_ - num_output_threads_; ++t) {
            struct bucketpartition_data_t &bp = bp_[t];
            bp.bp_id = t;
            bp.globals = g_;
        }

        // init op
        for (int t = 0; t < num_output_threads_; ++t) {
            op_[t].op_id = t;
            op_[t].globals = g_;
        }

        ori_bucket_id_ = (int *) MallocAndCheck(sizeof(int) * kNumBuckets, __FILE__, __LINE__);
        bucket_rank_ = (int *) MallocAndCheck(sizeof(int) * kNumBuckets, __FILE__, __LINE__);
        bucket_sizes_ = (int64_t *) MallocAndCheck(kNumBuckets * sizeof(int64_t), __FILE__, __LINE__);

        for (int i = 0; i < kNumBuckets; ++i) {
            ori_bucket_id_[i] = i;
            bucket_rank_[i] = i;
        }
    }

    inline void clean_() {
        for (int t = 0; t < num_cpu_threads_; ++t) {
            free(rp_[t].rp_bucket_sizes);
            free(rp_[t].rp_bucket_offsets);
        }

        free(rp_);
        free(bp_);
        free(op_);
        free(bucket_sizes_);
        free(ori_bucket_id_);
        free(bucket_rank_);
    }

    inline int find_end_buckets_(int start_bucket, int end_limit, int64_t item_limit, int64_t &num_items) {
        cur_lv1_buckets_.resize(kNumBuckets);
        std::fill(cur_lv1_buckets_.begin(), cur_lv1_buckets_.end(), false);

        num_items = 0;
        int end_bucket = start_bucket;

        while (end_bucket < end_limit) { // simple linear scan
            if (num_items + bucket_sizes_[end_bucket] > item_limit) {
                return end_bucket;
            }

            num_items += bucket_sizes_[end_bucket];
            cur_lv1_buckets_[ori_bucket_id_[end_bucket]] = true;
            end_bucket++;
        }

        return end_limit;
    }

    inline void reorder_buckets_() {
        std::vector<std::pair<int64_t, int> > tmp_v(kNumBuckets);

        for (int i = 0; i < kNumBuckets; ++i) {
            tmp_v[i] = std::make_pair(bucket_sizes_[i], i);
        }

        std::sort(tmp_v.rbegin(), tmp_v.rend());

        for (int i = 0; i < kNumBuckets; ++i) {
            bucket_sizes_[i] = tmp_v[i].first;
            ori_bucket_id_[i] = tmp_v[i].second;
            bucket_rank_[tmp_v[i].second] = i;
        }

        for (int tid = 0; tid < num_cpu_threads_; ++tid) {
            std::vector<int64_t> old_rp_bucket_sizes(rp_[tid].rp_bucket_sizes, rp_[tid].rp_bucket_sizes + kNumBuckets);

            for (int i = 0; i < kNumBuckets; ++i) {
                rp_[tid].rp_bucket_sizes[i] = old_rp_bucket_sizes[tmp_v[i].second];
            }
        }
    }

    inline int find_end_buckets_with_rank_(int start_bucket, int end_limit, int64_t mem_limit, int bytes_per_sorting_items, int64_t &num_items) {
        num_items = 0;
        int end_bucket = start_bucket;
        int used_threads = 0;
        int64_t mem_sorting_items = 0;

        cur_lv1_buckets_.resize(kNumBuckets);
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
        int64_t *offsets = rp_[0].rp_bucket_offsets;
        offsets[lv1_start_bucket_] = 0;

        for (int b = lv1_start_bucket_ + 1; b < lv1_end_bucket_; ++b) {
            offsets[b] = offsets[b - 1] + bucket_sizes_[b - 1]; // accumulate
        }

        // then for each read partition
        for (int t = 1; t < num_cpu_threads_; ++t) {
            int64_t *this_offsets = rp_[t].rp_bucket_offsets;
            int64_t *prev_offsets = rp_[t - 1].rp_bucket_offsets;
            int64_t *sizes = rp_[t - 1].rp_bucket_sizes;

            for (int b = lv1_start_bucket_; b < lv1_end_bucket_; ++b) {
                this_offsets[b] = prev_offsets[b] + sizes[b];
            }
        }
    }

    inline void lv2_distribute_bucket_partitions_() {
        int64_t average = lv2_num_items_ / (num_cpu_threads_ - num_output_threads_);
        // recall: we only have (num_cpu_threads_ - num_output_threads_) bucketpartitions

        int bucket = lv2_start_bucket_;

        for (int t = 0; t < num_cpu_threads_ - num_output_threads_ - 1; ++t) {
            int64_t num_items = 0;
            bp_[t].bp_start_bucket = bucket;

            while (bucket < lv2_end_bucket_) {
                num_items += bucket_sizes_[bucket++];

                if (num_items >= average) {
                    break;
                }
            }

            bp_[t].bp_end_bucket = bucket;
        }

        // last
        bp_[num_cpu_threads_ - num_output_threads_ - 1].bp_start_bucket = bucket;
        bp_[num_cpu_threads_ - num_output_threads_ - 1].bp_end_bucket = lv2_end_bucket_;
    }

    // === multi-thread wrappers ====
    inline void lv0_calc_bucket_size_mt_() {
        for (int t = 0; t < num_cpu_threads_; ++t) {
            pthread_create(&(rp_[t].thread), NULL, lv0_calc_bucket_size_func_, &rp_[t]);
        }

        for (int t = 0; t < num_cpu_threads_; ++t) {
            pthread_join(rp_[t].thread, NULL);
        }

        // sum up readpartitions bucketsizes to form global bucketsizes
        memset(bucket_sizes_, 0, kNumBuckets * sizeof(bucket_sizes_[0]));

        // the array accesses in this loop are optimized by the compiler??
        for (int t = 0; t < num_cpu_threads_; ++t) {
            for (int b = 0; b < kNumBuckets; ++b) {
                bucket_sizes_[b] += rp_[t].rp_bucket_sizes[b];
            }
        }
    }

    inline void lv1_fill_offset_mt_() {
        lv1_items_special_.clear();
        lv1_compute_offset_();

        // create threads
        for (int t = 0; t < num_cpu_threads_; ++t) {
            pthread_create(&(rp_[t].thread), NULL, lv1_fill_offset_func_, &rp_[t]);
        }

        for (int t = 0; t < num_cpu_threads_; ++t) {
            pthread_join(rp_[t].thread, NULL);
        }

        // revert rp_bucket_offsets
        lv1_compute_offset_();
    }

    inline void lv2_extract_substr_mt_() {
        lv2_distribute_bucket_partitions_();

        // create threads
        for (int t = 0; t < num_cpu_threads_ - num_output_threads_; ++t) {
            pthread_create(&(bp_[t].thread), NULL, lv2_extract_substr_func_, &bp_[t]);
        }

        for (int t = 0; t < num_cpu_threads_ - num_output_threads_; ++t) {
            pthread_join(bp_[t].thread, NULL);
        }
    }

    inline void lv2_output_mt_() {
        for (int t = 0; t < num_output_threads_; ++t) {
            op_[t].op_id = t;
            op_[t].globals = g_;
            pthread_create(&(op_[t].thread), NULL, lv2_output_func_, &op_[t]);
        }
    }

    inline void lv2_output_join_() {
        for (int t = 0; t < num_output_threads_; ++t) {
            pthread_join(op_[t].thread, NULL);
        }
    }

    // === go go go ===
    inline void run() {
        xtimer_t lv0_timer;

        // read input & prepare
        if (kCX1Verbose >= 2) {
            lv0_timer.reset();
            lv0_timer.start();
            xlog("Preparing data...\n");
        }

        prepare_func_(*g_);

        if (kCX1Verbose >= 2) {
            lv0_timer.stop();
            xlog("Preparing data... Done. Time elapsed: %.4f\n", lv0_timer.elapsed());
        }


        if (kCX1Verbose >= 2) {
            lv0_timer.reset();
            lv0_timer.start();
            xlog("Preparing partitions and initialing global data...\n");
        }

        // prepare rp bp and op
        prepare_rp_and_bp_();
        // calc bucket size
        lv0_calc_bucket_size_mt_();
        // init global datas
        init_global_and_set_cx1_func_(*g_);

        if (lv1_just_go_) {
            reorder_buckets_();
        }

        if (kCX1Verbose >= 2) {
            lv0_timer.stop();
            xlog("Preparing partitions and initialing global data... Done. Time elapsed: %.4f\n", lv0_timer.elapsed());
        }

        if (kCX1Verbose >= 2) {
            lv0_timer.reset();
            lv0_timer.start();
            xlog("Start main loop...\n");
        }

        // === start main loop ===
        bool output_thread_created = false;
        int lv1_iteration = 0;
        lv1_start_bucket_ = 0;

        while (lv1_start_bucket_ < kNumBuckets) {
            xtimer_t lv1_timer;

            lv1_iteration++;

            // --- finds the bucket range for this iteration ---
            if (lv1_just_go_) {
                lv1_end_bucket_ = find_end_buckets_with_rank_(lv1_start_bucket_, kNumBuckets, max_mem_remain_, bytes_per_sorting_item_, lv1_num_items_);
            }
            else {
                lv1_end_bucket_ = find_end_buckets_(lv1_start_bucket_, kNumBuckets, max_lv1_items_, lv1_num_items_);
            }

            if (lv1_num_items_ == 0) {
                fprintf(stderr, "Bucket %d too large for lv1: %lld > %lld\n", lv1_end_bucket_, (long long)bucket_sizes_[lv1_end_bucket_], (long long)max_lv1_items_);
                exit(1);
            }

            if (kCX1Verbose >= 3) {
                lv1_timer.reset();
                lv1_timer.start();
                xlog("Lv1 scanning from bucket %d to %d\n", lv1_start_bucket_, lv1_end_bucket_);
            }

            // --- scan to fill offset ---
            lv1_fill_offset_mt_();

            if (lv1_items_special_.size() > kSpDiffMaxNum) {
                fprintf(stderr, "Too many large diff items (%lu) from in buckets [%d, %d)\n", lv1_items_special_.size(), lv1_start_bucket_, lv1_end_bucket_);
                exit(1);
            }

            if (kCX1Verbose >= 3) {
                lv1_timer.stop();
                xlog("Lv1 scanning done. Large diff: %lu. Time elapsed: %.4f\n", lv1_items_special_.size(), lv1_timer.elapsed());
                lv1_timer.reset();
                lv1_timer.start();
            }

            if (lv1_just_go_) {
                lv1_sort_and_proc(*g_);
            }
            else {
                // --- lv2 loop ---
                int lv2_iteration = 0;
                lv2_start_bucket_ = lv1_start_bucket_;

                while (lv2_start_bucket_ < lv1_end_bucket_) {
                    xtimer_t lv2_timer;

                    lv2_iteration++;
                    lv2_end_bucket_ = find_end_buckets_(lv2_start_bucket_, lv1_end_bucket_, max_lv2_items_, lv2_num_items_);

                    if (lv2_num_items_ == 0) {
                        fprintf(stderr, "Bucket %d too large for lv2: %lld > %lld\n", lv2_end_bucket_, (long long)bucket_sizes_[lv2_end_bucket_], (long long)max_lv2_items_);
                        exit(1);
                    }

                    if (kCX1Verbose >= 4) {
                        lv2_timer.reset();
                        lv2_timer.start();
                        xlog("Lv2 fetching substrings from bucket %d to %d\n", lv2_start_bucket_, lv2_end_bucket_);
                    }

                    // --- extract lv2 substr and sort ---
                    lv2_extract_substr_mt_();

                    if (kCX1Verbose >= 4) {
                        lv2_timer.stop();
                        xlog("Lv2 fetching substrings done. Time elapsed: %.4f\n", lv2_timer.elapsed());
                        lv2_timer.reset();
                        lv2_timer.start();
                    }

                    lv2_sort_func_(*g_);

                    if (kCX1Verbose >= 4) {
                        lv2_timer.stop();
                        xlog("Lv2 sorting done. Time elapsed: %.4f\n", lv2_timer.elapsed());
                        lv2_timer.reset();
                        lv2_timer.start();
                    }

                    // --- the output is pipelined, join the previous one ---
                    if (output_thread_created) {
                        lv2_output_join_();
                        lv2_post_output_func_(*g_);
                    }

                    // --- the create new output threads ---
                    lv2_pre_output_partition_func_(*g_);
                    lv2_output_mt_();
                    output_thread_created = true;

                    lv2_start_bucket_ = lv2_end_bucket_;
                }
            }

            if (kCX1Verbose >= 3) {
                lv1_timer.stop();
                xlog("Lv1 fetching & sorting done. Time elapsed: %.4f\n", lv1_timer.elapsed());
            }

            lv1_start_bucket_ = lv1_end_bucket_;
        }

        if (output_thread_created) {
            lv2_output_join_();
            lv2_post_output_func_(*g_);
        }

        if (kCX1Verbose >= 2) {
            lv0_timer.stop();
            xlog("Main loop done. Time elapsed: %.4f\n", lv0_timer.elapsed());
        }

        if (kCX1Verbose >= 2) {
            lv0_timer.reset();
            lv0_timer.start();
            xlog("Postprocessing...\n");
        }

        post_proc_func_(*g_);
        clean_();

        if (kCX1Verbose >= 2) {
            lv0_timer.stop();
            xlog("Postprocess done. Time elapsed: %.4f\n", lv0_timer.elapsed());
        }
    }
};

#endif // CX1_H__