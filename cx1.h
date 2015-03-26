#ifndef CX1_H__
#define CX1_H__

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <pthread.h>
#include <vector>
#include <algorithm>

#include "mem_file_checker-inl.h"

 /**
  * @brief    an CX1 engine
  * @details  use CX1 algorithm to do all kinds of things related to substring sorting
  * 
  * @tparam global_data_type the type of global datas used in a specified CX1 engine
 *                           must contain a member `CX1* cx1`, where cx1->g_ is itself,
 *                           that enables interactions with CX1 functions
  * @tparam kNumBuckets      number of buckets 
  */
template <typename global_data_type, int kNumBuckets>
struct CX1
{
	typedef global_data_t global_data_type;
	static const int kCX1Verbose = 3; // tunable
	// other settings, don't change
	static const int kGPUBytePerItem = 16; // key & value, 4 byte each. double for radix sort internal buffer
	static const int kLv1BytePerItem = 4; // 32-bit differatial offset
  
	struct readpartition_data_t 
	{  // local data for each read partition (i.e. a subrange of input reads)
	    struct global_data_t* globals;
	    int rp_id; // ID of this read partition, in [ 0, num_cpu_threads ).
	    pthread_t thread;
	    int64_t rp_start_id, rp_end_id; // start and end IDs of this read partition (end is exclusive)
	    int64_t* rp_bucket_sizes; // bucket sizes for this read partition; len = 
	    int64_t* rp_bucket_offsets;
	    int64_t rp_lv1_differential_base; // the initial offset globals.lv1_items
	};

	struct bucketpartition_data_t 
	{ // local data for each bucket partition (i.e. a range of buckets), used in lv.2 (extract substring)
	    struct global_data_t* globals;
	    int bp_id;
	    pthread_t thread;
	    int bp_start_bucket, bp_end_bucket;
	};

	struct outputpartition_data_t
	{ // output data for each thread
	    struct global_data_t *globals;
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

	// other data
	int64_t *bucket_sizes_;
	readpartition_data_t *rp_;
	bucketpartition_data_t *bp_;
	outputpartition_data_t *op_;

	// may change as cx1 goes
	int64_t lv1_num_items_, lv2_num_items_;
	int64_t lv1_start_bucket_, lv1_end_bucket_;
	int64_t lv2_start_bucket_, lv2_end_bucket_;
	std::vector<int64_t> lv1_items_special_;

	// === functions to specify a CX1 instance ===
	int64_t (*encode_lv1_diff_base_func_) (int64_t, global_data_t &);
	void (*prepare_func_) (global_data_t &, CX1 &); // num_items_, num_cpu_threads_ and num_output_threads_ must be set here
	(void*) (*lv0_calc_bucket_size_func_) (void*);
	void (*init_global_and_set_cx1_func_) (global_data_t &, CX1 &);
	(void*) (*lv1_fill_offset_func_) (void*);
	(void*) (*lv2_extract_substr_func_) (void*);
	void (*lv1_sort_func_) (global_data_t &);
	void (*lv2_pre_output_partition_func_) (global_data_t &, CX1 &);
	(void*) (*lv2_output_func_) (void*);
	void (*lv2_post_output_func_) (global_data_t &);
	void (*post_proc_func_) (global_data_t &);

	// === single thread functions ===	
	inline void adjust_mem(int64_t mem_avail, int64_t lv2_bytes_per_item, int min_lv1_items, int min_lv2_items) {
	    // --- adjust max_lv2_items to fit memory ---
	    while (max_lv2_items_ >= min_lv2_items) {
	        int64_t mem_lv2 = lv2_bytes_per_item * max_lv2_items_;
	        if (mem_remained <= mem_lv2) {
	            max_lv2_items_ *= 0.9;
	            continue;
	        }

	        max_lv1_items_ = (mem_remained - mem_lv2) / LV1_BYTES_PER_ITEM;
	        if (max_lv1_items_ < min_lv1_items || 
	            max_lv1_items_ < max_lv2_items_) {
	            max_lv2_items_ *= 0.9;
	        } else {
	            break;
	        }
	    }

	    if (max_lv2_items_ < min_lv2_items) {
	        fprintf(stderr, "[CX1] No enough memory to process.\n");
	        exit(1);
	    }

	    // --- adjust max_lv2_items to fit more lv1 item ---
	    // TODO: 4 is arbitrary chosen, not fine tune
	    while (max_lv2_items_ * 4 > max_lv1_items_) {
	        if (max_lv2_items_ * 0.95 >= min_lv2_items) {
	            max_lv2_items_ *= 0.95;
	            max_lv1_items_ = (mem_remained - lv2_bytes_per_item * max_lv2_items_) / LV1_BYTES_PER_ITEM;
	        } else {
	            break;
	        }
	    } 
	}

	inline void prepare_rp_and_bp_() { // call after prepare_func_
		rp_ = (readpartition_data_t*) MallocAndCheck(sizeof(readpartition_data_t) * num_cpu_threads_, __FILE__, __LINE__);
		bp_ = (bucketpartition_data_t*) MallocAndCheck(sizeof(bucketpartition_data_t) * (num_cpu_threads_ - num_output_threads_), __FILE__, __LINE__);
		op_ = (outputpartition_data_t*) MallocAndCheck(sizeof(outputpartition_data_t) * num_output_threads_, __FILE__, __LINE__);

		for (int t = 0; t < num_cpu_threads_; ++t) {
	        struct readpartition_data_t &rp = rp_[t];
	        rp.rp_id = t;
	        rp.globals = &g_;
	        rp.rp_bucket_sizes = (int64_t *) MallocAndCheck(kNumBuckets * sizeof(int64_t), __FILE__, __LINE__);
	        rp.rp_bucket_offsets = (int64_t *) MallocAndCheck(kNumBuckets * sizeof(int64_t), __FILE__, __LINE__);
	        // distribute reads to partitions
	        int64_t average = num_items_ / num_cpu_threads_;
	        rp.rp_start_id = t * average;
	        rp.rp_end_id = t < num_cpu_threads_ - 1 ? (t + 1) * average : num_items_;
	        rp.rp_lv1_differential_base = encode_lv1_diff_base_func_(rp.rp_start_id, g_);
	    }

	    // init bucket partitions
	    for (int t = 0; t < num_cpu_threads_ - num_output_threads_; ++t) {
	        struct bucketpartition_data_t &bp = bp_[t];
	        bp.bp_id = t;
	        bp.g_ = &g_;
	    }

	    bucket_sizes_ = (int64_t *) MallocAndCheck(kNumBuckets * sizeof(int64_t), __FILE__, __LINE__);
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
	}

	inline int find_end_buckets_(int start_bucket, int end_limit, int64_t item_limit, int64_t &num_items) {    
	    num_items = 0;
	    int end_bucket = start_bucket;
	    while (end_bucket < end_limit) { // simple linear scan
	        if (num_items + bucket_sizes_[end_bucket] > item_limit) {
	            return end_bucket;
	        }
	        num_items += bucket_sizes_[end_bucket];
	        end_bucket++;
	    }
	    return end_limit;
	}

	inline void lv1_compute_offset_() {
	    // compute "global" (thread 0) offsets first
	    int64_t *offsets = rp_[0].rp_bucket_offsets;
	    offsets[lv1_start_bucket_] = 0;
	    for (int b = lv1_start_bucket_ + 1; b < lv1_end_bucket_; ++b) {
	        offsets[b] = offsets[b-1] + bucket_sizes_[b-1]; // accumulate
	    }
	    // then for each read partition
	    for (int t = 1; t < num_cpu_threads_; ++t) {
	        int64_t *this_offsets = rp_[t].rp_bucket_offsets;
	        int64_t *prev_offsets = rp_[t-1].rp_bucket_offsets;
	        int64_t *sizes = rp_[t-1].rp_bucket_sizes;
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
	    bp_[num_cpu_threads_ -num_output_threads_ -1].bp_start_bucket = bucket;
	    bp_[num_cpu_threads_ -num_output_threads_ -1].bp_end_bucket = lv2_end_bucket_;
	}

	// === multi-thread wrappers ====
	inline void lv0_calc_bucket_size_mt_() {
		for (int t = 1; t < num_cpu_threads_; ++t) {
        	pthread_create(&(rp_[t].thread), NULL, lv0_calc_bucket_size_func_, &rp_[t]);
    	}
    	lv0_calc_bucket_size_func_(&rp_[0]);
    	for (int t = 1; t < num_cpu_threads_; ++t) {
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
	    for (int t = 1; t < num_cpu_threads_; ++t) {
	        pthread_create(&(rp_[t].thread), NULL, lv1_fill_offset_func_, &rp_[t]);
	    }
	    lv1_fill_offset_func_(&rp_[0]);
	    for (int t = 1; t < num_cpu_threads_; ++t) {
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
	        op_[t].globals = &g_;
	        pthread_create(op_[t].thread, NULL, Lv2CountingThread, &op_[t]);
	    }
	}

	inline void lv2_output_join_() {
		for (int i = 0; t < num_output_threads_; ++t) {
			pthread_join(op_[t].thread, NULL);
		}
	}

	// === go go go ===
	inline void run() {
		// read input & prepare
		prepare_func_(g_, *this);
		// prepare rp bp and op
		prepare_rp_and_bp_();
		// calc bucket size
		lv0_calc_bucket_size_func_();
		// init global datas
		init_global_and_set_cx1_func_(g_, *this);

		// === start main loop ===
		bool output_thread_created = false;
		int lv1_iteration = 0;
		lv1_start_bucket_ = 0;
		while (lv1_start_bucket_ < kNumBuckets) {
        	lv1_iteration++;
			// --- finds the bucket range for this iteration ---
			lv1_end_bucket_ = find_end_buckets_(bucket_sizes_, lv1_start_bucket_, kNumBuckets, max_lv1_items, lv1_num_items_);
			if (lv1_num_items_ == 0) {
				fprintf(stderr, "[CX1] Bucket %d too large for lv1: %lld > %lld\n", lv1_end_bucket_, bucket_sizes_[lv1_end_bucket_], max_lv1_items_);
				exit(1);
			}

			// --- scan to fill offset ---
			lv1_fill_offset_mt_();
			if (lv1_items_special_.size() > kSpDiffLimit) {
				fprintf(stderr, "[CX1] Too many large diff items (%lu) from in buckets [%d, %d)\n", lv1_item_special_.size(), lv1_start_bucket_, lv1_end_bucket_);
				exit(1);
			}

			// --- lv2 loop ---
			int lv2_iteration = 0;
			lv2_start_bucket_ = lv1_start_bucket_;
			while (lv2_start_bucket_ < lv1_end_bucket_) {
				lv2_iteration++;
				lv2_end_bucket_ = find_end_buckets_(bucket_sizes_, lv2_start_bucket_, lv1_end_bucket_, max_lv2_items_, lv2_num_items_);
				if (lv2_num_items_ == 0) {
					fprintf(stderr, "[CX1] Bucket %d too large for lv2: %lld > %lld\n", lv2_end_bucket_, bucket_sizes_[lv2_end_bucket_], max_lv2_items_);
					exit(1);
				}

				// --- extract lv2 substr and sort ---
				lv2_extract_substr_mt_();
				lv2_sort_func_(g_);

				// --- the output is pipelined, join the previous one ---
				if (output_thread_created) {
					lv2_output_join_();
					lv2_post_output_func_(g_);
				}

				// --- the create new output threads ---
				lv2_pre_output_partition_func_(g_, *this);
				lv2_output_mt_();
				output_thread_created = true;
			}
		}

		if (output_thread_created) {
			lv2_output_join_();
			lv2_post_output_func_(g_);
		}

		post_proc_func_(g_);
		clean_();
	}
};

#endif // CX1_H__