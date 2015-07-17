#ifndef SDBG_MULTI_IO__
#define SDBG_MULTI_IO__

#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <string>
#include <vector>
#include <algorithm>

#include "definitions.h"
#include "utils.h"
#include "mem_file_checker-inl.h"

struct SdbgPartitionRecord {
	int thread_id;
	int64_t starting_offset;
	int64_t num_items;
	int64_t num_tips;
	int64_t num_large_mul;
	int64_t num_last1;
	int64_t num_w[9];

	SdbgPartitionRecord(): thread_id(-1), starting_offset(0), num_items(0),
	    num_tips(0), num_large_mul(0), num_last1(0) {
	    memset(num_w, 0, sizeof(num_w));
	}
};

class SdbgWriter {
  private:
  	std::string file_prefix_;
  	int num_threads_;
  	int num_buckets_;

  	std::vector<FILE*> files_;
  	std::vector<int> cur_bucket_;
  	std::vector<int64_t> cur_thread_offset_;	// offset in BYTE
  	std::vector<SdbgPartitionRecord> p_rec_;

  	bool is_opened_;
  	int kmer_size_;
  	int words_per_tip_label_;

  public:

	SdbgWriter() {}
	~SdbgWriter() {
		destroy();
	}

	void set_num_threads(int num_threads) { num_threads_ = num_threads; }
	void set_file_prefix(const std::string &file_prefix) { file_prefix_ = file_prefix; }
	void set_kmer_size(int k) {
		kmer_size_ = k;
		words_per_tip_label_ = DivCeiling(k * 2, 32);
	}
	void set_num_buckets(int num_buckets) { num_buckets_ = num_buckets; }

	void init_files() {
		files_.resize(num_threads_);
		cur_bucket_.resize(num_threads_, -1);
		cur_thread_offset_.resize(num_threads_, 0);
		p_rec_.resize(num_buckets_);

		for (int i = 0; i < num_threads_; ++i) {
			files_[i] = OpenFileAndCheck(FormatString("%s.sdbg.%d", file_prefix_.c_str(), i), "wb");
		}

		is_opened_ = true;
	}

	void write(int tid, int32_t bucket, int w, int last, int tip, multi_t multiplicity, uint32_t *packed_tip_label) {
		assert(tid < num_threads_);

		if (bucket != cur_bucket_[tid]) {
			cur_bucket_[tid] = bucket;
			assert(p_rec_[bucket].thread_id == -1);
			p_rec_[bucket].thread_id = tid;
			p_rec_[bucket].starting_offset = cur_thread_offset_[tid];
		}

		uint16_t packed_sdbg_item = w | (last << 4) | (tip << 5) | (std::min(multiplicity, (multi_t)kMulti2Sp) << 8);
		fwrite(&packed_sdbg_item, sizeof(uint16_t), 1, files_[tid]);
		++p_rec_[bucket].num_items;
		++p_rec_[bucket].num_w[w];
		p_rec_[bucket].num_last1 += last;
		cur_thread_offset_[tid] += sizeof(uint16_t);

		if (multiplicity > kMaxMulti2_t) {
			fwrite(&multiplicity, sizeof(multi_t), 1, files_[tid]);
			multiplicity = kMulti2Sp;
			++p_rec_[bucket].num_large_mul;
			cur_thread_offset_[tid] += sizeof(multi_t);
		}

		if (tip) {
			fwrite(&packed_tip_label, sizeof(uint32_t), words_per_tip_label_, files_[tid]);
			++p_rec_[bucket].num_tips;
			cur_thread_offset_[tid] += sizeof(uint32_t) * words_per_tip_label_;
		}
	}

	int64_t num_edges() {
		int64_t total_edges = 0;
		for (int i = 0; i < num_buckets_; ++i) {
			total_edges += p_rec_[i].num_items;
		}
		return total_edges;
	}

	int64_t num_w(int w) {
		int64_t ret = 0;
		for (int i = 0; i < num_buckets_; ++i) {
			ret += p_rec_[i].num_w[w];
		}
		return ret;
	}

	int64_t num_last1() {
		int64_t ret = 0;
		for (int i = 0; i < num_buckets_; ++i) {
			ret += p_rec_[i].num_last1;
		}
		return ret;
	}

	int64_t num_tips() {
		int64_t ret = 0;
		for (int i = 0; i < num_buckets_; ++i) {
			ret += p_rec_[i].num_tips;
		}
		return ret;
	}

	void destroy() {
		if (is_opened_) {
			for (int i = 0; i < num_threads_; ++i) {
				fclose(files_[i]);
			}

			FILE *sdbg_info = OpenFileAndCheck((file_prefix_ + ".sdbg_info").c_str(), "w");
			fprintf(sdbg_info, "k %d\n", kmer_size_);
			fprintf(sdbg_info, "words_per_tip_label %d\n", words_per_tip_label_);

			int64_t total_edges = 0;
			int64_t total_tips = 0;
			int64_t total_large_mul = 0;

			for (int i = 0; i < num_buckets_; ++i) {
				total_edges += p_rec_[i].num_items;
				total_tips += p_rec_[i].num_tips;
				total_large_mul += p_rec_[i].num_large_mul;
			}

			fprintf(sdbg_info, "total_size %lld\n", (long long)total_edges);
			fprintf(sdbg_info, "num_tips %lld\n", (long long)total_tips);
			fprintf(sdbg_info, "large_multi %lld\n", (long long)total_large_mul);

			for (int i = 0; i < num_buckets_; ++i) {
				fprintf(sdbg_info, "%d %d %lld %lld %lld %lld\n", i,
					p_rec_[i].thread_id,
					(long long)p_rec_[i].starting_offset,
					(long long)p_rec_[i].num_items, 
					(long long)p_rec_[i].num_tips, 
					(long long)p_rec_[i].num_large_mul);
			}

			files_.clear();
		  	cur_bucket_.clear();
		  	cur_thread_offset_.clear();	// offset in BYTE
		  	p_rec_.clear();

			is_opened_ = false;
		}
	}
};

#endif