#ifndef EDGE_IO_H__
#define EDGE_IO_H__

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include <string>
#include <vector>

#include "definitions.h"
#include "mem_file_checker-inl.h"
#include "utils.h"

struct PartitionRecord {
	int thread_id;
	long long starting_offset;
	long long total_number;

	PartitionRecord(): thread_id(-1), starting_offset(0), total_number(0) {}		
};

class EdgeWriter {
  private:
  	int32_t kmer_size_;
  	int32_t words_per_edge_;
  	int32_t num_threads_;
  	int32_t num_buckets_;

  	bool unsorted_;
  	int64_t num_unsorted_edges_;

  	std::string file_prefix_;
  	std::vector<FILE*> files_;
  	std::vector<int32_t> cur_bucket_;
  	std::vector<int64_t> cur_num_edges_;
  	std::vector<PartitionRecord> p_rec_;

  	bool is_opened_;

  public:
  	EdgeWriter(): unsorted_(false), is_opened_(false) {};
  	~EdgeWriter() { destroy(); }

  	void set_kmer_size(int32_t k) {
  		kmer_size_ = k;
  		words_per_edge_ = DivCeiling((k + 1) * 2 + 16, 32);
  	}

  	void set_num_threads(int32_t num_threads) {
  		num_threads_ = num_threads;
  	}

  	void set_file_prefix(const std::string &prefix) {
  		file_prefix_ = prefix;
  	}

  	void set_num_buckets(int num_buckets) {
  		num_buckets_ = num_buckets;
  	}

  	void set_unsorted() {
  		num_buckets_ = 0;
  		p_rec_.clear();
  		unsorted_ = true;
  		num_unsorted_edges_ = 0;
  	}

  	void init_files() {
  		assert(!is_opened_);

  		files_.resize(num_threads_);
  		cur_bucket_.resize(num_threads_, -1);
  		cur_num_edges_.resize(num_threads_, 0);
  		p_rec_.resize(num_buckets_, PartitionRecord());

  		for (int i = 0; i < num_threads_; ++i) {
  			files_[i] = OpenFileAndCheck(FormatString("%s.edges.%d", file_prefix_.c_str(), i), "wb");
  		}

  		is_opened_ = true;
  	}

  	void write(uint32_t *edge_ptr, int32_t bucket, int tid) {
  		// assert(bucket >= cur_bucket_[tid]);
  		if (bucket != cur_bucket_[tid]) {
  			assert(p_rec_[bucket].thread_id == -1);
  			p_rec_[bucket].thread_id = tid;
  			p_rec_[bucket].starting_offset = cur_num_edges_[tid];
  			cur_bucket_[tid] = bucket;
  		}

  		fwrite(edge_ptr, sizeof(uint32_t), words_per_edge_, files_[tid]);
  		++cur_num_edges_[tid];
  		++p_rec_[bucket].total_number;
  	}

  	void write_unsorted(uint32_t *edge_ptr, int tid) {
  		fwrite(edge_ptr, sizeof(uint32_t), words_per_edge_, files_[tid]);
  		++num_unsorted_edges_;
  	}

  	void destroy() {
  		if (is_opened_) {
  			for (int i = 0; i < num_threads_; ++i) {
	  			fclose(files_[i]);
	  		}

	  		int64_t num_edges = 0;
	  		if (!unsorted_) {
		  		for (unsigned i = 0; i < p_rec_.size(); ++i) {
		  			num_edges += p_rec_[i].total_number;
		  		}
	  		} else {
	  			num_edges = num_unsorted_edges_;
	  		}

	  		FILE *info = OpenFileAndCheck(FormatString("%s.edges.info", file_prefix_.c_str()), "w");
	  		fprintf(info, "kmer_size %d\n", (int)kmer_size_);
	  		fprintf(info, "words_per_edge %d\n", (int)words_per_edge_);
	  		fprintf(info, "num_threads %d\n", (int)num_threads_);
	  		fprintf(info, "num_bucket %d\n", (int)num_buckets_);
	  		fprintf(info, "num_edges %lld\n", (long long)num_edges);
	  		for (unsigned i = 0; i < p_rec_.size(); ++i) {
	  			fprintf(info, "%u %d %lld %lld\n", i, p_rec_[i].thread_id, (long long)p_rec_[i].starting_offset, (long long)p_rec_[i].total_number);
	  		}
	  		fclose(info);

			files_.clear();
	  		cur_bucket_.clear();
	  		cur_num_edges_.clear();
	  		p_rec_.clear();

	  		is_opened_ = false;
  		}
  	}
};

class EdgeReader {
  private:
  	int kmer_size_;
  	int words_per_edge_;
  	int num_files_;
  	int num_buckets_;
  	long long num_edges_;

  	std::string file_prefix_;
    std::vector<FILE*> files_;
  	std::vector<PartitionRecord> p_rec_;
  	std::vector<uint32_t> buf_;
  	int buf_idx_;
  	int buf_size_;

  	int cur_bucket_;
  	int cur_file_num_; // used for unsorted edges
  	long long cur_bucket_cnt_;
  	FILE *cur_file_;

    bool is_opened_;

  	static const int kEdgesInBuf = 4096;

  public:
  	EdgeReader(): is_opened_(false) {
  		cur_file_num_ = -1;
  		cur_bucket_ = -1;
  		cur_bucket_cnt_ = 0;
  		cur_file_ = NULL;
  	}
  	~EdgeReader() {
  		destroy();
  	}

  	void set_file_prefix(const std::string &prefix) {
  		file_prefix_ = prefix;
  	}

  	void read_info() {
  		FILE *info = OpenFileAndCheck(FormatString("%s.edges.info", file_prefix_.c_str()), "r");
  		assert(fscanf(info, "kmer_size %d\n", &kmer_size_) == 1);
  		assert(fscanf(info, "words_per_edge %d\n", &words_per_edge_) == 1);
  		assert(fscanf(info, "num_threads %d\n", &num_files_) == 1);
  		assert(fscanf(info, "num_bucket %d\n", &num_buckets_) == 1);
  		assert(fscanf(info, "num_edges %lld\n", &num_edges_) == 1);
  		p_rec_.resize(num_buckets_);
  		for (int i = 0; i < num_buckets_; ++i) {
  			unsigned dummy;
  			assert(fscanf(info, "%u %d %lld %lld\n", &dummy, &p_rec_[i].thread_id, &p_rec_[i].starting_offset, &p_rec_[i].total_number) == 4);
  		}

  		buf_.resize(words_per_edge_ * kEdgesInBuf);
  		buf_size_ = 0;
  		buf_idx_ = 0;
  		fclose(info);
  	}

    void init_files() {
      assert(!is_opened_);
      files_.resize(num_files_);
      for (int i = 0; i < num_files_; ++i) {
        files_[i] = OpenFileAndCheck(FormatString("%s.edges.%d", file_prefix_.c_str(), i), "rb");
      }
      is_opened_ = true;
    }

  	bool is_unsorted() {
  		return num_buckets_ == 0;
  	}

  	int kmer_size() { return kmer_size_; }
  	int words_per_edge() { return words_per_edge_; }
  	int num_edges() { return num_edges_; }

  	uint32_t *NextSortedEdge() {
  		if (cur_bucket_ >= num_buckets_) { return NULL; }

  		while (cur_bucket_ == -1 || cur_bucket_cnt_ >= p_rec_[cur_bucket_].total_number) {

  			buf_idx_ = 0;
  			buf_size_ = 0;

  			++cur_bucket_;
  			while (cur_bucket_ < num_buckets_ && p_rec_[cur_bucket_].thread_id < 0) {
  				++cur_bucket_;
  			}

  			if (cur_bucket_ >= num_buckets_) { return NULL; }
  			cur_bucket_cnt_ = 0;

  			cur_file_ = files_[p_rec_[cur_bucket_].thread_id];
  			fseek(cur_file_, sizeof(uint32_t) * words_per_edge_ * p_rec_[cur_bucket_].starting_offset, SEEK_SET);
  			assert(!ferror(cur_file_));
  		}

  		if (buf_idx_ == buf_size_) {
  			buf_size_ = kEdgesInBuf;
  			if (buf_size_ > p_rec_[cur_bucket_].total_number - cur_bucket_cnt_) {
  				buf_size_ = p_rec_[cur_bucket_].total_number - cur_bucket_cnt_;
  			}
  			buf_size_ *= words_per_edge_;
  			assert(buf_size_ != 0);

  			assert(fread(&buf_[0], sizeof(uint32_t), buf_size_, cur_file_) == (unsigned)buf_size_);
  			buf_idx_ = 0;
  		}

  		++cur_bucket_cnt_;
  		buf_idx_ += words_per_edge_;
  		return &buf_[buf_idx_ - words_per_edge_];
  	}

  	uint32_t *NextUnsortedEdge() {
  		if (cur_file_num_ < 0) {
  			cur_file_num_ = 0;
  		}

  		if (cur_file_num_ >= num_files_) {
  			return NULL;
  		}

  		while (buf_idx_ >= buf_size_) {
  			buf_idx_ = 0;
  			buf_size_ = fread(&buf_[0], sizeof(uint32_t), kEdgesInBuf * words_per_edge_, cur_file_);

  			if (buf_size_ != 0) {
  				break;
  			}

  			cur_file_ = NULL;
  			cur_file_num_++;
  			if (cur_file_num_ >= num_files_) {
  				return NULL;
  			}

  			cur_file_ = files_[cur_file_num_];
  		}

  		buf_idx_ += words_per_edge_;
  		return &buf_[buf_idx_ - words_per_edge_];
  	}

  	void destroy() {
      if (is_opened_) {
        for (int i = 0; i < num_files_; ++i) {
          fclose(files_[i]);
        }
        is_opened_ = false;
        cur_file_ = NULL;
      }
  	}
};

#endif