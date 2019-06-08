#ifndef MEGAHIT_EDGE_IO_H
#define MEGAHIT_EDGE_IO_H

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <utils/buffered_reader.h>

#include "definitions.h"
#include "utils/utils.h"

struct PartitionRecord {
  int thread_id;
  int64_t starting_offset;
  int64_t total_number;

  PartitionRecord() : thread_id(-1), starting_offset(0), total_number(0) {}
};

class EdgeWriter {
 private:
  uint32_t kmer_size_{};
  uint32_t words_per_edge_{};
  uint32_t num_threads_{};
  uint32_t num_buckets_{};

  bool unsorted_;
  std::vector<int64_t> num_unsorted_edges_;

  std::string file_prefix_;
  std::vector<std::unique_ptr<std::ofstream>> files_;
  std::vector<int32_t> cur_bucket_;
  std::vector<int64_t> cur_num_edges_;
  std::vector<PartitionRecord> p_rec_;

  bool is_opened_;

 public:
  class Snapshot {
   private:
    PartitionRecord p_rec;
    int bucket_id{-1};
    friend class EdgeWriter;
  };

  EdgeWriter() : unsorted_(false), is_opened_(false) {};
  ~EdgeWriter() { Finalize(); }

  void SetKmerSize(uint32_t k) {
    kmer_size_ = k;
    words_per_edge_ = DivCeiling((k + 1) * 2 + 16, 32);
  }

  void SetNumThreads(int32_t num_threads) { num_threads_ = num_threads; }
  void SetFilePrefix(const std::string &prefix) { file_prefix_ = prefix; }
  void SetNumBuckets(int num_buckets) { num_buckets_ = num_buckets; }
  void SetUnsorted() {
    num_buckets_ = 0;
    p_rec_.clear();
    unsorted_ = true;
    num_unsorted_edges_.clear();
    num_unsorted_edges_.resize(num_threads_, 0);
  }

  void InitFiles() {
    assert(!is_opened_);

    cur_bucket_.resize(num_threads_, -1);
    cur_num_edges_.resize(num_threads_, 0);
    p_rec_.resize(num_buckets_, PartitionRecord());

    for (unsigned i = 0; i < num_threads_; ++i) {
      files_.emplace_back(new std::ofstream((file_prefix_ + ".edges." + std::to_string(i)).c_str(),
                                            std::ofstream::binary | std::ofstream::out));
    }

    is_opened_ = true;
  }

  void Write(uint32_t *edge_ptr, int32_t bucket, int tid, Snapshot *snapshot) const {
    if (bucket != snapshot->bucket_id) {
      assert(snapshot->bucket_id == -1);
      assert(snapshot->p_rec.thread_id == -1);
      snapshot->bucket_id = bucket;
      snapshot->p_rec.thread_id = tid;
      snapshot->p_rec.starting_offset = cur_num_edges_[tid];
    }
    assert(snapshot->bucket_id == bucket);
    assert(snapshot->p_rec.thread_id == tid);

    files_[tid]->write(reinterpret_cast<const char *>(edge_ptr), sizeof(uint32_t) * words_per_edge_);
    ++snapshot->p_rec.total_number;
  }

  void SaveSnapshot(const Snapshot &snapshot) {
    if (snapshot.bucket_id != -1) {
      p_rec_[snapshot.bucket_id] = snapshot.p_rec;
      cur_num_edges_[snapshot.p_rec.thread_id] = snapshot.p_rec.total_number + snapshot.p_rec.starting_offset;
    }
  }

  void WriteUnsorted(uint32_t *edge_ptr, int tid) {
    files_[tid]->write(reinterpret_cast<const char *>(edge_ptr), sizeof(uint32_t) * words_per_edge_);
    ++num_unsorted_edges_[tid];
  }

  void Finalize() {
    if (is_opened_) {
      for (auto &file: files_) {
        file->close();
      }

      int64_t num_edges = 0;

      if (!unsorted_) {
        for (unsigned i = 0; i < p_rec_.size(); ++i) {
          num_edges += p_rec_[i].total_number;
        }
      } else {
        for (unsigned i = 0; i < num_unsorted_edges_.size(); ++i) {
          num_edges += num_unsorted_edges_[i];
        }
      }

      std::ofstream info_file(file_prefix_ + ".edges.info");
      info_file << "kmer_size " << kmer_size_ << '\n'
                << "words_per_edge " << words_per_edge_ << '\n'
                << "num_threads " << num_threads_ << '\n'
                << "num_buckets " << num_buckets_ << '\n'
                << "item_count " << num_edges << '\n';

      for (unsigned i = 0; i < p_rec_.size(); ++i) {
        info_file << i << ' '
                  << p_rec_[i].thread_id << ' '
                  << p_rec_[i].starting_offset << ' '
                  << p_rec_[i].total_number << '\n';
      }

      for (unsigned i = 0; i < num_unsorted_edges_.size(); ++i) {
        info_file << i << ' ' << num_unsorted_edges_[i] << '\n';
      }

      info_file.close();
      files_.clear();
      cur_bucket_.clear();
      cur_num_edges_.clear();
      p_rec_.clear();
      is_opened_ = false;
    }
  }
};

class MegahitEdgeReader {
 private:
  uint32_t kmer_size_{};
  uint32_t words_per_edge_{};
  int32_t num_files_{};
  int32_t num_buckets_{};
  int64_t num_edges_{};

  std::string file_prefix_;
  std::vector<BufferedReader> readers_;
  std::vector<std::unique_ptr<std::ifstream>> in_streams_;
  std::ifstream unsorted_stream_;
  BufferedReader unsorted_reader_;
  std::vector<PartitionRecord> p_rec_;
  std::vector<int64_t> file_sizes_;
  std::vector<uint32_t> buffer_;

  int cur_bucket_{};
  int cur_file_num_{};  // used for unsorted edges
  int64_t cur_cnt_{};
  int64_t cur_vol_{};
  BufferedReader *cur_reader_{};

  bool is_opened_;

 public:
  MegahitEdgeReader() : is_opened_(false) {}
  ~MegahitEdgeReader() { Close(); }

  void SetFilePrefix(const std::string &prefix) { file_prefix_ = prefix; }

  template<typename T>
  static void ScanField(std::ifstream &in, const std::string &field, T &out) {
    std::string s;
    in >> s >> out;
    assert(s == field);
  }

  void ReadInfo() {
    std::ifstream meta_file(file_prefix_ + ".edges.info");
    ScanField(meta_file, "kmer_size", kmer_size_);
    ScanField(meta_file, "words_per_edge", words_per_edge_);
    ScanField(meta_file, "num_threads", num_files_);
    ScanField(meta_file, "num_buckets", num_buckets_);
    ScanField(meta_file, "item_count", num_edges_);

    p_rec_.resize(num_buckets_);
    file_sizes_.resize(num_files_, 0);

    for (int i = 0; i < num_buckets_; ++i) {
      int b_id;
      meta_file >> b_id >> p_rec_[i].thread_id >> p_rec_[i].starting_offset
                >> p_rec_[i].total_number;
      if (b_id != i) {
        xfatal("Invalid format: bucket id not matched!\n");
      }
      if (p_rec_[i].thread_id >= num_files_) {
        xfatal("Record ID %d is greater than number of files %d\n", p_rec_[i].thread_id, num_files_);
      }
      if (p_rec_[i].thread_id >= 0) {
        file_sizes_[p_rec_[i].thread_id] += p_rec_[i].total_number;
      }
    }

    if (num_buckets_ == 0) {
      for (int i = 0; i < num_files_; ++i) {
        int b_id;
        if (!(meta_file >> b_id >> file_sizes_[i]) || b_id != i) {
          xfatal("Invalid format\n");
        }
      }
    }

    meta_file.close();
  }

  void InitFiles() {
    assert(!is_opened_);
    ReadInfo();
    buffer_.resize(words_per_edge_);

    for (int i = 0; i < num_files_; ++i) {
      in_streams_.emplace_back(new std::ifstream(file_prefix_ + ".edges." + std::to_string(i),
                                                 std::ifstream::binary | std::ifstream::in));
      readers_.emplace_back(BufferedReader());
      readers_.back().reset(in_streams_.back().get());
    }

    cur_cnt_ = 0;
    cur_vol_ = 0;
    cur_bucket_ = -1;

    // for unsorted
    if (IsUnsorted()) {
      cur_file_num_ = -1;
    }

    is_opened_ = true;
  }

  bool IsUnsorted() const { return num_buckets_ == 0; }
  uint32_t k() const { return kmer_size_; }
  uint32_t words_per_edge() const { return words_per_edge_; }
  int64_t num_edges() const { return num_edges_; }

  uint32_t *NextSortedEdge() {
    if (cur_bucket_ >= num_buckets_) {
      return nullptr;
    }

    while (cur_cnt_ >= cur_vol_) {
      ++cur_bucket_;

      while (cur_bucket_ < num_buckets_ && p_rec_[cur_bucket_].thread_id < 0) {
        ++cur_bucket_;
      }

      if (cur_bucket_ >= num_buckets_) {
        return nullptr;
      }

      cur_cnt_ = 0;
      cur_vol_ = p_rec_[cur_bucket_].total_number;
      cur_reader_ = &readers_[p_rec_[cur_bucket_].thread_id];
    }

    ++cur_cnt_;
    auto n_read = cur_reader_->read(buffer_.data(), words_per_edge_);
    assert(n_read == words_per_edge_ * sizeof(uint32_t));
    (void)(n_read);
    return buffer_.data();
  }

  uint32_t *NextUnsortedEdge() {
    if (cur_file_num_ >= num_files_) {
      return nullptr;
    }

    while (cur_cnt_ >= cur_vol_) {
      ++cur_file_num_;
      if (cur_file_num_ >= num_files_) {
        return nullptr;
      }
      unsorted_stream_.open(file_prefix_ + ".edges." + std::to_string(cur_file_num_),
                            std::ifstream::binary | std::ifstream::in);
      unsorted_reader_.reset(&unsorted_stream_);
      cur_cnt_ = 0;
      cur_vol_ = file_sizes_[cur_file_num_];
    }

    ++cur_cnt_;
    unsorted_reader_.read(buffer_.data(), words_per_edge_);
    return buffer_.data();
  }

  void Close() {
  }
};

#endif // MEGAHIT_EDGE_IO_H