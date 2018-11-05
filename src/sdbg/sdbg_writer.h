//
// Created by vout on 11/5/18.
//

#ifndef MEGAHIT_SDBG_WRITER_H
#define MEGAHIT_SDBG_WRITER_H

#include "sdbg_def.h"

#include <string>
#include <vector>
#include <fstream>
#include <memory>

#include "sdbg_meta.h"

/**
 * A SDBG writer is used to write partitioned SDBG to files
 */
class SdbgWriter {
 private:
  std::string file_prefix_;
  size_t num_threads_;
  size_t num_buckets_;

  std::vector<std::shared_ptr<std::ofstream>> files_;
  std::vector<int> cur_bucket_;
  std::vector<int64_t> cur_thread_offset_;    // offset in BYTE
  std::vector<SdbgBucketRecord> bucket_rec_;

  bool is_opened_;
  unsigned k_;
  size_t words_per_tip_label_;

 public:

  SdbgWriter() {}
  ~SdbgWriter() {
    Finalize();
  }

  void set_num_threads(size_t num_threads) {
    num_threads_ = num_threads;
  }
  void set_file_prefix(const std::string &file_prefix) {
    file_prefix_ = file_prefix;
  }
  void set_kmer_size(unsigned k) {
    k_ = k;
    words_per_tip_label_ = (k + kCharsPerLabelWord - 1) / kCharsPerLabelWord;
  }
  void set_num_buckets(size_t num_buckets) {
    num_buckets_ = num_buckets;
  }

  void InitFiles();
  void Write(unsigned tid, int32_t bucket, int w, int last, int tip, mul_t multiplicity,
             LabelWordType *packed_tip_label);
  void Finalize();
  int64_t num_edges() {
    int64_t total_edges = 0;

    for (size_t i = 0; i < num_buckets_; ++i) {
      total_edges += bucket_rec_[i].num_items;
    }

    return total_edges;
  }

  int64_t num_w(int w) {
    int64_t ret = 0;
    for (size_t i = 0; i < num_buckets_; ++i) {
      ret += bucket_rec_[i].num_w[w];
    }
    return ret;
  }

  int64_t num_last1() {
    int64_t ret = 0;
    for (size_t i = 0; i < num_buckets_; ++i) {
      ret += bucket_rec_[i].num_last1;
    }
    return ret;
  }

  int64_t num_tips() {
    int64_t ret = 0;
    for (size_t i = 0; i < num_buckets_; ++i) {
      ret += bucket_rec_[i].num_tips;
    }
    return ret;
  }

};

#endif //MEGAHIT_SDBG_WRITER_H
