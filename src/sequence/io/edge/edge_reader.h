//
// Created by vout on 5/11/19.
//

#ifndef MEGAHIT_EDGE_READER_H
#define MEGAHIT_EDGE_READER_H

#include "edge_io_meta.h"
#include "sequence/io/base_reader.h"

class EdgeReader : public BaseSequenceReader {
 public:
  EdgeReader(const std::string &file_prefix) {
    file_prefix_ = file_prefix;
    std::ifstream is(file_prefix + ".edges.info");
    metadata_.Deserialize(is);
    InitFiles();
  }

  EdgeReader *SetMultiplicityVec(std::vector<mul_t> *mul) {
    mul_ = mul;
    return this;
  }

  int64_t ReadUnsorted(SeqPackage *pkg, std::vector<mul_t> *mul,
                       int64_t max_num) {
    for (int64_t i = 0; i < max_num; ++i) {
      uint32_t *next_edge = NextUnsortedEdge();
      if (next_edge == nullptr) {
        return i;
      }
      pkg->AppendCompactSequence(next_edge, metadata_.kmer_size + 1);
      if (mul) {
        mul->push_back(next_edge[metadata_.words_per_edge - 1] & kMaxMul);
      }
    }
    return max_num;
  }

  int64_t ReadSorted(SeqPackage *pkg, std::vector<mul_t> *mul,
                     int64_t max_num) {
    for (int64_t i = 0; i < max_num; ++i) {
      uint32_t *next_edge = NextSortedEdge();
      if (next_edge == nullptr) {
        return i;
      }
      pkg->AppendCompactSequence(next_edge, metadata_.kmer_size + 1);
      if (mul) {
        mul->push_back(next_edge[metadata_.words_per_edge - 1] & kMaxMul);
      }
    }
    return max_num;
  }

  int64_t Read(SeqPackage *pkg, int64_t max_num, int64_t max_num_bases,
               bool reverse = false) override {
    if (metadata_.is_sorted) {
      return ReadSorted(pkg, mul_, max_num);
    } else {
      return ReadUnsorted(pkg, mul_, max_num);
    }
  }

 private:
  std::vector<mul_t> *mul_{nullptr};

  std::string file_prefix_;
  std::vector<std::unique_ptr<std::ifstream>> in_streams_;
  BufferedReader cur_reader_;
  std::vector<uint32_t> buffer_;

  int cur_bucket_{};
  int64_t cur_cnt_{};
  int64_t cur_vol_{};
  bool is_opened_{false};

  EdgeIoMetadata metadata_;

 private:
  void InitFiles() {
    assert(!is_opened_);
    buffer_.resize(metadata_.words_per_edge);

    for (unsigned i = 0; i < metadata_.num_files; ++i) {
      in_streams_.emplace_back(
          new std::ifstream(file_prefix_ + ".edges." + std::to_string(i),
                            std::ifstream::binary | std::ifstream::in));
    }

    cur_cnt_ = 0;
    cur_vol_ = 0;
    cur_bucket_ = -1;

    if (!metadata_.is_sorted) {
      cur_reader_.reset(in_streams_[0].get());
    }

    is_opened_ = true;
  }

 public:
  const EdgeIoMetadata &GetMetadata() const { return metadata_; }

 private:
  uint32_t *NextSortedEdge() {
    if (cur_bucket_ >= static_cast<int>(metadata_.buckets.size())) {
      return nullptr;
    }

    while (cur_cnt_ >= cur_vol_) {
      ++cur_bucket_;

      while (cur_bucket_ < static_cast<int>(metadata_.buckets.size()) &&
             metadata_.buckets[cur_bucket_].file_id < 0) {
        ++cur_bucket_;
      }

      if (cur_bucket_ >= static_cast<int>(metadata_.buckets.size())) {
        return nullptr;
      }

      const auto &bucket = metadata_.buckets[cur_bucket_];
      cur_cnt_ = 0;
      cur_vol_ = bucket.total_number;
      auto is = in_streams_[bucket.file_id].get();
      is->clear();
      is->seekg(bucket.file_offset * sizeof(uint32_t) *
                metadata_.words_per_edge);
      cur_reader_.reset(is, bucket.total_number * sizeof(uint32_t) *
                                metadata_.words_per_edge);
    }

    ++cur_cnt_;
    auto n_read = cur_reader_.read(buffer_.data(), metadata_.words_per_edge);
    assert(n_read == metadata_.words_per_edge * sizeof(uint32_t));
    (void)n_read;
    return buffer_.data();
  }

  uint32_t *NextUnsortedEdge() {
    if (cur_cnt_ >= metadata_.num_edges) {
      return nullptr;
    }

    ++cur_cnt_;
    auto n_read = cur_reader_.read(buffer_.data(), metadata_.words_per_edge);
    assert(n_read == metadata_.words_per_edge * sizeof(uint32_t));
    (void)n_read;
    return buffer_.data();
  }
};

#endif  // MEGAHIT_EDGE_READER_H
