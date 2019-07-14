#ifndef MEGAHIT_EDGE_WRITER_H
#define MEGAHIT_EDGE_WRITER_H

#include <cassert>
#include <cstdint>

#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "definitions.h"
#include "edge_io_meta.h"
#include "utils/buffered_reader.h"
#include "utils/utils.h"

class EdgeWriter {
 private:
  EdgeIoMetadata metadata_{};
  std::string file_prefix_;
  std::vector<std::unique_ptr<std::ofstream>> files_;
  std::vector<int64_t> n_edges_at_thread_;

  bool is_opened_;

 public:
  class Snapshot {
   private:
    EdgeIoBucketInfo bucket_info;
    int bucket_id{-1};
    friend class EdgeWriter;
  };

  EdgeWriter() : is_opened_(false){};
  ~EdgeWriter() { Finalize(); }

  void SetKmerSize(uint32_t k) {
    metadata_.kmer_size = k;
    metadata_.words_per_edge = DivCeiling((k + 1) * 2 + 16, 32);
  }

  void SetNumThreads(int32_t num_threads) { metadata_.num_files = num_threads; }
  void SetFilePrefix(const std::string &prefix) { file_prefix_ = prefix; }
  void SetNumBuckets(int num_buckets) {
    metadata_.buckets.clear();
    metadata_.buckets.resize(num_buckets);
  }
  void SetUnordered() {
    metadata_.buckets.clear();
    metadata_.is_sorted = false;
    metadata_.num_files = 1;
    metadata_.num_edges = 0;
  }

  void InitFiles() {
    assert(!is_opened_);
    n_edges_at_thread_.resize(metadata_.num_files, 0);

    for (unsigned i = 0; i < metadata_.num_files; ++i) {
      files_.emplace_back(new std::ofstream(
          (file_prefix_ + ".edges." + std::to_string(i)).c_str(),
          std::ofstream::binary | std::ofstream::out));
    }

    is_opened_ = true;
  }

  void Write(uint32_t *edge_ptr, int32_t bucket, int tid,
             Snapshot *snapshot) const {
    assert(metadata_.is_sorted);
    if (bucket != snapshot->bucket_id) {
      assert(snapshot->bucket_id == -1);
      assert(snapshot->bucket_info.file_id == -1);
      snapshot->bucket_id = bucket;
      snapshot->bucket_info.file_id = tid;
      snapshot->bucket_info.file_offset = n_edges_at_thread_[tid];
    }
    assert(snapshot->bucket_id == bucket);
    assert(snapshot->bucket_info.file_id == tid);

    files_[tid]->write(reinterpret_cast<const char *>(edge_ptr),
                       sizeof(uint32_t) * metadata_.words_per_edge);
    ++snapshot->bucket_info.total_number;
  }

  void SaveSnapshot(const Snapshot &snapshot) {
    if (snapshot.bucket_id != -1) {
      metadata_.buckets[snapshot.bucket_id] = snapshot.bucket_info;
      n_edges_at_thread_[snapshot.bucket_info.file_id] =
          snapshot.bucket_info.total_number + snapshot.bucket_info.file_offset;
    }
  }

  void WriteUnordered(uint32_t *edge_ptr) {
    assert(!metadata_.is_sorted);
    files_[0]->write(reinterpret_cast<const char *>(edge_ptr),
                     sizeof(uint32_t) * metadata_.words_per_edge);
    ++metadata_.num_edges;
  }

  void Finalize() {
    if (is_opened_) {
      for (auto &file : files_) {
        file->close();
      }
      std::ofstream info_file(file_prefix_ + ".edges.info");
      metadata_.Serialize(info_file);
      info_file.close();
      is_opened_ = false;
    }
  }
};

#endif  // MEGAHIT_EDGE_WRITER_H