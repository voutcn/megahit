//
// Created by vout on 11/4/18.
//

#ifndef MEGAHIT_SDBG_META_H
#define MEGAHIT_SDBG_META_H

#include "sdbg_def.h"

#include <cstdint>
#include <cstddef>
#include <cstring>
#include <string>
#include <vector>
#include <limits>

/**
 * @brief bucket record of a SDBG
 * Edges with the same L-prefix (L=8) are grouped in one bucket
 * A bucket record contains the the statistics and metadata of the partition
 */
struct SdbgBucketRecord {
  static const size_t kUninitializedFileID = std::numeric_limits<size_t>::max();
  size_t file_id{kUninitializedFileID};
  size_t bucket_id{0};
  size_t accumulate_item_count{0};
  size_t accumulate_tip_count{0};
  size_t starting_offset{0};
  size_t num_items{0};
  size_t num_tips{0};
  size_t num_large_mul{0};
  size_t num_last1{0};
  size_t num_w[9]{};

  SdbgBucketRecord() {
    memset(num_w, 0, sizeof(num_w));
  }
};

/**
 * Metadata of a SDBG, including bucket profile
 */
class SdbgMetadata {
 private:
  size_t size_{};
  size_t k_{};
  size_t words_per_tip_label_{};
  size_t num_buckets_{};
  size_t num_files_{};
  size_t num_tips_{};
  size_t num_large_mul_{};
  std::vector<SdbgBucketRecord> bucket_rec_;
  std::vector<SdbgBucketRecord> bucket_rec_sorted_;
 public:
  SdbgMetadata() = default;
  explicit SdbgMetadata(const std::string &file_name);
  explicit SdbgMetadata(const std::vector<SdbgBucketRecord> &bucket_rec);
  ~SdbgMetadata() = default;
  size_t k() const {
    return k_;
  }
  size_t bucket_size() const {
    return bucket_rec_.size();
  }
  size_t size() const {
    return size_;
  }
  size_t num_tips() const {
    return num_tips_;
  }
  size_t num_large_mul() const {
    return num_large_mul_;
  }
  size_t words_per_tip_label() const {
    return words_per_tip_label_;
  }
  std::vector<SdbgBucketRecord>::const_iterator sorted_begin() {
    return bucket_rec_sorted_.begin();
  }
  std::vector<SdbgBucketRecord>::const_iterator sorted_end() {
    return bucket_rec_sorted_.end();
  }
};

#endif //MEGAHIT_SDBG_META_H
