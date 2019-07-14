//
// Created by vout on 11/4/18.
//

#ifndef MEGAHIT_SDBG_META_H
#define MEGAHIT_SDBG_META_H

#include "sdbg_def.h"

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <vector>

/**
 * @brief bucket record of a SDBG
 * Edges with the same L-prefix (L=8) are grouped in one bucket
 * A bucket record contains the the statistics and metadata of the partition
 */
struct SdbgBucketRecord {
  static const size_t kNullID = std::numeric_limits<size_t>::max();
  size_t file_id{kNullID};
  size_t bucket_id{kNullID};
  size_t accumulate_item_count{0};
  size_t accumulate_tip_count{0};
  size_t starting_offset{0};
  size_t num_items{0};
  size_t num_tips{0};
  size_t num_large_mul{0};
  size_t ones_in_last{0};
  size_t num_w[kWAlphabetSize]{};
};

/**
 * Metadata of a SDBG, including bucket profile
 */
class SdbgMeta {
 private:
  uint32_t k_{};
  uint32_t words_per_tip_label_{};
  size_t item_count_{};
  size_t num_files_{};
  size_t tip_count_{};
  size_t large_mul_count_{};
  size_t ones_in_last_{};
  size_t w_count_[kWAlphabetSize]{};
  std::vector<SdbgBucketRecord> bucket_rec_;

 public:
  SdbgMeta() = default;
  void Serialize(std::ofstream &os);
  SdbgMeta &Deserialize(std::ifstream &is);
  SdbgMeta &FromBucketRecord(const std::vector<SdbgBucketRecord> &bucket_rec,
                             uint32_t k, uint32_t words_per_tip_label);
  ~SdbgMeta() = default;
  uint32_t k() const { return k_; }
  size_t bucket_count() const { return bucket_rec_.size(); }
  size_t item_count() const { return item_count_; }
  size_t tip_count() const { return tip_count_; }
  size_t w_count(unsigned w) const { return w_count_[w]; }
  size_t ones_in_last() const { return ones_in_last_; }
  size_t large_mul_count() const { return large_mul_count_; }
  size_t words_per_tip_label() const { return words_per_tip_label_; }
  std::vector<SdbgBucketRecord>::const_iterator begin_bucket() const {
    return bucket_rec_.begin();
  }
  std::vector<SdbgBucketRecord>::const_iterator end_bucket() const {
    return bucket_rec_.end();
  }
};

#endif  // MEGAHIT_SDBG_META_H
