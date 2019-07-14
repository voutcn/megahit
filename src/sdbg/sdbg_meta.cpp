//
// Created by vout on 11/5/18.
//

#include "sdbg_meta.h"

#include <algorithm>
#include <cassert>
#include <fstream>
#include "utils/utils.h"

SdbgMeta &SdbgMeta::FromBucketRecord(
    const std::vector<SdbgBucketRecord> &bucket_rec, uint32_t k,
    uint32_t words_per_tip_label) {
  bucket_rec_ = bucket_rec;
  k_ = k;
  words_per_tip_label_ = words_per_tip_label;
  item_count_ = 0;
  tip_count_ = 0;
  ones_in_last_ = 0;
  num_files_ = 0;
  large_mul_count_ = 0;
  std::fill(w_count_, w_count_ + kWAlphabetSize, 0);
  std::sort(bucket_rec_.begin(), bucket_rec_.end(),
            [](const SdbgBucketRecord &a, const SdbgBucketRecord &b) {
              return a.bucket_id < b.bucket_id;
            });
  for (auto &b : bucket_rec_) {
    if (b.file_id == SdbgBucketRecord::kNullID ||
        b.bucket_id == SdbgBucketRecord::kNullID) {
      continue;
    }
    b.accumulate_item_count = item_count_;
    b.accumulate_tip_count = tip_count_;
    item_count_ += b.num_items;
    tip_count_ += b.num_tips;
    large_mul_count_ += b.num_large_mul;
    ones_in_last_ += b.ones_in_last;
    num_files_ = std::max(num_files_, b.file_id + 1);
    for (unsigned w = 0; w < kWAlphabetSize; ++w) {
      w_count_[w] += b.num_w[w];
    }
  }
  std::sort(bucket_rec_.begin(), bucket_rec_.end(),
            [](const SdbgBucketRecord &a, const SdbgBucketRecord &b) {
              if (a.file_id != b.file_id) return a.file_id < b.file_id;
              return a.starting_offset < b.starting_offset;
            });
  return *this;
}
void SdbgMeta::Serialize(std::ofstream &os) {
  os << "k " << k_ << "\n"
     << "words_per_tip_label " << words_per_tip_label_ << "\n"
     << "num_buckets " << bucket_rec_.size() << "\n"
     << "num_files " << num_files_ << "\n";
  for (auto &b : bucket_rec_) {
    os << b.bucket_id << ' ' << b.file_id << ' ' << b.starting_offset << ' '
       << b.num_items << ' ' << b.num_tips << ' ' << b.num_large_mul << '\n';
  }
  os.close();
}

SdbgMeta &SdbgMeta::Deserialize(std::ifstream &is) {
  size_t num_buckets;
  ScanField(is, "k", k_);
  ScanField(is, "words_per_tip_label", words_per_tip_label_);
  ScanField(is, "num_buckets", num_buckets);
  ScanField(is, "num_files", num_files_);
  std::vector<SdbgBucketRecord> bucket_rec(num_buckets);
  for (size_t i = 0; i < num_buckets; ++i) {
    is >> bucket_rec[i].bucket_id >> bucket_rec[i].file_id >>
        bucket_rec[i].starting_offset >> bucket_rec[i].num_items >>
        bucket_rec[i].num_tips >> bucket_rec[i].num_large_mul;
  }
  return FromBucketRecord(bucket_rec, k_, words_per_tip_label_);
}