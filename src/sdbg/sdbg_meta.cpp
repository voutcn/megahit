//
// Created by vout on 11/5/18.
//

#include "sdbg_meta.h"

#include <cassert>
#include <fstream>
#include <algorithm>

/**
 * Helper function to scan name and value
 * @tparam T
 * @param in
 * @param field
 * @param out
 */
template<typename T>
static void ScanField(std::ifstream &in, const std::string &field, T &out) {
  std::string s;
  in >> s >> out;
  assert(s == field);
}

SdbgMetadata &SdbgMetadata::FromBucketRecord(
    const std::vector<SdbgBucketRecord> &bucket_rec, size_t k,
    size_t words_per_tip_label) {
  bucket_rec_ = bucket_rec;
  k_ = k;
  words_per_tip_label_ = words_per_tip_label;
  size_ = 0;
  num_tips_ = 0;
  num_files_ = 0;
  num_large_mul_ = 0;
  std::sort(bucket_rec_.begin(), bucket_rec_.end(),
            [](const SdbgBucketRecord &a, const SdbgBucketRecord &b) {
              return a.bucket_id < b.bucket_id;
            }
  );
  for (auto &b : bucket_rec_) {
    b.accumulate_item_count = size_;
    b.accumulate_tip_count = num_tips_;
    size_ += b.num_items;
    num_tips_ += b.num_tips;
    num_large_mul_ += b.num_large_mul;
    num_files_ = std::max(num_files_, b.file_id);
  }
  std::sort(bucket_rec_.begin(), bucket_rec_.end(),
            [](const SdbgBucketRecord &a, const SdbgBucketRecord &b) {
              if (a.file_id != b.file_id) return a.file_id < b.file_id;
              return a.starting_offset < b.starting_offset;
            }
  );
  return *this;
}
void SdbgMetadata::Serialize(std::ofstream &os) {
  os << "k " << k_ << "\n"
     << "words_per_tip_label " << words_per_tip_label_ << "\n"
     << "num_buckets " << bucket_rec_.size() << "\n"
     << "num_files " << num_files_ << "\n";
  for (auto &b : bucket_rec_) {
    os << b.bucket_id << ' '
       << b.file_id << ' '
       << b.starting_offset << ' '
       << b.num_items << ' '
       << b.num_tips << ' '
       << b.num_large_mul << '\n';
  }
  os.close();
}

SdbgMetadata &SdbgMetadata::Deserialize(std::ifstream &is) {
  size_t num_buckets;
  ScanField(is, "k", k_);
  ScanField(is, "words_per_tip_label", words_per_tip_label_);
  ScanField(is, "num_buckets", num_buckets);
  ScanField(is, "num_files", num_files_);
  bucket_rec_.resize(num_buckets);
  for (size_t i = 0; i < num_buckets; ++i) {
    is >> bucket_rec_[i].bucket_id
       >> bucket_rec_[i].file_id
       >> bucket_rec_[i].starting_offset
       >> bucket_rec_[i].num_items
       >> bucket_rec_[i].num_tips
       >> bucket_rec_[i].num_large_mul;
  }
  this->FromBucketRecord(bucket_rec_, k_, words_per_tip_label_);
  return *this;
}