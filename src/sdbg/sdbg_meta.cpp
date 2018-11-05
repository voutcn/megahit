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

SdbgMetadata::SdbgMetadata(const std::string &file_name) {
  std::ifstream meta_file(file_name);
  ScanField(meta_file, "k", k_);
  ScanField(meta_file, "words_per_tip_label", words_per_tip_label_);
  ScanField(meta_file, "num_buckets", num_buckets_);
  ScanField(meta_file, "num_files", num_files_);
  bucket_rec_.resize(num_buckets_);
  size_ = 0;
  num_tips_ = 0;
  for (size_t i = 0; i < num_buckets_; ++i) {
    bucket_rec_[i].accumulate_item_count = size_;
    bucket_rec_[i].accumulate_tip_count = num_tips_;
    meta_file >> bucket_rec_[i].bucket_id
              >> bucket_rec_[i].file_id
              >> bucket_rec_[i].starting_offset
              >> bucket_rec_[i].num_items
              >> bucket_rec_[i].num_tips
              >> bucket_rec_[i].num_large_mul;
    assert(bucket_rec_[i].bucket_id == i);
    size_ += bucket_rec_[i].num_items;
    num_tips_ += bucket_rec_[i].num_tips;
    num_large_mul_ += bucket_rec_[i].num_large_mul;
  }
  bucket_rec_sorted_ = bucket_rec_;
  std::sort(bucket_rec_sorted_.begin(), bucket_rec_sorted_.end(),
            [](const SdbgBucketRecord &a, const SdbgBucketRecord &b) {
              if (a.file_id != b.file_id) return a.file_id < b.file_id;
              return a.starting_offset < b.starting_offset;
            }
  );
}