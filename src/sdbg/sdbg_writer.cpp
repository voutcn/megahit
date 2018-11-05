//
// Created by vout on 11/5/18.
//

#include "sdbg_writer.h"

#include <cassert>

void SdbgWriter::InitFiles() {
  files_.resize(num_threads_);
  cur_bucket_.resize(num_threads_, -1);
  cur_thread_offset_.resize(num_threads_, 0);
  bucket_rec_.resize(num_buckets_);

  for (size_t i = 0; i < num_threads_; ++i) {
    files_[i] = std::make_shared<std::ofstream>((file_prefix_ + ".sdbg." + std::to_string(i)).c_str(),
                                                std::ofstream::binary | std::ofstream::out);
    assert(files_[i]->is_open());
  }
  is_opened_ = true;
}

void SdbgWriter::Write(unsigned tid,
                       int32_t bucket,
                       int w,
                       int last,
                       int tip,
                       mul_t multiplicity,
                       LabelWordType *packed_tip_label) {
  assert(tid < num_threads_);

  if (bucket != cur_bucket_[tid]) {
    cur_bucket_[tid] =
        bucket;
    assert(bucket_rec_[bucket].file_id == SdbgBucketRecord::kUninitializedFileID);
    bucket_rec_[bucket].file_id = tid;
    bucket_rec_[bucket].starting_offset = cur_thread_offset_[tid];
  }

  uint16_t packed_item = w | (last << 4) | (tip << 5) |
      (std::min(multiplicity, mul_t(kSmallMulSentinel)) << 8);
  files_[tid]->write(reinterpret_cast
                         <const char *>(&packed_item),
                     sizeof(packed_item));
  ++bucket_rec_[bucket].
      num_items;
  ++bucket_rec_[bucket].num_w[w];
  bucket_rec_[bucket].num_last1 +=
      last;
  cur_thread_offset_[tid] += sizeof(uint16_t);

  if (multiplicity > kMaxSmallMul) {
    files_[tid]->write(reinterpret_cast
                           <const char *>(&multiplicity),
                       sizeof(multiplicity));
    ++bucket_rec_[bucket].
        num_large_mul;
    cur_thread_offset_[tid] += sizeof(mul_t);
  }

  if (tip) {
    files_[tid]->write(reinterpret_cast
                           <const char *>(packed_tip_label),
                       sizeof(LabelWordType) * words_per_tip_label_);
    ++bucket_rec_[bucket].num_tips;
    cur_thread_offset_[tid] += sizeof(uint32_t) *
        words_per_tip_label_;
  }
}


void SdbgWriter::Finalize() {
  if (is_opened_) {
    for (size_t i = 0; i < num_threads_; ++i) {
      files_[i]->close();
    }

    std::ofstream meta_file((file_prefix_ + ".sdbg_info").c_str());
    meta_file << "k " << k_ << "\n"
              << "words_per_tip_label " << words_per_tip_label_ << "\n"
              << "num_buckets " << num_buckets_ << "\n"
              << "num_files " << num_threads_ << "\n";
    for (size_t i = 0; i < num_buckets_; ++i) {
      meta_file << i << ' '
                << bucket_rec_[i].file_id << ' '
                << bucket_rec_[i].starting_offset << ' '
                << bucket_rec_[i].num_items << ' '
                << bucket_rec_[i].num_tips << ' '
                << bucket_rec_[i].num_large_mul << '\n';
    }
    meta_file.close();
    files_.clear();
    cur_bucket_.clear();
    cur_thread_offset_.clear();
    bucket_rec_.clear();
    is_opened_ = false;
  }
}