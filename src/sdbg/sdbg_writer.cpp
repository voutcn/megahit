//
// Created by vout on 11/5/18.
//

#include "sdbg_writer.h"

#include <cassert>
#include <algorithm>
#include "sdbg_item.h"

void SdbgWriter::InitFiles() {
  files_.resize(num_threads_);
  cur_bucket_.resize(num_threads_, -1);
  cur_thread_offset_.resize(num_threads_, 0);
  bucket_rec_.resize(num_buckets_);

  for (size_t i = 0; i < num_threads_; ++i) {
    files_[i] = std::make_shared<std::ofstream>(
        (file_prefix_ + ".sdbg." + std::to_string(i)).c_str(),
        std::ofstream::binary | std::ofstream::out);
    assert(files_[i]->is_open());
  }
  is_opened_ = true;
}

void SdbgWriter::Write(unsigned tid,
                       uint32_t bucket_id,
                       uint8_t w,
                       uint8_t last,
                       uint8_t tip,
                       mul_t multiplicity,
                       label_word_t *packed_tip_label) {
  assert(tid < num_threads_);
  if (bucket_id != cur_bucket_[tid]) {
    cur_bucket_[tid] = bucket_id;
    assert(bucket_rec_[bucket_id].file_id == SdbgBucketRecord::kUninitializedFileID);
    bucket_rec_[bucket_id].file_id = tid;
    bucket_rec_[bucket_id].bucket_id = bucket_id;
    bucket_rec_[bucket_id].starting_offset = cur_thread_offset_[tid];
  }

  SdbgItem item(w, last, tip, std::min(multiplicity, mul_t{kSmallMulSentinel}));
  files_[tid]->write(reinterpret_cast<const char *>(&item),
                     sizeof(item));
  ++bucket_rec_[bucket_id].num_items;
  ++bucket_rec_[bucket_id].num_w[w];
  bucket_rec_[bucket_id].ones_in_last += last;
  cur_thread_offset_[tid] += sizeof(uint16_t);

  if (multiplicity > kMaxSmallMul) {
    files_[tid]->write(reinterpret_cast<const char *>(&multiplicity),
                       sizeof(multiplicity));
    ++bucket_rec_[bucket_id].num_large_mul;
    cur_thread_offset_[tid] += sizeof(mul_t);
  }

  if (tip) {
    files_[tid]->write(reinterpret_cast<const char *>(packed_tip_label),
                       sizeof(label_word_t) * words_per_tip_label_);
    ++bucket_rec_[bucket_id].num_tips;
    cur_thread_offset_[tid] += sizeof(label_word_t) * words_per_tip_label_;
  }
}

void SdbgWriter::Finalize() {
  if (is_opened_) {
    for (size_t i = 0; i < num_threads_; ++i) {
      files_[i]->close();
    }

    std::ofstream os((file_prefix_ + ".sdbg_info").c_str());
    final_meta_.FromBucketRecord(bucket_rec_, k_, words_per_tip_label_).Serialize(os);
    os.close();
    is_opened_ = false;
  }
}