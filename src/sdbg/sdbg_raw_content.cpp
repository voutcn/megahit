//
// Created by vout on 11/5/18.
//

#include "sdbg_raw_content.h"
#include "sdbg_meta.h"

#include <fstream>
#include "kmlib/kmbit.h"
#include "sdbg_item.h"
#include "utils/buffered_reader.h"

/**
 * load SDBG raw content from disk
 * @param raw_content the pointer the raw_content to be stored to
 * @param file_prefix the prefix of the SDBG files
 */
void LoadSdbgRawContent(SdbgRawContent *raw_content,
                        const std::string &file_prefix) {
  std::ifstream is(file_prefix + ".sdbg_info");
  raw_content->meta.Deserialize(is);
  const auto &metadata = raw_content->meta;
  raw_content->w.resize(metadata.item_count());
  raw_content->last.resize(metadata.item_count());
  raw_content->tip.resize(metadata.item_count());
  raw_content->tip_lables.resize(metadata.words_per_tip_label() *
                                 metadata.tip_count());
  bool use_full_mul =
      metadata.large_mul_count() >= metadata.item_count() * 0.08;
  if (use_full_mul) {
    raw_content->full_mul.resize(metadata.item_count());
  } else {
    raw_content->small_mul.resize(metadata.item_count());
    raw_content->large_mul.reserve(metadata.large_mul_count());
  }

  BufferedReader in;
  size_t last_file_id = SdbgBucketRecord::kNullID;
  size_t file_offset = 0;

  for (auto bucket_it = metadata.begin_bucket();
       bucket_it != metadata.end_bucket() &&
       bucket_it->bucket_id != bucket_it->kNullID;
       ++bucket_it) {
    if (bucket_it->file_id != last_file_id) {
      if (is.is_open()) {
        is.close();
      }
      file_offset = 0;
      assert(in.read<char>(nullptr) == 0);
      is.open((file_prefix + ".sdbg." + std::to_string(bucket_it->file_id))
                  .c_str());
      in.reset(&is);
      last_file_id = bucket_it->file_id;
    }
    assert(file_offset == bucket_it->starting_offset);
    SdbgItem item;
    label_word_t *tip_label_ptr =
        raw_content->tip_lables.data() +
        bucket_it->accumulate_tip_count * metadata.words_per_tip_label();

    for (size_t i = 0; i < bucket_it->num_items; ++i) {
      size_t index = i + bucket_it->accumulate_item_count;
      file_offset += in.read(&item);
      raw_content->w[index] = item.w;
      raw_content->last[index] = item.last;
      raw_content->tip[index] = item.tip;
      mul_t mul = item.mul;
      if (mul == kSmallMulSentinel) {
        file_offset += in.read(&mul);
      }
      if (use_full_mul) {
        raw_content->full_mul[index] = mul;
      } else {
        if (mul < kMaxSmallMul) {
          raw_content->small_mul[index] = mul;
        } else {
          raw_content->small_mul[index] = kSmallMulSentinel;
          raw_content->large_mul[index] = mul;
        }
      }
      if (item.tip) {
        file_offset += in.read(tip_label_ptr, metadata.words_per_tip_label());
        for (unsigned j = 0; j < metadata.words_per_tip_label(); ++j) {
          tip_label_ptr[j] =
              kmlib::bit::Reverse<kBitsPerChar>(tip_label_ptr[j]);
        }
        tip_label_ptr += metadata.words_per_tip_label();
      }
    }
    assert(tip_label_ptr - raw_content->tip_lables.data() ==
           ptrdiff_t((bucket_it->accumulate_tip_count + bucket_it->num_tips) *
                     metadata.words_per_tip_label()));
  }
  assert(in.read<char>(nullptr) == 0);
}