//
// Created by vout on 11/5/18.
//

#include "sdbg_raw_content.h"
#include "sdbg_meta.h"

#include <fstream>
#include "utils/buffered_reader.h"
#include "sdbg_item.h"

SdbgRawContent ReadSdbgFromFile(const std::string &file_prefix) {
  SdbgRawContent raw_content;
  std::ifstream is((file_prefix + ".sdbg_info").c_str());
  raw_content.meta.Deserialize(is);
  const auto &metadata = raw_content.meta;
  raw_content.w.resize(metadata.size());
  raw_content.last.resize(metadata.size());
  raw_content.tip.resize(metadata.size());
  raw_content.tip_lables.resize(metadata.words_per_tip_label() * metadata.num_tips());
  bool use_full_mul = metadata.num_large_mul() >= metadata.size() * 0.08;
  if (use_full_mul) {
    raw_content.full_mul.resize(metadata.size());
  } else {
    raw_content.small_mul.resize(metadata.size());
    raw_content.large_mul.reserve(metadata.num_large_mul());
  }

  BufferedReader in;
  size_t last_file_id = SdbgBucketRecord::kUninitializedFileID;
  size_t file_offset = 0;

  for (auto bucket_it = metadata.sorted_begin();
       bucket_it != metadata.sorted_end() && bucket_it->file_id != SdbgBucketRecord::kUninitializedFileID;
       ++bucket_it) {
    if (bucket_it->file_id != last_file_id) {
      if (is.is_open()) {
        is.close();
      }
      file_offset = 0;
      is.open((file_prefix + ".sdbg." + std::to_string(bucket_it->file_id)).c_str());
      in.reset(&is);
      last_file_id = bucket_it->file_id;
    }
    assert(file_offset == bucket_it->starting_offset);
    SdbgItem item;
    label_word_t *tip_label_ptr = &raw_content.tip_lables[
        bucket_it->accumulate_tip_count * metadata.words_per_tip_label()];

    for (size_t i = 0; i < bucket_it->num_items; ++i) {
      size_t index = i + bucket_it->accumulate_item_count;
      file_offset += in.read(&item);
      raw_content.w[index] = item.w;
      raw_content.last[index] = item.last;
      raw_content.tip[index] = item.tip;
      mul_t mul = item.mul;
      if (mul == kSmallMulSentinel) {
        file_offset += in.read(&mul);
      }
      if (use_full_mul) {
        raw_content.full_mul[index] = mul;
      } else {
        if (mul < kMaxSmallMul) {
          raw_content.small_mul[index] = mul;
        } else {
          raw_content.small_mul[index] = kSmallMulSentinel;
          raw_content.large_mul[index] = mul;
        }
      }
      if (item.tip) {
        file_offset += in.read(tip_label_ptr, metadata.words_per_tip_label());
        tip_label_ptr += metadata.words_per_tip_label();
      }
    }
    assert(tip_label_ptr - raw_content.tip_lables.data() ==
        ptrdiff_t((bucket_it->accumulate_tip_count + bucket_it->num_tips) * metadata.words_per_tip_label()));
  }
  return raw_content;
}