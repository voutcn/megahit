//
// Created by vout on 11/5/18.
//

#include "sdbg_raw_content.h"
#include "sdbg_meta.h"

#include <fstream>

SdbgRawContent ReadSdbgFromFile(const std::string &file_prefix) {
  SdbgMetadata metadata(file_prefix + ".sdbg_info");
  SdbgRawContent raw_content;
  raw_content.meta = metadata;
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

  std::ifstream in;
  char buffer[4096];
  size_t last_file_id = SdbgBucketRecord::kUninitializedFileID;
  size_t file_offset = 0;

  for (auto bucket_it = metadata.sorted_begin();
       bucket_it != metadata.sorted_end() && bucket_it->file_id != SdbgBucketRecord::kUninitializedFileID;
       ++bucket_it) {
    if (bucket_it->file_id != last_file_id) {
      if (in.is_open()) {
        in.close();
      }
      file_offset = 0;
      in.rdbuf()->pubsetbuf(buffer, 4096);
      in.open((file_prefix + ".sdbg." + std::to_string(bucket_it->file_id)).c_str());
      last_file_id = bucket_it->file_id;
    }
    assert(file_offset == bucket_it->starting_offset);
    uint16_t packed_item;
    LabelWordType *tip_label_ptr = &raw_content.tip_lables[
        bucket_it->accumulate_tip_count * metadata.words_per_tip_label()];
    for (size_t i = 0; i < bucket_it->num_items; ++i) {
      size_t index = i + bucket_it->accumulate_item_count;
      in.read(reinterpret_cast<char *>(&packed_item), sizeof(packed_item));
      file_offset += sizeof(packed_item);
      auto w = packed_item & 0xFU;
      auto last = (packed_item >> 4) & 1U;
      auto tip = (packed_item >> 5) & 1U;
      mul_t mul = packed_item >> 8;
      raw_content.w[index] = w;
      raw_content.last[index] = last;
      raw_content.tip[index] = tip;
      if (mul == kSmallMulSentinel) {
        in.read(reinterpret_cast<char *>(&mul), sizeof(mul_t));
        file_offset += sizeof(mul_t);
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
      if (tip) {
        in.read(reinterpret_cast<char *>(tip_label_ptr), sizeof(LabelWordType) * metadata.words_per_tip_label());
        tip_label_ptr += metadata.words_per_tip_label();
        file_offset += sizeof(LabelWordType) * metadata.words_per_tip_label();
      }
    }
    assert(tip_label_ptr - raw_content.tip_lables.data() ==
        ptrdiff_t((bucket_it->accumulate_tip_count + bucket_it->num_tips) * metadata.words_per_tip_label()));
  }
  return raw_content;
}