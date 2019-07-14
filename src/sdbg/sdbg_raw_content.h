//
// Created by vout on 11/5/18.
//

#ifndef MEGAHIT_SDBG_RAW_CONTENT_H
#define MEGAHIT_SDBG_RAW_CONTENT_H

#include "kmlib/kmcompactvector.h"
#include "parallel_hashmap/phmap.h"
#include "sdbg_def.h"
#include "sdbg_meta.h"

/**
 * The raw (non-indexed) data of a SDBG
 */
struct SdbgRawContent {
  SdbgRawContent() = default;
  ~SdbgRawContent() = default;
  SdbgMeta meta;
  kmlib::CompactVector<kAlphabetSize, uint64_t> w;
  kmlib::CompactVector<1, uint64_t> last, tip;
  std::vector<small_mul_t> small_mul;
  std::vector<label_word_t> tip_lables;
  phmap::parallel_flat_hash_map<uint64_t, mul_t> large_mul;
  std::vector<mul_t> full_mul;
};

void LoadSdbgRawContent(SdbgRawContent *raw_content,
                        const std::string &file_prefix);

#endif  // MEGAHIT_SDBG_RAW_CONTENT_H
