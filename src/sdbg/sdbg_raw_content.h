//
// Created by vout on 11/5/18.
//

#ifndef MEGAHIT_SDBG_RAW_CONTENT_H
#define MEGAHIT_SDBG_RAW_CONTENT_H

#include "sdbg_def.h"
#include "sdbg_meta.h"
#include <kmlib/kmcompactvector.h>
#include <sparsepp/sparsepp/spp.h>

/**
 * The raw (non-indexed) data of a SDBG
 */
struct SdbgRawContent {
  SdbgRawContent() = default;
  ~SdbgRawContent() = default;
  SdbgMetadata meta;
  kmlib::CompactVector<kAlphabetSize, uint64_t> w;
  kmlib::CompactVector<1, uint64_t> last, tip;
  std::vector<small_mul_t> small_mul;
  std::vector<label_word_t> tip_lables;
  spp::sparse_hash_map<uint64_t, mul_t> large_mul;
  std::vector<mul_t> full_mul;
};


void ReadSdbgFromFile(const std::string &file_prefix, SdbgRawContent *raw_content);

#endif //MEGAHIT_SDBG_RAW_CONTENT_H
