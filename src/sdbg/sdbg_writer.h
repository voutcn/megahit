//
// Created by vout on 11/5/18.
//

#ifndef MEGAHIT_SDBG_WRITER_H
#define MEGAHIT_SDBG_WRITER_H

#include "sdbg_def.h"
#include "sdbg_meta.h"

#include <fstream>
#include <memory>
#include <string>
#include <vector>

/**
 * A SDBG writer is used to write partitioned SDBG to files
 */
class SdbgWriter {
 public:
  struct Snapshot {
   private:
    uint64_t cur_thread_offset{0};
    SdbgBucketRecord bucket_record{};

    friend class SdbgWriter;
  };

  SdbgWriter() = default;
  ~SdbgWriter() { Finalize(); }

  void set_num_threads(size_t num_threads) { num_threads_ = num_threads; }
  void set_file_prefix(const std::string &file_prefix) {
    file_prefix_ = file_prefix;
  }
  void set_kmer_size(unsigned k) {
    k_ = k;
    words_per_tip_label_ = (k + kCharsPerLabelWord - 1) / kCharsPerLabelWord;
  }
  void set_num_buckets(size_t num_buckets) { num_buckets_ = num_buckets; }

  void InitFiles();
  void Write(unsigned tid, uint32_t bucket_id, uint8_t w, uint8_t last,
             uint8_t tip, mul_t multiplicity, label_word_t *packed_tip_label,
             Snapshot *snapshot);
  void SaveSnapshot(const Snapshot &snapshot);
  void Finalize();
  const SdbgMeta &final_meta() const { return final_meta_; }

 private:
  std::string file_prefix_;
  size_t num_threads_{};
  size_t num_buckets_{};

  std::vector<std::unique_ptr<std::ofstream>> files_;
  std::vector<uint64_t> cur_thread_offset_;  // offset in BYTE
  std::vector<SdbgBucketRecord> bucket_rec_;

  bool is_opened_{};
  unsigned k_{};
  size_t words_per_tip_label_{};
  SdbgMeta final_meta_;
};

#endif  // MEGAHIT_SDBG_WRITER_H
