//
// Created by vout on 5/11/19.
//

#ifndef MEGAHIT_PAIR_END_FASTX_READER_H
#define MEGAHIT_PAIR_END_FASTX_READER_H

#include <memory>
#include "fastx_reader.h"

class PairedFastxReader : public BaseSequenceReader {
 public:
  PairedFastxReader(const std::string &file1, const std::string &file2) {
    readers_[0].reset(new FastxReader(file1));
    readers_[1].reset(new FastxReader(file2));
  }

  int64_t Read(SeqPackage *pkg, int64_t max_num, int64_t max_num_bases,
               bool reverse) override;

 private:
  std::unique_ptr<FastxReader> readers_[2];
  bool trim_n_{true};
};

#endif  // MEGAHIT_PAIR_END_FASTX_READER_H
