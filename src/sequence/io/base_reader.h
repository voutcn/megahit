//
// Created by vout on 5/11/19.
//

#ifndef MEGAHIT_BASE_READER_H
#define MEGAHIT_BASE_READER_H

#include "sequence/sequence_package.h"

class BaseSequenceReader {
 public:
  BaseSequenceReader() = default;
  virtual ~BaseSequenceReader() = default;
  static const int64_t kMaxNumSeq = 1LL << 60;
  static const int64_t kMaxNumBases = 1LL << 60;
  virtual int64_t Read(SeqPackage *pkg, int64_t max_num, int64_t max_num_bases,
                       bool reverse = false) = 0;
  int64_t ReadAll(SeqPackage *pkg, bool reverse) {
    return Read(pkg, kMaxNumSeq, kMaxNumBases, reverse);
  }
};

#endif  // MEGAHIT_BASE_READER_H
