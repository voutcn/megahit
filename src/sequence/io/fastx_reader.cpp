//
// Created by vout on 4/28/19.
//

#include "fastx_reader.h"
#include <cassert>
#include <stdexcept>

FastxReader::FastxReader(const std::string &file_name) {
  fp_ = file_name == "-" ? gzdopen(fileno(stdin), "r")
                         : gzopen(file_name.c_str(), "r");
  if (fp_ == nullptr) {
    throw std::invalid_argument("Cannot open file " + file_name);
  }
  kseq_reader_ = kseq_init(fp_);
  assert(kseq_reader_ != nullptr);
}

FastxReader::~FastxReader() {
  if (kseq_reader_) {
    kseq_destroy(kseq_reader_);
  }
  if (fp_) {
    gzclose(fp_);
  }
}

int64_t FastxReader::Read(SeqPackage *pkg, int64_t max_num,
                          int64_t max_num_bases, bool reverse) {
  int64_t num_bases = 0;
  for (int64_t i = 0; i < max_num; ++i) {
    auto record = ReadNext();
    if (record) {
      int b = 0, e = record->seq.l;
      if (trim_n_) {
        TrimN(record->seq.s, record->seq.l, &b, &e);
      }

      if (reverse) {
        pkg->AppendReversedStringSequence(record->seq.s + b, e - b);
      } else {
        pkg->AppendStringSequence(record->seq.s + b, e - b);
      }

      num_bases += e - b;
      if (num_bases >= max_num_bases && i % 2 == 1) {
        return i + 1;
      }
    } else {
      return i;
    }
  }
  return max_num;
}

void FastxReader::TrimN(const char *s, int len, int *out_bpos, int *out_epos) {
  *out_bpos = *out_epos = len;
  int i;
  for (i = 0; i < len; ++i) {
    if (s[i] == 'N' || s[i] == 'n') {
      if (*out_bpos < len) {
        break;
      }
    } else {
      if (*out_bpos == len) {
        *out_bpos = i;
      }
    }
  }
  *out_epos = i;
}
