//
// Created by vout on 4/28/19.
//

#ifndef MEGAHIT_FASTX_READER_H
#define MEGAHIT_FASTX_READER_H

#include <vector>
#include <string>
#include "kseq.h"
#include <zlib.h>
#include <assert.h>
#include <stdexcept>
#include <sequence/sequence_package.h>
#include "definitions.h"

#ifndef KSEQ_INITED
#define KSEQ_INITED
KSEQ_INIT(gzFile, gzread)
#endif

class FastxReader {
 public:
  FastxReader(const std::vector<std::string> &file_names);

  virtual ~FastxReader() {
    Clear();
  }

  void Clear();
  virtual int64_t Read(SeqPackage *pkg, int64_t max_num, int64_t max_num_bases, bool reverse, bool trim_n);

 protected:

  kseq_t *ReadNext() {
    while (cur_reader_ && kseq_read(cur_reader_) < 0) {
      ++file_idx_;
      if (file_idx_ == kseq_readers_.size()) {
        cur_reader_ = nullptr;
      } else {
        cur_reader_ = kseq_readers_[file_idx_];
      }
    }
    return cur_reader_;
  }

  virtual void ParseFastxComment(const char *comment) {}
  static void TrimN(const char *s, int len, int *out_bpos, int *out_epos);

  unsigned file_idx_{0};
  kseq_t *cur_reader_{nullptr};
  std::vector<gzFile> fps_;
  std::vector<kseq_t *> kseq_readers_;
};

class PairEndFastxReader : public FastxReader {
 public:
  explicit PairEndFastxReader(const std::vector<std::string> &file_names);

  int64_t Read(SeqPackage *pkg,
               int64_t max_num,
               int64_t max_num_bases,
               bool reverse,
               bool trim_n) override;

 protected:
  std::pair<kseq_t *, kseq_t *> ReadNextPair() {
    while (cur_readers_[0] && cur_readers_[1]) {
      auto r1 = kseq_read(cur_readers_[0]);
      auto r2 = kseq_read(cur_readers_[1]);
      if (r1 < 0 && r2 < 0) {
        file_idx_ += 2;
        if (file_idx_ == kseq_readers_.size()) {
          cur_readers_[0] = cur_readers_[1] = nullptr;
        } else {
          cur_readers_[0] = kseq_readers_[file_idx_];
          cur_readers_[1] = kseq_readers_[file_idx_ + 1];
        }
      } else {
        if (r1 < 0 || r2 < 0) {
          xfatal("Number of PE reads not match");
        }
        assert(r1 >= 0 && r2 >= 0);
        break;
      }
    }
    return std::make_pair(cur_readers_[0], cur_readers_[1]);
  }

 private:
  kseq_t *cur_readers_[2]{nullptr, nullptr};
};

#endif //MEGAHIT_FASTX_READER_H
