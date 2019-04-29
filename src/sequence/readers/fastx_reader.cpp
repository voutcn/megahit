//
// Created by vout on 4/28/19.
//

#include "fastx_reader.h"

FastxReader::FastxReader(const std::vector<std::string> &file_names) {
  for (const auto &file_name : file_names) {
    auto fp = file_name == "-" ? gzdopen(fileno(stdin), "r") : gzopen(file_name.c_str(), "r");
    if (fp == nullptr) {
      Clear();
      throw std::invalid_argument("Cannot open file " + file_name);
    }
    auto kseq_reader = kseq_init(fp);
    assert(kseq_reader != nullptr);
    fps_.push_back(fp);
    kseq_readers_.push_back(kseq_reader);
  }
  if (!kseq_readers_.empty()) {
    cur_reader_ = kseq_readers_[0];
  }
}

void FastxReader::Clear() {
  for (auto reader: kseq_readers_) {
    kseq_destroy(reader);
  }
  kseq_readers_.clear();
  for (auto fp : fps_) {
    gzclose(fp);
  }
  fps_.clear();
  file_idx_ = 0;
}

int64_t FastxReader::Read(SeqPackage *pkg, int64_t max_num, int64_t max_num_bases, bool reverse, bool trim_n) {
  int64_t num_bases = 0;
  for (int64_t i = 0; i < max_num; ++i) {
    auto record = ReadNext();
    if (record) {
      ParseFastxComment(record->comment.s);

      int b = 0, e = record->seq.l;
      if (trim_n) {
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

PairEndFastxReader::PairEndFastxReader(const std::vector<std::string> &file_names) : FastxReader(file_names) {
  if (kseq_readers_.size() % 2 != 0) {
    Clear();
    throw std::invalid_argument("Files must be paired with PairEndFastxReader");
  }
  if (!kseq_readers_.empty()) {
    cur_readers_[0] = kseq_readers_[0];
    cur_readers_[1] = kseq_readers_[1];
  }
}

int64_t PairEndFastxReader::Read(SeqPackage *pkg, int64_t max_num, int64_t max_num_bases, bool reverse, bool trim_n) {
  int64_t num_bases = 0;
  for (int64_t i = 0; i < max_num; i += 2) {
    auto record = ReadNextPair();
    auto r0 = record.first;
    auto r1 = record.second;
    if (r0 && r1) {
      int b0 = 0, e0 = r0->seq.l;
      int b1 = 0, e1 = r1->seq.l;

      if (trim_n) {
        TrimN(r0->seq.s, r0->seq.l, &b0, &e0);
        TrimN(r1->seq.s, r1->seq.l, &b1, &e1);
      }

      if (reverse) {
        pkg->AppendReversedStringSequence(r0->seq.s + b0, e0 - b0);
        pkg->AppendReversedStringSequence(r1->seq.s + b1, e1 - b1);
      } else {
        pkg->AppendStringSequence(r0->seq.s + b0, e0 - b0);
        pkg->AppendStringSequence(r1->seq.s + b1, e1 - b1);
      }

      num_bases += e0 - b0 + e1 - b1;

      if (num_bases >= max_num_bases) {
        return i + 2;
      }
    } else {
      return i;
    }
  }
  return max_num;
}
