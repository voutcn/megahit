//
// Created by vout on 5/11/19.
//

#ifndef MEGAHIT_BINARY_READER_H
#define MEGAHIT_BINARY_READER_H

#include <fstream>
#include "base_reader.h"

#include "utils/buffered_reader.h"

class BinaryReader : public BaseSequenceReader {
 public:
  explicit BinaryReader(const std::string &filename)
      : is_(filename), buf_(120) {
    if (is_.bad()) {
      throw std::invalid_argument("Failed to open file " + filename);
    }
    reader_.reset(&is_);
  }

  int64_t Read(SeqPackage *pkg, int64_t max_num, int64_t max_num_bases,
               bool reverse) override {
    max_num = (max_num + 1) / 2 * 2;
    int64_t num_bases = 0;
    uint32_t read_len;

    for (int64_t i = 0; i < max_num; ++i) {
      if (reader_.read(&read_len) == 0) {
        return i;
      }
      auto num_words = DivCeiling(read_len, SeqPackage::kBasesPerWord);
      if (buf_.size() < num_words) {
        buf_.resize(num_words);
      }
      auto bytes_read = reader_.read(buf_.data(), num_words);
      assert(bytes_read == num_words * sizeof(buf_[0]));
      (void)(bytes_read);

      if (!reverse) {
        pkg->AppendCompactSequence(buf_.data(), read_len);
      } else {
        pkg->AppendReversedCompactSequence(buf_.data(), read_len);
      }

      num_bases += read_len;
      if (read_len >= max_num_bases && i % 2 == 1) {
        return i + 1;
      }
    }
    return max_num;
  }

 private:
  std::ifstream is_;
  BufferedReader reader_;
  std::vector<SeqPackage::TWord> buf_;
};

#endif  // MEGAHIT_BINARY_READER_H
