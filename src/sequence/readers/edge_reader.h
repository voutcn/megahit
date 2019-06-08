//
// Created by vout on 5/11/19.
//

#ifndef MEGAHIT_EDGE_READER_H
#define MEGAHIT_EDGE_READER_H

#include <sequence/readers/edge_io.h>
#include "base_reader.h"

class EdgeReader : public BaseSequenceReader {
 public:
  EdgeReader(const std::string &file_prefix) {
    edge_reader_.SetFilePrefix(file_prefix);
    edge_reader_.ReadInfo();
    edge_reader_.InitFiles();
  }

  int64_t ReadUnsorted(SeqPackage *pkg, std::vector<mul_t> *mul, int64_t max_num) {
    for (int64_t i = 0; i < max_num; ++i) {
      uint32_t *next_edge = edge_reader_.NextUnsortedEdge();
      if (next_edge == nullptr) {
        return i;
      }
      pkg->AppendCompactSequence(next_edge, edge_reader_.k() + 1);
      if (mul) {
        mul->push_back(next_edge[edge_reader_.words_per_edge() - 1] & kMaxMul);
      }
    }
    return max_num;
  }

  int64_t ReadSorted(SeqPackage *pkg, std::vector<mul_t> *mul, int64_t max_num) {
    for (int64_t i = 0; i < max_num; ++i) {
      uint32_t *next_edge = edge_reader_.NextSortedEdge();
      if (next_edge == nullptr) {
        return i;
      }
      pkg->AppendCompactSequence(next_edge, edge_reader_.k() + 1);
      if (mul) {
        mul->push_back(next_edge[edge_reader_.words_per_edge() - 1] & kMaxMul);
      }
    }
    return max_num;
  }

  int64_t Read(SeqPackage *pkg, int64_t max_num, int64_t max_num_bases, bool reverse = false) override {
    if (edge_reader_.IsUnsorted()) {
      return ReadSorted(pkg, mul_, max_num);
    } else {
      return ReadUnsorted(pkg, mul_, max_num);
    }
  }

 private:
  MegahitEdgeReader edge_reader_;
  std::vector<mul_t> *mul_{nullptr};
};

#endif  // MEGAHIT_EDGE_READER_H
