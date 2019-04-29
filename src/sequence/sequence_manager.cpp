/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong & L3 Bioinformatics Limited
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* contact: Dinghua Li <dhli@cs.hku.hk> */

#include "sequence_manager.h"

#include <assert.h>
#include <zlib.h>
#include <string>

#include "utils/utils.h"
#include "kmlib/kmbit.h"

void SequenceManager::set_file(const std::string &file_name) {
  assert(f_type != kMegahitEdges && f_type != kSortedEdges);
  assert(files_.size() == 0);
  assert(kseq_readers_.size() == 0);

  files_.resize(1);
  files_[0] = file_name == "-" ? gzdopen(fileno(stdin), "r") : gzopen(file_name.c_str(), "r");
  assert(files_[0] != NULL);

  if (f_type == kMegahitContigs) {
    kseq_readers_.resize(1);
    kseq_readers_[0] = kseq_init(files_[0]);
    assert(kseq_readers_[0] != NULL);
  }
}

void SequenceManager::set_edge_files(const std::string &file_prefix) {
  assert(f_type == kSortedEdges || f_type == kMegahitEdges);
  assert(files_.size() == 0);
  assert(kseq_readers_.size() == 0);

  assert(!edge_reader_inited_);
  edge_reader_.set_file_prefix(file_prefix);
  edge_reader_.read_info();
  edge_reader_.init_files();

  if (edge_reader_.is_unsorted()) {
    f_type = kMegahitEdges;
  } else {
    f_type = kSortedEdges;
  }

  edge_reader_inited_ = true;
}

int64_t SequenceManager::ReadShortReads(int64_t max_num,
                                        int64_t max_num_bases,
                                        bool append,
                                        bool reverse,
                                        bool trimN,
                                        std::string file_name) {

  assert(f_type == kBinaryReads);
  if (!append) {
    package_->clear();
  }

  max_num = (max_num + 1) / 2 * 2;
  int64_t num_bases = 0;
  uint32_t read_len;

  for (int64_t i = 0; i < max_num; ++i) {
    if (gzread(files_[0], &read_len, sizeof(read_len)) == 0) {
      return i;
    }

    int num_words = DivCeiling(read_len, 16);

    if (buf_.size() < (unsigned) num_words) {
      buf_.resize(num_words);
    }

    auto n_read = gzread(files_[0], &buf_[0], sizeof(uint32_t) * num_words);
    assert(static_cast<size_t>(n_read) == num_words * sizeof(uint32_t));

    if (!reverse) {
      package_->AppendCompactSequence(&buf_[0], read_len);
    } else {
      package_->AppendReversedCompactSequence(&buf_[0], read_len);
    }

    num_bases += read_len;

    if (read_len >= max_num_bases && i % 2 == 1) {
      return i + 1;
    }
  }

  return max_num;
}

int64_t SequenceManager::ReadEdgesWithFixedLen(int64_t max_num, bool append) {
  if (!append) {
    multi_->clear();
    package_->clear();
  }

  if (f_type == kMegahitEdges) {
    for (int64_t i = 0; i < max_num; ++i) {
      uint32_t *next_edge = edge_reader_.NextUnsortedEdge();

      if (next_edge == nullptr) {
        return i;
      }

      package_->AppendCompactSequence(next_edge, edge_reader_.kmer_size() + 1);
      multi_->push_back(next_edge[edge_reader_.words_per_edge() - 1] & kMaxMul);
    }

    return max_num;
  } else if (f_type == kSortedEdges) {
    for (int64_t i = 0; i < max_num; ++i) {
      uint32_t *next_edge = edge_reader_.NextSortedEdge();

      if (next_edge == NULL) {
        return i;
      }

      package_->AppendCompactSequence(next_edge, edge_reader_.kmer_size() + 1);
      multi_->push_back(next_edge[edge_reader_.words_per_edge() - 1] & kMaxMul);
    }

    return max_num;
  }

  assert(false);
}