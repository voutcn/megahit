#include <utility>

#include <utility>

/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong & L3 Bioinformatics
 * Limited
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

#ifndef LIB_INFO_H
#define LIB_INFO_H

#include <cstdint>
#include <string>
#include "binary_reader.h"
#include "paired_fastx_reader.h"
#include "sequence/sequence_package.h"
#include "utils/utils.h"

class SequenceLib {
 public:
  SequenceLib(SeqPackage *data_holder, int64_t held_begin, int64_t held_end,
              unsigned max_read_length, bool is_paired, std::string description)
      : data_holder_(data_holder),
        held_begin_(held_begin),
        held_end_(held_end),
        max_read_len_(max_read_length),
        is_paired_(is_paired),
        description_(std::move(description)) {}

  SeqPackage::SeqView GetSequenceView(int64_t seq_id) const {
    assert(seq_id >= 0 && seq_id + held_begin_ < held_end_);
    return data_holder_->GetSeqView(seq_id + held_begin_);
  }

  unsigned GetMaxLength() const { return max_read_len_; }

  bool IsPaired() const { return is_paired_; }

  void DumpMetadata(std::ostream &os) const {
    os << description_ << '\n';
    os << held_begin_ << ' ' << held_end_ << ' ' << max_read_len_ << ' '
       << is_paired_ << '\n';
  }

  void LoadMetadata(std::istream &is) {
    is >> description_ >> held_begin_ >> held_end_ >> max_read_len_ >>
        is_paired_;
  }

  size_t seq_count() const { return held_end_ - held_begin_; }

 private:
  SeqPackage *data_holder_;
  int64_t held_begin_;
  int64_t held_end_;
  unsigned max_read_len_;
  bool is_paired_;
  std::string description_;
};

class SequenceLibCollection {
 public:
  SequenceLibCollection() = default;
  explicit SequenceLibCollection(const std::string &path) : path_(path) {}

  static void Build(const std::string &lib_file, const std::string &out_prefix);

  void SetPath(const std::string &path) { path_ = path; }
  std::pair<int64_t, int64_t> GetSize() const;
  void Read(SeqPackage *pkg, bool reverse_seq = false);
  size_t size() const { return libs_.size(); }
  const SequenceLib &GetLib(size_t id) const { return libs_[id]; }

 private:
  std::string path_;
  std::vector<SequenceLib> libs_;
};

#endif