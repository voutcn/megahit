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

#ifndef MEGAHIT_READ_TO_SDBG_H
#define MEGAHIT_READ_TO_SDBG_H

#include <cstdint>
#include <string>
#include <vector>
#include "base_engine.h"
#include "definitions.h"
#include "edge_counter.h"
#include "kmlib/kmbitvector.h"
#include "sdbg/sdbg_writer.h"
#include "sequence/io/sequence_lib.h"
#include "sequence/sequence_package.h"

struct Read2SdbgOption {
  unsigned k{21};
  int solid_threshold{2};
  double host_mem{0};
  int n_threads{0};
  std::string read_lib_file;
  std::string output_prefix{"out"};
  int mem_flag{1};
  bool need_mercy{false};
};

struct SeqPkgWithSolidMarker {
  SeqPackage package;
  AtomicBitVector is_solid;  // mark <read_id, offset> is solid
  int n_mercy_files;
};

class Read2SdbgS1 : public BaseSequenceSortingEngine {
 public:
  static const unsigned kSentinelValue = 4;
  static const unsigned kBWTCharNumBits = 3;

  Read2SdbgS1(const Read2SdbgOption &opt, SeqPkgWithSolidMarker *pkg)
      : BaseSequenceSortingEngine(opt.host_mem, opt.mem_flag, opt.n_threads),
        opt_(opt),
        seq_pkg_(pkg) {}

 public:
  MemoryStat Initialize() override;

 protected:
  int64_t Lv0EncodeDiffBase(int64_t) override;
  void Lv0CalcBucketSize(int64_t seq_from, int64_t seq_to,
                         std::array<int64_t, kNumBuckets> *out) override;
  void Lv1FillOffsets(OffsetFiller &filler, int64_t seq_from,
                      int64_t seq_to) override;
  void Lv2ExtractSubString(OffsetFetcher &fetcher,
                           SubstrPtr substr_ptr) override;
  void Lv2Postprocess(int64_t start_index, int64_t end_index, int thread,
                      uint32_t *substr_ptr) override;
  void Lv0Postprocess() override;

 private:
  Read2SdbgOption opt_;
  SeqPkgWithSolidMarker *seq_pkg_;

  int64_t words_per_substr_;  // substrings to be sorted by GPU
  // stat-stage1
  EdgeMultiplicityRecorder edge_counter_;
  // output-stage1
  std::vector<FILE *> mercy_files_;
};

class Read2SdbgS2 : public BaseSequenceSortingEngine {
 public:
  static const unsigned kSentinelValue = 4;
  static const unsigned kBWTCharNumBits = 3;

  Read2SdbgS2(const Read2SdbgOption &opt, SeqPkgWithSolidMarker *pkg)
      : BaseSequenceSortingEngine(opt.host_mem, opt.mem_flag, opt.n_threads),
        opt_(opt),
        seq_pkg_(pkg) {}

 public:
  MemoryStat Initialize() override;

 protected:
  int64_t Lv0EncodeDiffBase(int64_t) override;
  void Lv0CalcBucketSize(int64_t seq_from, int64_t seq_to,
                         std::array<int64_t, kNumBuckets> *out) override;
  void Lv1FillOffsets(OffsetFiller &filler, int64_t seq_from,
                      int64_t seq_to) override;
  void Lv2ExtractSubString(OffsetFetcher &fetcher,
                           SubstrPtr substr_ptr) override;
  void Lv2Postprocess(int64_t start_index, int64_t end_index, int thread_id,
                      uint32_t *substr_ptr) override;
  void Lv0Postprocess() override;

 private:
  Read2SdbgOption opt_;
  SeqPkgWithSolidMarker *seq_pkg_;

  int64_t words_per_substr_;
  int words_per_dummy_node_;

  SdbgWriter sdbg_writer_;
};

#endif  // MEGAHIT_READ_TO_SDBG_H