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

#ifndef MEGAHIT_SEQ_TO_SDBG_H
#define MEGAHIT_SEQ_TO_SDBG_H

#include <sequence/sequence_package.h>
#include <stdint.h>
#include <mutex>
#include <string>
#include <vector>
#include "base_engine.h"
#include "sdbg/sdbg_writer.h"

struct Seq2SdbgOption {
  double host_mem{0};
  int n_threads{0};
  unsigned k{0};
  unsigned k_from{0};
  std::string contig;
  std::string bubble_seq;
  std::string addi_contig;
  std::string local_contig;
  std::string input_prefix;
  std::string output_prefix;
  int mem_flag{1};
  bool need_mercy{false};
};

class SeqToSdbg : public BaseSequenceSortingEngine {
 public:
  // binary search look up table
  static const unsigned kLookUpPrefixLength = 12;
  static const unsigned kLookUpShift = 32 - kLookUpPrefixLength * 2;
  static const unsigned kLookUpSize = 1u << (2 * kLookUpPrefixLength);
  static const unsigned kSentinelValue = 4;
  static const unsigned kBWTCharNumBits = 3;

  explicit SeqToSdbg(const Seq2SdbgOption &opt)
      : BaseSequenceSortingEngine(opt.host_mem, opt.mem_flag, opt.n_threads),
        opt_(opt) {}

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
  // input options
  Seq2SdbgOption opt_;
  int64_t words_per_substr_{};

  // big arrays
  SeqPackage seq_pkg_;
  std::vector<mul_t> multiplicity;

  // output
  SdbgWriter sdbg_writer_;
  void GenMercyEdges();
};
#endif  // MEGAHIT_SEQ_TO_SDBG_H