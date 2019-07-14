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

#ifndef MEGAHIT_KMER_COUNTER_H
#define MEGAHIT_KMER_COUNTER_H

#include <stdint.h>
#include <mutex>
#include <string>
#include <vector>
#include "base_engine.h"
#include "definitions.h"
#include "edge_counter.h"
#include "sequence/io/edge/edge_writer.h"
#include "sequence/sequence_package.h"
#include "utils/atomic_wrapper.h"

struct KmerCounterOption {
  unsigned k{21};
  int solid_threshold{2};
  double host_mem{0};
  int n_threads{0};
  std::string read_lib_file{};
  std::string output_prefix{"out"};
  int mem_flag{1};
};

class KmerCounter : public BaseSequenceSortingEngine {
 public:
  static const unsigned kSentinelValue = 4;
  static const uint32_t kSentinelOffset = 4294967295U;
  explicit KmerCounter(const KmerCounterOption &opt)
      : BaseSequenceSortingEngine(opt.host_mem, opt.mem_flag, opt.n_threads),
        opt_(opt) {}
  ~KmerCounter() final = default;

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
  void PackEdge(uint32_t *dest, uint32_t *item, int64_t counting);

 private:
  KmerCounterOption opt_;

  int words_per_edge_{};  // number of (32-bit) words needed to represent a
                          // (k+1)-mer
  int64_t words_per_substr_{};  // substrings to be sorted by GPU
  SeqPackage seq_pkg_;
  std::vector<AtomicWrapper<uint32_t>> first_0_out_;
  std::vector<AtomicWrapper<uint32_t>> last_0_in_;
  // stat
  EdgeMultiplicityRecorder edge_counter_;
  // output
  EdgeWriter edge_writer_;
};

#endif  // MEGAHIT_KMER_COUNTER_H