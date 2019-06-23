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

#ifndef MEGAHIT_READ_TO_SDBG_H
#define MEGAHIT_READ_TO_SDBG_H

#include <stdint.h>
#include <mutex>
#include <string>
#include <vector>
#include "base_sequence_sorting_engine.h"
#include "definitions.h"
#include "kmlib/kmbitvector.h"
#include "sdbg/sdbg_writer.h"
#include "sequence/lib_info.h"
#include "sequence/sequence_package.h"

struct read2sdbg_opt_t {
  int kmer_k{21};
  int kmer_freq_threshold{2};
  double host_mem{0};
  int num_cpu_threads{0};
  std::string read_lib_file;
  std::string assist_seq_file;
  std::string output_prefix{"out"};
  int mem_flag{1};
  bool need_mercy{false};
};

struct SeqPkgWithSolidMarker {
  SeqPackage package;
  AtomicBitVector is_solid;  // mark <read_id, offset> is solid
  int n_mercy_files;
};

class ReadToSdbg : public BaseSequenceSortingEngine {
 public:
  static const int kSentinelValue = 4;
  static const int64_t kMaxDummyEdges = 4294967294LL;
  static const int kBWTCharNumBits = 3;
};

class Read2SdbgS1 : public ReadToSdbg {
 public:

  Read2SdbgS1(const read2sdbg_opt_t &opt, SeqPkgWithSolidMarker *pkg)
      : opt(opt), seq_pkg_(pkg) {}

  int64_t encode_lv1_diff_base_func_(int64_t) override;
  void prepare_func_() override;  // num_items_, num_cpu_threads_ and num_cpu_threads_ must be set here
  void lv0_calc_bucket_size_func_(ReadPartition *) override;
  void init_global_and_set_cx1_func_() override;  // xxx set here
  void lv1_fill_offset_func_(ReadPartition *) override;
  void lv2_extract_substr_(unsigned bucket_from, unsigned bucket_to, uint32_t *substr_ptr) override;
  void output_(int64_t start_index, int64_t end_index, int thread_id, uint32_t *substrings) override;
  void post_proc_func_() override;
 private:
  read2sdbg_opt_t opt;
  SeqPkgWithSolidMarker *seq_pkg_;

  int64_t words_per_substring;  // substrings to be sorted by GPU
  // stat-stage1
  std::vector<std::vector<int64_t>> thread_edge_counting;
  // output-stage1
  std::vector<FILE *> mercy_files;
};

struct Read2SdbgS2 : public ReadToSdbg {
  Read2SdbgS2(const read2sdbg_opt_t &opt, SeqPkgWithSolidMarker *pkg)
      :opt(opt), seq_pkg_(pkg) {}
  int64_t encode_lv1_diff_base_func_(int64_t) override;
  void prepare_func_() override;  // num_items_, num_cpu_threads_ and num_cpu_threads_ must be set here
  void lv0_calc_bucket_size_func_(ReadPartition *) override;
  void init_global_and_set_cx1_func_() override;  // xxx set here
  void lv1_fill_offset_func_(ReadPartition *) override;
  void lv2_extract_substr_(unsigned bucket_from, unsigned bucket_to, uint32_t *substr_ptr) override;
  void output_(int64_t start_index, int64_t end_index, int thread_id, uint32_t *substrings) override;
  void post_proc_func_() override;
 private:
  read2sdbg_opt_t opt;
  SeqPkgWithSolidMarker *seq_pkg_;

  int64_t words_per_substring;  // substrings to be sorted by GPU
  int words_per_dummy_node;

  SdbgWriter sdbg_writer;
};

#endif  // MEGAHIT_READ_TO_SDBG_H