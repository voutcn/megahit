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

#ifndef MEGAHIT_SEQ_TO_SDBG_H
#define MEGAHIT_SEQ_TO_SDBG_H

#include <sequence/sequence_package.h>
#include <stdint.h>
#include <mutex>
#include <string>
#include <vector>
#include "base_sequence_sorting_engine.h"
#include "sdbg/sdbg_writer.h"

struct seq2sdbg_opt_t {
  double host_mem{0};
  int num_cpu_threads{0};
  int kmer_k{0};
  int kmer_from{0};
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
  static const int kBucketBase = 4;
// binary search look up table
  static const int kLookUpPrefixLength = 12;
  static const int kLookUpShift = 32 - kLookUpPrefixLength * 2;
  static const int kLookUpSize = 1 << (2 * kLookUpPrefixLength);
  static const int kSentinelValue = 4;
  static const int kBWTCharNumBits = 3;

  explicit SeqToSdbg(const seq2sdbg_opt_t &opt) :
      host_mem(opt.host_mem), mem_flag(opt.mem_flag), num_cpu_threads(opt.num_cpu_threads),
      kmer_k(opt.kmer_k), kmer_from(opt.kmer_from),
      input_prefix(opt.input_prefix), output_prefix(opt.output_prefix),
      contig(opt.contig), bubble_seq(opt.bubble_seq), addi_contig(opt.addi_contig),
      local_contig(opt.local_contig), need_mercy(opt.need_mercy) {}

  int64_t encode_lv1_diff_base_func_(int64_t) override;
  void prepare_func_() override;  // num_items_, num_cpu_threads_ and num_cpu_threads_ must be set here
  void lv0_calc_bucket_size_func_(ReadPartition *) override;
  void init_global_and_set_cx1_func_() override;  // xxx set here
  void lv1_fill_offset_func_(ReadPartition *) override;
  void lv2_extract_substr_(unsigned bucket_from, unsigned bucket_to, uint32_t *substr_ptr) override;
  void output_(int64_t start_index, int64_t end_index, int thread_id, uint32_t *substrings) override;
  void post_proc_func_() override;

 private:
  // input options
  int64_t host_mem;
  int mem_flag;
  int num_cpu_threads;

  int kmer_k;
  int kmer_from;

  std::string input_prefix;
  std::string output_prefix;
  std::string contig;
  std::string bubble_seq;
  std::string addi_contig;
  std::string local_contig;
  bool need_mercy;

  int64_t words_per_substring{};
  int words_per_dummy_node{};

  // big arrays
  SeqPackage package;
  std::vector<mul_t> multiplicity;

  // output
  SdbgWriter sdbg_writer;
  void GenMercyEdges();
};
#endif  // MEGAHIT_SEQ_TO_SDBG_H