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

#ifndef CX1_SEQUENCES2SDBG_H__
#define CX1_SEQUENCES2SDBG_H__

#include <sequence/sequence_package.h>
#include <stdint.h>
#include <mutex>
#include <string>
#include <vector>
#include "cx1.h"
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

namespace cx1_seq2sdbg {

static const int kBucketPrefixLength = 8;  // less than 16 (chars per word)
static const int kBucketBase = 4;
static const int kNumBuckets = 65536;  // pow(4, 8)
// binary search look up table
static const int kLookUpPrefixLength = 12;
static const int kLookUpShift = 32 - kLookUpPrefixLength * 2;
static const int kLookUpSize = 1 << (2 * kLookUpPrefixLength);

static const int64_t kDefaultLv1ScanTime = 8;
static const int64_t kMaxLv1ScanTime = 64;
static const int kSentinelValue = 4;
static const int kBWTCharNumBits = 3;

struct seq2sdbg_global_t;

struct CX1Seq2Sdbg: public CX1<seq2sdbg_global_t, kNumBuckets> {
  int64_t encode_lv1_diff_base_func_(int64_t, global_data_t &) override;
  void prepare_func_(global_data_t &) override;  // num_items_, num_cpu_threads_ and num_cpu_threads_ must be set here
  void lv0_calc_bucket_size_func_(ReadPartition *) override;
  void init_global_and_set_cx1_func_(global_data_t &) override;  // xxx set here
  void lv1_fill_offset_func_(ReadPartition *) override;
  void lv2_extract_substr_(unsigned bucket_from, unsigned bucket_to, global_data_t &g, uint32_t *substr_ptr) override;
  void output_(int64_t start_index, int64_t end_index, int thread_id, global_data_t &g, uint32_t *substrings) override;
  void post_proc_func_(global_data_t &) override;
};

struct seq2sdbg_global_t {
  std::unique_ptr<CX1Seq2Sdbg> cx1;

  // input options
  int kmer_k;
  int kmer_from;
  int num_cpu_threads;
  int64_t host_mem;
  int mem_flag;
  bool need_mercy;

  std::string contig;
  std::string bubble_seq;
  std::string addi_contig;
  std::string local_contig;
  std::string input_prefix;
  std::string output_prefix;

  int64_t num_seq;
  int64_t words_per_substring;  // substrings to be sorted by GPU
  int64_t max_bucket_size;
  int64_t tot_bucket_size;
  int words_per_dummy_node;

  // big arrays
  SeqPackage package;
  std::vector<mul_t> multiplicity;

  int64_t max_sorting_items;
  // memory usage
  int64_t mem_packed_seq;

  // output
  SdbgWriter sdbg_writer;
};

}  // end of namespace cx1_seq2sdbg
#endif  // CX1_SEQUENCES2SDBG_H__