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

#ifndef CX1_KMER_COUNT_H__
#define CX1_KMER_COUNT_H__

#include <stdint.h>
#include <mutex>
#include <string>
#include <vector>
#include "cx1.h"
#include "definitions.h"
#include "sequence/lib_info.h"
#include "sequence/readers/edge_io.h"
#include "sequence/sequence_package.h"
#include "utils/atomic_wrapper.h"

struct count_opt_t {
  int kmer_k{21};
  int kmer_freq_threshold{2};
  double host_mem{0};
  int num_cpu_threads{0};
  std::string read_lib_file{};
  std::string assist_seq_file{};
  std::string output_prefix{"out"};
  int mem_flag{1};
  bool need_mercy{true};
};

static const int kBucketPrefixLength = 8;
static const int kNumBuckets = 65536;  // pow(4, 8)
static const int64_t kDefaultLv1ScanTime = 8;
static const int64_t kMaxLv1ScanTime = 64;
static const int kSentinelValue = 4;
static const uint32_t kSentinelOffset = 4294967295U;

struct CX1KmerCount : public CX1<int, kNumBuckets> {
 public:
  CX1KmerCount(int kmer_k, int kmer_freq_threshold, int64_t host_mem, int mem_flag,  int num_cpu_threads,
              const std::string &read_lib_file, const std::string &assist_seq_file, const std::string &output_prefix):
              kmer_k(kmer_k), kmer_freq_threshold(kmer_freq_threshold), host_mem(host_mem), mem_flag(mem_flag),
              num_cpu_threads(num_cpu_threads), read_lib_file(read_lib_file), assist_seq_file(assist_seq_file),
              output_prefix(output_prefix) {}
  ~CX1KmerCount() final = default;

  int64_t encode_lv1_diff_base_func_(int64_t, global_data_t &) override;
  void prepare_func_(global_data_t &) override;  // num_items_, num_cpu_threads_ and num_cpu_threads_ must be set here
  void lv0_calc_bucket_size_func_(ReadPartition *) override;
  void init_global_and_set_cx1_func_(global_data_t &) override;  // xxx set here
  void lv1_fill_offset_func_(ReadPartition *) override;
  void lv2_extract_substr_(unsigned bucket_from, unsigned bucket_to, global_data_t &g, uint32_t *substr_ptr) override;
  void output_(int64_t start_index, int64_t end_index, int thread_id, global_data_t &g, uint32_t *substrings) override;
  void post_proc_func_(global_data_t &) override;

 private:
  // input options
  int kmer_k;
  int kmer_freq_threshold;
  int64_t host_mem;
  int mem_flag;
  int num_cpu_threads;
  std::string read_lib_file;
  std::string assist_seq_file;
  std::string output_prefix;

  int words_per_edge{};           // number of (32-bit) words needed to represent a (k+1)-mer
  int64_t words_per_substring{};  // substrings to be sorted by GPU

  SeqPackage package;
  std::vector<AtomicWrapper<uint32_t>> first_0_out;
  std::vector<AtomicWrapper<uint32_t>> last_0_in;
  // stat
  std::vector<std::vector<int64_t>> thread_edge_counting;
  // output
  EdgeWriter edge_writer;
  void PackEdge(uint32_t *dest, uint32_t *item, int64_t counting);
};

#endif  // CX1_KMER_COUNT_H__