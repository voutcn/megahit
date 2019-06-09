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
  int kmer_k;
  int kmer_freq_threshold;
  double host_mem;
  double gpu_mem;
  int num_cpu_threads;
  int num_output_threads;
  std::string read_lib_file;
  std::string assist_seq_file;
  std::string output_prefix;
  int mem_flag;
  bool need_mercy;

  count_opt_t() {
    kmer_k = 21;
    kmer_freq_threshold = 2;
    host_mem = 0;
    gpu_mem = 0;
    num_cpu_threads = 0;
    num_output_threads = 0;
    read_lib_file = "";
    output_prefix = "out";
    mem_flag = 1;
    need_mercy = true;
  }
};

namespace cx1_kmer_count {

static const int kBucketPrefixLength = 8;
static const int kNumBuckets = 65536;  // pow(4, 8)
static const int64_t kDefaultLv1ScanTime = 8;
static const int64_t kMaxLv1ScanTime = 64;
static const int kSentinelValue = 4;
static const uint32_t kSentinelOffset = 4294967295U;

struct count_global_t;

struct CX1KmerCount : public CX1<count_global_t, kNumBuckets> {
 public:
  int64_t encode_lv1_diff_base_func_(int64_t, global_data_t &) override;
  void prepare_func_(global_data_t &) override;  // num_items_, num_cpu_threads_ and num_output_threads_ must be set here
  void lv0_calc_bucket_size_func_(void *) override;
  void init_global_and_set_cx1_func_(global_data_t &) override;  // xxx set here
  void lv1_fill_offset_func_(void *) override;
  void lv1_sort_and_proc(global_data_t &) override;
  void post_proc_func_(global_data_t &) override;
};

struct count_global_t {
  CX1KmerCount cx1;

  // input options
  int max_read_length;
  int kmer_k;
  int kmer_freq_threshold;
  int num_cpu_threads;
  int num_output_threads;
  int64_t host_mem;
  int mem_flag;
  std::string read_lib_file;
  std::string assist_seq_file;
  std::string output_prefix;

  int words_per_edge;           // number of (32-bit) words needed to represent a (k+1)-mer
  int64_t words_per_substring;  // substrings to be sorted by GPU
  int64_t max_bucket_size;
  int64_t tot_bucket_size;
  int64_t num_reads;  // total number of reads

  // big arrays
  SeqPackage package;
  std::vector<lib_info_t> lib_info;

  // lv1 new sorting scheme
  int64_t max_sorting_items;

  std::vector<AtomicWrapper<uint32_t>> first_0_out;
  std::vector<AtomicWrapper<uint32_t>> last_0_in;
  std::vector<int32_t> lv1_items;
  std::mutex lv1_items_scanning_lock;

  // memory usage
  int64_t mem_packed_reads;

  // stat
  std::vector<std::vector<int64_t>> thread_edge_counting;

  // output
  EdgeWriter edge_writer;
};

}  // end of namespace cx1_kmer_count
#endif  // CX1_KMER_COUNT_H__