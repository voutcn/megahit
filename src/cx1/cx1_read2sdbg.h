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

#ifndef CX1_READ2SDBG_H__
#define CX1_READ2SDBG_H__

#include <stdint.h>
#include <mutex>
#include <string>
#include <vector>
#include "cx1.h"
#include "definitions.h"
#include "kmlib/kmbitvector.h"
#include "sdbg/sdbg_writer.h"
#include "sequence/lib_info.h"
#include "sequence/sequence_package.h"

struct read2sdbg_opt_t {
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

  read2sdbg_opt_t() {
    kmer_k = 21;
    kmer_freq_threshold = 2;
    host_mem = 0;
    gpu_mem = 0;
    num_cpu_threads = 0;
    num_output_threads = 0;
    read_lib_file = "";
    output_prefix = "out";
    mem_flag = 1;
    need_mercy = false;
  }
};

namespace cx1_read2sdbg {

static const int kBucketPrefixLength = 8;
static const int kNumBuckets = 65536;  // pow(4, 8)
static const int64_t kDefaultLv1ScanTime = 8;
static const int64_t kMaxLv1ScanTime = 64;
static const int kSentinelValue = 4;
static const int64_t kMaxDummyEdges = 4294967294LL;
static const int kBWTCharNumBits = 3;

struct read2sdbg_global_t {
  CX1<read2sdbg_global_t, kNumBuckets> cx1;

  // input options
  int max_read_length;
  int kmer_k;
  int kmer_freq_threshold;
  int num_cpu_threads;
  int num_output_threads;
  int64_t host_mem;
  int mem_flag;
  bool need_mercy;
  std::string read_lib_file;
  std::string assist_seq_file;
  std::string output_prefix;

  int num_k1_per_read;  // max_read_len - kmer_k
  int num_mercy_files;
  int words_per_edge;           // number of (32-bit) words needed to represent a (k+1)-mer
  int64_t words_per_substring;  // substrings to be sorted by GPU
  int words_per_dummy_node;
  int offset_num_bits;  // the number of bits needed to store the offset of a base in the
                        // read/(k+1)-mer (i.e. log(read_length))
  int64_t max_bucket_size;
  int64_t tot_bucket_size;
  int words_per_read;  // number of (32-bit) words needed to represent a read in 2-bit-per-char
                       // format
  int read_length_mask;
  int64_t num_reads;             // total number of reads (including assist seq)
  int64_t num_short_read_bases;  // number of bases in short reads excl assist seq
  int64_t num_short_reads;       // total number of short reads
  // new sorting
  int64_t max_bucket_size_for_dynamic_sort;
  int64_t max_sorting_items;

  // big arrays
  SeqPackage package;
  std::vector<lib_info_t> lib_info;
  AtomicBitVector is_solid;  // mark <read_id, offset> is solid
  int32_t *lv1_items;        // each item is an offset (read ID and position) in differential representation
  std::mutex lv1_items_scanning_lock;

  // memory usage
  int64_t mem_packed_reads;

  // stat-stage1
  int64_t *edge_counting;  // count the number of (k+1)mer with occurs i times
  int64_t *thread_edge_counting;

  // output-stage1
  std::vector<FILE *> mercy_files;

  // output-stage2
  SdbgWriter sdbg_writer;
};

namespace s1 {
// stage1 cx1 core functions
int64_t s1_encode_lv1_diff_base(int64_t read_id, read2sdbg_global_t &g);
void s1_read_input_prepare(read2sdbg_global_t &g);  // num_items_, num_cpu_threads_ and
                                                    // num_output_threads_ must be set here
void *s1_lv0_calc_bucket_size(void *);              // pthread working function
void s1_init_global_and_set_cx1(read2sdbg_global_t &g);
void *s1_lv1_fill_offset(void *);  // pthread working function
void s1_lv1_direct_sort_and_count(read2sdbg_global_t &g);
void s1_post_proc(read2sdbg_global_t &g);
}  // namespace s1

namespace s2 {
// stage2 cx1 core functions
int64_t s2_encode_lv1_diff_base(int64_t read_id, read2sdbg_global_t &g);
void s2_read_mercy_prepare(read2sdbg_global_t &g);  // num_items_, num_cpu_threads_ and
                                                    // num_output_threads_ must be set here
void *s2_lv0_calc_bucket_size(void *);              // pthread working function
void s2_init_global_and_set_cx1(read2sdbg_global_t &g);
void *s2_lv1_fill_offset(void *);  // pthread working function
void s2_lv1_direct_sort_and_proc(read2sdbg_global_t &g);
void s2_post_proc(read2sdbg_global_t &g);
}  // namespace s2

}  // namespace cx1_read2sdbg

#endif  // CX1_READ2SDBG_H__