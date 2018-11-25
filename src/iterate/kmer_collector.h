//
// Created by vout on 11/24/18.
//

#ifndef MEGAHIT_SPANNING_KMER_COLLECTOR_H
#define MEGAHIT_SPANNING_KMER_COLLECTOR_H

#include <mutex>
#include <omp.h>
#include "sdbg/sdbg_def.h"
#include "sequence/kmer_plus.h"
#include "sparsepp/spp.h"
#include "edge_io.h"

template<class KmerType>
class KmerCollector {
 public:
  using kmer_type = KmerType;
  using kmer_plus = KmerPlus<KmerType, mul_t>;
  using hash_set = spp::sparse_hash_set<kmer_plus, KmerHash>;

  KmerCollector(unsigned k, const std::string &out_prefix, unsigned n_threads)
      : k_(k), output_prefix_(out_prefix), n_threads_(n_threads) {
    last_shift_ = k_ % 16;
    last_shift_ = (last_shift_ == 0 ? 0 : 16 - last_shift_) * 2;
    words_per_kmer_ = DivCeiling(k_ * 2 + kBitsPerMul, 32);
    buffer_.resize(n_threads * words_per_kmer_);

    writer_.set_num_threads(n_threads);
    writer_.set_file_prefix(out_prefix);
    writer_.set_unsorted();
    writer_.set_kmer_size(k_ - 1);
    writer_.init_files();
  }

  void Insert(const KmerType &kmer, mul_t mul) {
    std::lock_guard<std::mutex> lk(lock_);
    collection_.emplace(kmer, mul);
  }
  const hash_set &collection() const {
    return collection_;
  }
  void FlushToFile() {
    for (auto &item: collection_) {
      WriteToFile(item.kmer, item.aux);
    }
  }
 private:
  void WriteToFile(const KmerType &kmer, mul_t mul) {
    int tid = omp_get_thread_num();
    uint32_t *start_ptr = buffer_.data() + tid * words_per_kmer_;
    auto ptr = start_ptr;
    uint32_t w = 0;

    for (unsigned j = 0; j < k_; ++j) {
      w = (w << 2) | kmer.GetBase(k_ - 1 - j);
      if (j % 16 == 15) {
        *ptr = w;
        w = 0;
        ++ptr;
      }
    }

    *ptr = (w << last_shift_);
    assert((start_ptr[words_per_kmer_ - 1] & kMaxMul) == 0);
    start_ptr[words_per_kmer_ - 1] |= mul;
    writer_.write_unsorted(ptr, tid);
  }
 private:
  unsigned k_;
  std::string output_prefix_;
  unsigned n_threads_;
  std::mutex lock_;
  hash_set collection_;

  EdgeWriter writer_;
  unsigned last_shift_;
  unsigned words_per_kmer_;
  std::vector<uint32_t> buffer_;
};

#endif //MEGAHIT_SPANNING_KMER_COLLECTOR_H
