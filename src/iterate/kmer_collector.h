//
// Created by vout on 11/24/18.
//

#ifndef MEGAHIT_SPANNING_KMER_COLLECTOR_H
#define MEGAHIT_SPANNING_KMER_COLLECTOR_H

#include <mutex>
#include "sdbg/sdbg_def.h"
#include "sequence/kmer_plus.h"
#include "sparsepp/spp.h"

template<class KmerType>
class KmerCollector {
 public:
  using kmer_type = KmerType;
  using kmer_plus = KmerPlus<KmerType::kNumWords, typename KmerType::word_type, mul_t>;
  using hash_set = spp::sparse_hash_set<kmer_plus, KmerHash>;
  void insert(const KmerType &kmer, mul_t mul) {
    std::lock_guard<std::mutex> lk(lock_);
    collection_.emplace(kmer, mul);
  }
  const hash_set &collection() const {
    return collection_;
  }
 private:
  std::mutex lock_;
  hash_set collection_;
};

#endif //MEGAHIT_SPANNING_KMER_COLLECTOR_H
