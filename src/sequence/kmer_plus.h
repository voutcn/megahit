//
// Created by vout on 4/8/18.
//

#ifndef MEGAHIT_KMER_PLUS_H
#define MEGAHIT_KMER_PLUS_H

#include "kmer.h"

/**
 * @brief a kmer plus any annotation
 */
template<unsigned NumWords, class WordType, class Auxiliary>
struct KmerPlus {
  using KmerType = Kmer<NumWords, WordType>;
  using AuxType = Auxiliary;
  static const unsigned n_bytes = KmerType::n_bytes;

  KmerPlus(const KmerType &kmer = KmerType(), const AuxType &ann = AuxType())
      : kmer(kmer), aux(ann) {}

  bool operator<(const KmerPlus &rhs) const {
    return kmer < rhs.kmer;
  }
  bool operator==(const KmerPlus &rhs) const {
    return kmer == rhs.kmer;
  }
  unsigned kth_byte(int k) const {
    return kmer.kth_byte(k);
  }
  KmerType kmer;
  AuxType aux;
};

#include "xxHash/xxhash.h"

struct KmerHash
{
  template <unsigned NumWords, class WordType>
  size_t operator()(Kmer<NumWords, WordType> const &kmer) const {
    return XXH64(static_cast<const void*>(kmer.data()), sizeof(WordType) * NumWords, 0);
  }
  template <unsigned NumWords, class WordType, class Auxiliary>
  size_t operator()(KmerPlus<NumWords, WordType, Auxiliary> const &kmer_plus) const {
    return operator()(kmer_plus.kmer);
  }
};

#endif //MEGAHIT_KMER_PLUS_H