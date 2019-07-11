//
// Created by vout on 4/8/18.
//

#ifndef MEGAHIT_KMER_PLUS_H
#define MEGAHIT_KMER_PLUS_H

#include "kmer.h"
#include "xxhash/xxh3.h"

/**
 * @brief a kmer plus any annotation
 */
template <class KmerType, class Auxiliary>
struct KmerPlus {
  using kmer_type = KmerType;
  using aux_type = Auxiliary;
  static const unsigned n_bytes = kmer_type::n_bytes;

  KmerPlus(const kmer_type &kmer = kmer_type(),
           const aux_type &ann = aux_type())
      : kmer(kmer), aux(ann) {}

  bool operator<(const KmerPlus &rhs) const { return kmer < rhs.kmer; }
  bool operator==(const KmerPlus &rhs) const { return kmer == rhs.kmer; }
  unsigned kth_byte(int k) const { return kmer.kth_byte(k); }
  kmer_type kmer;
  aux_type aux;
};

struct KmerHash {
  template <unsigned NumWords, class WordType>
  size_t operator()(const Kmer<NumWords, WordType> &kmer) const {
    return XXH3_64bits(static_cast<const void *>(kmer.data()),
                       sizeof(WordType) * NumWords);
  }
  template <class KmerType, class Auxiliary>
  size_t operator()(const KmerPlus<KmerType, Auxiliary> &kmer_plus) const {
    return operator()(kmer_plus.kmer);
  }
};

#endif  // MEGAHIT_KMER_PLUS_H