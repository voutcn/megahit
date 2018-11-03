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

#endif //MEGAHIT_KMER_PLUS_H