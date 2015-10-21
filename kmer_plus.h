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

#ifndef KMER_PLUS_H__
#define KMER_PLUS_H__

#include <algorithm>
#include "kmer.h"

template <uint32_t kNumWords, typename kmer_word_t, typename ann_t>
struct KmerPlus {
    typedef Kmer<kNumWords, kmer_word_t> kmer_t;
    kmer_t kmer;
    ann_t ann;

    explicit KmerPlus(const kmer_t &kmer = kmer_t(), const ann_t &ann = ann_t()): kmer(kmer), ann(ann) {}

    KmerPlus(const KmerPlus &rhs): kmer(rhs.kmer), ann(rhs.ann) {}

    const KmerPlus &operator =(const KmerPlus &rhs) {
        kmer = rhs.kmer;
        ann = rhs.ann;
        return *this;
    }

    const kmer_t &key() const {
        return kmer;
    }
    void swap(KmerPlus &rhs) {
        if (this != &rhs) {
            kmer.swap(rhs.kmer);
            std::swap(ann, rhs.ann);
        }
    }
} __attribute__((packed));

#endif // KMER_PLUS_H__