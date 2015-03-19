#ifndef KMER_PLUS_H__
#define KMER_PLUS_H__

#include <algorithm>
#include "kmer.h"

template <uint32_t kNumUint64, typename kmer_word_t, typename ann_t>
struct KmerPlus
{
	typedef Kmer<kNumUint64, kmer_word_t> kmer_t;
	kmer_t kmer;
	ann_t ann;

	explicit KmerPlus(const kmer_t &kmer = kmer_t(), const ann_t &ann = ann_t()): kmer(kmer), ann(ann) {}

    KmerPlus(const KmerPlus &rhs): kmer(rhs.kmer), ann(rhs.ann) {}

    const KmerPlus &operator =(const KmerPlus &rhs)
    {
    	kmer = rhs.kmer;
    	ann = rhs.ann;
    	return *this;
    }

    const kmer_t &key() const { return kmer; }
    void swap(KmerPlus &rhs)
    {
    	if (this != &rhs)
    	{
    		kmer.swap(rhs.kmer);
    		std::swap(ann, rhs.ann);
    	}
    }
} __attribute__((packed));

#endif // KMER_PLUS_H__