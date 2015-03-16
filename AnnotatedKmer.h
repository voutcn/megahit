#ifndef ANNOTATED_KMER_H__
#define ANNOTATED_KMER_H__

/**
 * @brief: A kmer struct composed by 64bit words
 *  	   the last bits in the last word are reserved for annotation (e.g. counting)
 * @tparam: kNumUint64 number of 64bit words
 * @tparam: kNumAnnBits number of bits reserved for annotation (< 64)
 */
template <int kNumUint64, int kNumAnnBits>
struct AnnotatedKmer
{
	uint64_t data_[kNumUint64];
	static const uint32_t kBitsForKmer = (kNumUint64 * 64 - kNumAnnBits);
	static const uint32_t kMaxSize = kBitsForKmer / 2;
	static const uint64_t kMaxAnn = (uint64_t(1) << kNumAnnBits) - 1;
	static const uint64_t kAnnMask = kMaxAnn << (64 - kNumAnnBits);
	static const uint64_t kLastWordMask = ~kMaxAnn;

	AnnotatedKmer()
	{ std::memset(data_, 0, sizeof(uint64_t) * kNumUint64); }

	AnnotatedKmer(const AnnotatedKmer &ann_kmer)
    { std::memcpy(data_, ann_kmer.data_, sizeof(uint64_t) * kNumUint64); }

    ~AnnotatedKmer() {}

    const AnnotatedKmer &operator = (const AnnotatedKmer &ann_kmer)
    { std::memcpy(data_, ann_kmer.data_, sizeof(uint64_t) * kNumUint64); return *this; }

    uint8_t operator [] (uint32_t index) const 
    { return (data_[index>>5] >> ((index & 31) << 1)) & 3; }

    uint8_t get_base(uint32_t index) const
    { return (data_[index>>5] >> ((index & 31) << 1)) & 3; }

    void set_base(uint32_t index, uint8_t ch)
    {
        ch &= 3;
        unsigned offset = (index & 31) << 1;
        data_[index>>5] = (data_[index>>5] & ~(3ULL << offset)) | (uint64_t(ch) << offset);
    }

    void ShiftAppend(uint8_t ch, int k) 
    { // 0 <= ch < 4; WARNING: the annotation will be changed
    	unsigned used_words = (k + 31) >> 5;

    	for (unsigned i = 0; i+1 < used_words; ++i)
            data_[i] = (data_[i] >> 2) | (data_[i+1] << 62);

        data_[used_words-1] = (data_[used_words-1] >> 2) | (uint64_t(ch) << (((k - 1) & 31) << 1));
    }

    void ShiftPreappend(uint8_t ch, int k)
    { // 0 <= ch < 4; WARNING: the annotation will be changed
        ch &= 3;
        unsigned used_words = (k + 31) >> 5;

        for (int i = used_words-1; i > 0; --i)
            data_[i] = (data_[i] << 2) | (data_[i-1] >> 62);
        data_[0] = (data_[0] << 2) | ch;

        if ((k & 31) != 0)
            data_[used_words-1] &= (1ULL << ((k & 31) << 1)) - 1;
    }

    bool KmerEqualTo(const AnnotatedKmer &ann_kmer) const
    {
    	for (unsigned i = 0; i+1 < kNumUint64; ++i)
    		if (data_[i] != ann_kmer[i])
    			return false;
    	if ((data_[kNumUint64-1] & kLastWordMask) != (ann_kmer.data_[kNumUint64-1] & kLastWordMask))
    		return false;

    	return true;
    }

    void UpdateAnn(int64_t x)
    {
    	data_[kNumUint64-1] = (data_[kNumUint64-1] & kLastWordMask) | (x << (64 - kNumAnnBits));
    }

    void AnnAtomicIncrement()
    { // i.e. if (ann <= kMaxAnn) { ann++; }
    	while (true) {
    		uint64_t old_value = data_[kNumUint64-1];
    		if ((old_value >> (64 - kNumAnnBits)) >= kMaxAnn) {
    			return;
    		}

    		uint64_t new_value = old_value + (1 << (64 - kNumAnnBits));
    		if (__sync_bool_compare_and_swap(data_[kNumUint64-1], old_value, new_value)) {
    			return;
    		}
    	}
    }
};

template <class T>
struct KmerEqualTo: binary_function <T,T,bool>
{
	bool operator(const T &x, const T &y) const {
		return x.KmerEqualTo(y);
	}
};

#endif # ANNOTATED_KMER_H__