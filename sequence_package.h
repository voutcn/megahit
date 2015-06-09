#ifndef SEQUENCE_PACKAGE_H__
#define SEQUENCE_PACKAGE_H__

#include <stdint.h>
#include <vector>

/**
 * @brief hold a set of sequences
 */

struct SequencePackage {

	typedef uint32_t word_t;	// do not change
	const static unsigned kBitsPerWord = 8 * sizeof(word_t);
	const static unsigned kCharsPerWord = kBitsPerWord / 2;
	char dna_map_[256];

	std::vector<word_t> packed_seq; // packed all 
	std::vector<uint64_t> start_idx; // the index of the starting position of a sequence

	uint8_t unused_bits_; // the number of unused bits in the last word
	int max_read_len_;

	SequencePackage() {
		start_idx.push_back(0);
		packed_seq.push_back(word_t(0));
		unused_bits_ = kBitsPerWord;
		max_read_len_ = 0;
		for (int i = 0; i < 10; ++i) {
			dna_map_[(int)("ACGTNacgtn"[i])] = "0123201232"[i] - '0';
		}
	}

	~SequencePackage() {}

	void clear() {
		packed_seq.clear();
		packed_seq.push_back(word_t(0));
		start_idx.clear();
		start_idx.push_back(word_t(0));
		unused_bits_ = kBitsPerWord;
	}

	size_t size() {
		return start_idx.size() - 1;
	}

	size_t base_size() {
		return start_idx.back();
	}

	size_t size_in_byte() {
		return sizeof(uint32_t) * packed_seq.size() + sizeof(uint64_t) * start_idx.size();
	}

	size_t max_read_len() {
		return max_read_len_;
	}

	void shrink_to_fit() {
		// packed_seq.shrink_to_fit();
		// start_idx.shrink_to_fit();
	}

	size_t length(size_t seq_id) {
		return start_idx[seq_id + 1] - start_idx[seq_id];
	}

	uint8_t get_base(size_t seq_id, size_t offset) {
		uint64_t where = start_idx[seq_id] + offset;
		return packed_seq[where / kCharsPerWord] >> (kCharsPerWord - 1 - where % kCharsPerWord) * 2 & 3;
	}

	void AppendSeq(const char *s, int len) {
		for (int i = 0; i < len; ++i) {
			unused_bits_ -= 2;
			packed_seq.back() |= dna_map_[(int)s[i]] << unused_bits_;
			if (unused_bits_ == 0) {
				unused_bits_ = kBitsPerWord;
				packed_seq.push_back(word_t(0));
			}
		}
		if (len > max_read_len_) { max_read_len_ = len; }
		uint64_t end = start_idx.back() + len;
		start_idx.push_back(end);
	}

	void AppendReverseSeq(const char *s, int len) {
		for (int i = len - 1; i >= 0; --i) {
			unused_bits_ -= 2;
			packed_seq.back() |= dna_map_[(int)s[i]] << unused_bits_;
			if (unused_bits_ == 0) {
				unused_bits_ = kBitsPerWord;
				packed_seq.push_back(word_t(0));
			}
		}
		if (len > max_read_len_) { max_read_len_ = len; }
		uint64_t end = start_idx.back() + len;
		start_idx.push_back(end);
	}

	void AppendSeq(const word_t *s, int len) {
		if (len > max_read_len_) { max_read_len_ = len; }
		uint64_t end = start_idx.back() + len;
		start_idx.push_back(end);

		if (len * 2 <= unused_bits_) {
			unused_bits_ -= len * 2;
			packed_seq.back() |= s[0] >> (kCharsPerWord - len) * 2 << unused_bits_;
			if (unused_bits_ == 0) {
				unused_bits_ = kBitsPerWord;
				packed_seq.push_back(word_t(0));
			}
		} else {
			int num_words = (len + kCharsPerWord - 1) / kCharsPerWord;

			// append to unused bits
			packed_seq.back() |= s[0] >> (kBitsPerWord - unused_bits_);

			if (unused_bits_ == kBitsPerWord) {
				packed_seq.back() = s[0];
				for (int i = 1; i < num_words; ++i) {
					packed_seq.push_back(s[i]);
				}
			} else {
				int num_words_to_append = (len - unused_bits_ / 2 + kCharsPerWord - 1) / kCharsPerWord;
				for (int i = 0; i < num_words_to_append; ++i) {
					if (i + 1 < num_words)
						packed_seq.push_back((s[i] << unused_bits_) | (s[i+1] >> (kBitsPerWord - unused_bits_)));
					else
						packed_seq.push_back(s[i] << unused_bits_);
				}
			}

			int bits_in_last_word = end % kCharsPerWord * 2;
			unused_bits_ = kBitsPerWord - bits_in_last_word;
			if (unused_bits_ == kBitsPerWord) {
				packed_seq.push_back(word_t(0));
			} else {
				packed_seq.back() >>= unused_bits_;
				packed_seq.back() <<= unused_bits_;
			}
		}
	}

	void get_seq(std::vector<word_t> &s, size_t seq_id, int begin = 0, int end = -1) {
		if (end == -1) {
			end = length(seq_id) - 1;
		}

		size_t first_word = (start_idx[seq_id] + begin) / kCharsPerWord;
		size_t last_word = (start_idx[seq_id] + end) / kCharsPerWord;
		int first_shift = (start_idx[seq_id] + begin) % kCharsPerWord * 2;

		s.clear();
		if (end < begin) {
			return;
		}

		if (first_shift == 0) {
			for (size_t i = first_word; i <= last_word; ++i) {
				s.push_back(packed_seq[i]);
			}
		} else {
			for (size_t i = first_word; i < last_word; ++i) {
				s.push_back((packed_seq[i] << first_shift) | (packed_seq[i+1] >> (kBitsPerWord - first_shift)));
			}
			if (kCharsPerWord * s.size() < (unsigned)end - begin + 1) {
				s.push_back(packed_seq[last_word] << first_shift);
			}
		}

		unsigned shift_clean = kBitsPerWord - (end - begin + 1) * 2 % kBitsPerWord;
		if (shift_clean != kBitsPerWord) {
			s.back() >>= shift_clean;
			s.back() <<= shift_clean;	
		}

		// for (int j = 0; j < end - begin + 1; ++j) {
		// 	if ((s[j/16] >> (15-j%16) * 2 & 3) != get_base(seq_id, begin + j)) {
		// 		printf("j: %d, first_shift: %d\n", j, first_shift);
		// 	}
		// }
	}
};

#endif