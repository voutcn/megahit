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
	uint32_t max_read_len_;

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

	size_t max_read_len() {
		return max_read_len_;
	}

	// void shrink_to_fit() {
	// 	packed_seq.shrink_to_fit();
	// 	start_idx.shrink_to_fit();
	// }

	size_t length(size_t seq_idx) {
		return start_idx[seq_idx + 1] - start_idx[seq_idx];
	}

	uint8_t get_base(size_t seq_idx, size_t offset) {
		uint64_t where = start_idx[seq_idx] + offset;
		return packed_seq[where / kCharsPerWord] >> (kCharsPerWord - 1 - where % kCharsPerWord) * 2 & 3;
	}

	void AppendSeq(const char *s, unsigned len) {
		for (unsigned i = 0; i < len; ++i) {
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

	void AppendSeq(word_t *s, unsigned len) {
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
			unsigned num_words = (len + kCharsPerWord - 1) / kCharsPerWord;
			unsigned bits_in_last_word = len * 2 % kBitsPerWord;

			// clear last word
			s[num_words-1] >>= kBitsPerWord - bits_in_last_word;

			// append to unused bits
			packed_seq.back() |= s[0] >> (kBitsPerWord - unused_bits_);

			for (unsigned i = 0; i + 1 < num_words; ++i) {
				s[i] = (s[i] << unused_bits_) | (s[i+1] >> (kBitsPerWord - unused_bits_));
			}
			s[num_words-1] <<= unused_bits_;

			unsigned num_words_to_append = (len - unused_bits_ / 2 + kCharsPerWord - 1) / kCharsPerWord;
			for (unsigned i = 0; i < num_words_to_append; ++i) {
				packed_seq.push_back(s[i]);
			}

			bits_in_last_word = end % kCharsPerWord * 2;
			unused_bits_ = kBitsPerWord - bits_in_last_word;
			if (unused_bits_ == kBitsPerWord) {
				packed_seq.push_back(word_t(0));
			}
		}
	}

	void get_seq(std::vector<word_t> &s, size_t seq_idx, int begin, int end = -1) {
		if (end == -1) {
			end = begin + length(seq_idx) - 1;
		}

		size_t first_word = (start_idx[seq_idx] + begin) / kCharsPerWord;
		size_t last_word = (start_idx[seq_idx+1] + end) / kCharsPerWord;
		int first_shift = start_idx[seq_idx] % kCharsPerWord * 2;

		s.clear();

		for (size_t i = first_word; i < last_word; ++i) {
			s.push_back((packed_seq[i] << first_shift) | (packed_seq[i+1] >> (kBitsPerWord - first_shift)));
		}

		int shift_clean = kBitsPerWord - (end - begin + 1) * 2 % kBitsPerWord;
		s.back() >>= shift_clean;
		s.back() <<= shift_clean;
	}
};

#endif