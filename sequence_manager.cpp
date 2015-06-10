#include "sequence_manager.h"

#include <assert.h>
#include <zlib.h>

#include "utils.h"
#include "edge_reader.h"
#include "bit_operation.h"

void SequenceManager::set_file(const std::string &file_name) {
	assert(f_type != kMegahitEdges && f_type != kSortedEdges);
	assert(files_.size() == 0);
	assert(kseq_readers_.size() == 0);

	files_.resize(1);
	kseq_readers_.resize(1);

	files_[0] = file_name == "-" ? gzdopen(fileno(stdin), "r") : gzopen(file_name.c_str(), "r");
	assert(files_[0] != NULL);
	if (f_type == kFastxReads || f_type == kMegahitContigs) {
		kseq_readers_[0] = kseq_init(files_[0]);
		assert(kseq_readers_[0] != NULL);
	}
}

void SequenceManager::set_pe_files(const std::string &file_name1, const std::string &file_name2) {
	assert(f_type == kFastxReads && r_type == kPaired);
	assert(files_.size() == 0);
	assert(kseq_readers_.size() == 0);

	files_.resize(2);
	kseq_readers_.resize(2);

	files_[0] = file_name1 == "-" ? gzdopen(fileno(stdin), "r") : gzopen(file_name1.c_str(), "r");
	files_[1] = file_name2 == "-" ? gzdopen(fileno(stdin), "r") : gzopen(file_name2.c_str(), "r");

	assert(files_[0] != NULL);
	if (f_type == kFastxReads || f_type == kMegahitContigs) {
		for (int i = 0; i < 2; ++i) {
			kseq_readers_[i] = kseq_init(files_[i]);
			assert(kseq_readers_[i] != NULL);
		}
	}
}

void SequenceManager::set_edge_files(const std::string &file_prefix, int num_files) {
	assert(f_type == kSortedEdges || f_type == kMegahitEdges);
	assert(files_.size() == 0);
	assert(kseq_readers_.size() == 0);

	assert(!edge_reader_inited_);
	edge_reader_.init(file_prefix, num_files);
	edge_reader_inited_ = true;
}

int64_t SequenceManager::ReadShortReads(int64_t max_num, int64_t max_num_bases, bool append, bool reverse) {
	if (!append) {
		package_->clear();
	}

	max_num = (max_num + 1) / 2 * 2;
	int64_t num_bases = 0;
	if (f_type == kFastxReads) {
		if (r_type == kPaired) {
			for (int64_t i = 0; i < max_num; i += 2) {
				if (kseq_read(kseq_readers_[0]) >= 0) {
					assert(kseq_read(kseq_readers_[1]) >= 0);

					if (reverse) {
						package_->AppendReverseSeq(kseq_readers_[0]->seq.s, kseq_readers_[0]->seq.l);
						package_->AppendReverseSeq(kseq_readers_[1]->seq.s, kseq_readers_[1]->seq.l);
					} else {
						package_->AppendSeq(kseq_readers_[0]->seq.s, kseq_readers_[0]->seq.l);
						package_->AppendSeq(kseq_readers_[1]->seq.s, kseq_readers_[1]->seq.l);
					}

					num_bases += kseq_readers_[0]->seq.l + kseq_readers_[1]->seq.l;
					if (num_bases >= max_num_bases) {
						return i + 2;
					}
				} else {
					assert(kseq_read(kseq_readers_[1]) < 0);
					return i;
				}
			}
		} else {
			for (int64_t i = 0; i < max_num; ++i) {
				if (kseq_read(kseq_readers_[0]) >= 0) {
					if (reverse) {
						package_->AppendReverseSeq(kseq_readers_[0]->seq.s, kseq_readers_[0]->seq.l);
					} else {
						package_->AppendSeq(kseq_readers_[0]->seq.s, kseq_readers_[0]->seq.l);
					}

					num_bases += kseq_readers_[0]->seq.l;
					if (num_bases >= max_num_bases && i % 2 == 1) {
						return i + 1;
					}
				} else {
					return i;
				}
			}
		}
		return max_num;
	} 
	else if (f_type == kBinaryReads) {
		assert(!reverse);
		uint32_t read_len;
		for (int64_t i = 0; i < max_num; ++i) {
			if (gzread(files_[0], &read_len, sizeof(read_len)) == 0) {
				return i;
			}

			int num_words = DivCeiling(read_len, 16);
			if (buf_.size() < (unsigned)num_words) { buf_.resize(num_words); }
			assert((unsigned)gzread(files_[0], &buf_[0], sizeof(uint32_t) * num_words) == num_words * sizeof(uint32_t));
			package_->AppendSeq(&buf_[0], read_len);
			num_bases += read_len;
			if (read_len >= max_num_bases && i % 2 == 1) {
				return i + 1;
			}
		}
		return max_num;
	}

	assert(false);
}

int64_t SequenceManager::ReadEdges(int64_t max_num, bool append) {
	if (!append) {
		multi_->clear();
		package_->clear();
	}

	if (f_type == kMegahitEdges) {
		for (int64_t i = 0; i < max_num; ++i) {
			uint32_t *next_edge = edge_reader_.NextUnsortedEdge();
			if (next_edge == NULL) { return i; }
			package_->AppendSeq(next_edge, edge_reader_.kmer_k + 1);
			multi_->push_back(next_edge[edge_reader_.words_per_edge - 1] & kMaxMulti_t);
		}
		return max_num;
	} else if (f_type == kSortedEdges) {
		for (int64_t i = 0; i < max_num; ++i) {
			uint32_t *next_edge = edge_reader_.NextSortedEdge();
			if (next_edge == NULL) { return i; }
			package_->AppendSeq(next_edge, edge_reader_.kmer_k + 1);
			multi_->push_back(next_edge[edge_reader_.words_per_edge - 1] & kMaxMulti_t);
		}
		return max_num;
	}

	assert(false);
}


int64_t SequenceManager::ReadEdgesWithFixedLen(int64_t max_num, bool append) {
	if (!append) {
		multi_->clear();
		package_->clear();
	}

	if (f_type == kMegahitEdges) {
		for (int64_t i = 0; i < max_num; ++i) {
			uint32_t *next_edge = edge_reader_.NextUnsortedEdge();
			if (next_edge == NULL) { return i; }
			package_->AppendFixedLenSeq(next_edge, edge_reader_.kmer_k + 1);
			multi_->push_back(next_edge[edge_reader_.words_per_edge - 1] & kMaxMulti_t);
		}
		return max_num;
	} else if (f_type == kSortedEdges) {
		for (int64_t i = 0; i < max_num; ++i) {
			uint32_t *next_edge = edge_reader_.NextSortedEdge();
			if (next_edge == NULL) { return i; }
			package_->AppendFixedLenSeq(next_edge, edge_reader_.kmer_k + 1);
			multi_->push_back(next_edge[edge_reader_.words_per_edge - 1] & kMaxMulti_t);
		}
		return max_num;
	}

	assert(false);
}

int64_t SequenceManager::ReadMegahitContigs(int64_t max_num, int64_t max_num_bases, bool append, bool reverse, 
										    int discard_flag, bool extend_loop, bool calc_depth) {
	assert(f_type == kMegahitContigs);
	assert(!(calc_depth && multi_ == NULL));
	assert(!((discard_flag & (contig_flag::kLoop | contig_flag::kIsolated)) && extend_loop)); // loop must be isolated

	if (!append) {
		if (multi_ != NULL) { multi_->clear(); }
		package_->clear();
	}

	int64_t num_bases = 0;

	for (int64_t i = 0; i < max_num; ++i) {
		if (kseq_read(kseq_readers_[0]) >= 0) {
			if ((int)kseq_readers_[0]->seq.l < min_len_ ) { --i; continue; }

			// comment = "flag=x multi=xx.xxxx"
			if (discard_flag & (kseq_readers_[0]->comment.s[5] - '0')) {
				--i;
				continue;
			}

			if (extend_loop && kseq_readers_[0]->seq.l >= k_to_ + 1U &&
				((kseq_readers_[0]->comment.s[5] - '0') & contig_flag::kLoop)) {

				std::string ss(kseq_readers_[0]->seq.s);
				for (int i = 0; i < k_to_ + 1 - k_from_; ++i) {
					ss.push_back(ss[i + k_from_]);
				}

				if (reverse) {
					package_->AppendReverseSeq(ss.c_str(), kseq_readers_[0]->seq.l);
				} else {
					package_->AppendSeq(ss.c_str(), kseq_readers_[0]->seq.l);
				}
			} else {
				if (reverse) {
					package_->AppendReverseSeq(kseq_readers_[0]->seq.s, kseq_readers_[0]->seq.l);
				} else {
					package_->AppendSeq(kseq_readers_[0]->seq.s, kseq_readers_[0]->seq.l);
				}
			}

			if (multi_ != NULL) {
				if (calc_depth) {
					double depth_from = atof(kseq_readers_[0]->comment.s + 13);

					int num_kmer = kseq_readers_[0]->seq.l - k_from_ + 1;
	                int num_nextk1 = kseq_readers_[0]->seq.l - (k_to_ + 1) + 1;
	                int internal_max = std::min(k_to_ + 1 - k_from_ + 1, num_nextk1);
	                int num_external = internal_max - 1;
	                int num_internal = num_kmer - num_external * 2;

	                double exp_num_kmer = (double)num_external * (num_external + 1) / (k_to_ + 1 - k_from_ + 1)
	                                      + (double)internal_max / (k_to_ + 1 - k_from_ + 1) * num_internal;
	                exp_num_kmer *= depth_from;
	                multi_->push_back(std::min(int(exp_num_kmer * k_from_ / (k_to_ + 1) / num_nextk1 + 0.5), kMaxMulti_t));
				} else {
					multi_->push_back(1);
				}
			}

			num_bases += kseq_readers_[0]->seq.l;
			if (num_bases >= max_num_bases) {
				return i + 1;
			}
		} else {
			return i;
		}
	}

	return max_num;
}

void SequenceManager::WriteBinarySequences(FILE *file, bool reverse, int64_t from, int64_t to) {
	if (to == -1) {
		to = package_->size() - 1;
	}
	uint32_t len;
	std::vector<uint32_t> s;

	for (int64_t i = from; i <= to; ++i) {
		len = package_->length(i);
		package_->get_seq(s, i);

		if (reverse) {
			for (int j = 0; j < (int)s.size(); ++j) {
				s[j] = bit_operation::Reverse(s[j]);
			}
			for (int j = 0, k = s.size() - 1; j < k; ++j, --k) {
				std::swap(s[j], s[k]);
			}

			int shift = (16 - len % 16)  * 2;
			if (shift != 32) {
				for (int j = 0; j < (int)s.size() - 1; ++j) {
					s[j] = (s[j] << shift) | (s[j+1] >> (32 - shift));
				}
				s.back() <<= shift;
			}
		}

		fwrite(&len, sizeof(uint32_t), 1, file);
		fwrite(&s[0], sizeof(uint32_t), s.size(), file);
	}
}

