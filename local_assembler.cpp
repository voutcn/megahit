#include "local_assembler.h"

#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <zlib.h>
#include <vector>
#include <string>
#include <algorithm>

#include "utils.h"
#include "mem_file_checker-inl.h"
#include "kseq.h"
#include "histgram.h"
#include "bit_operation.h"

#ifndef KSEQ_INITED
#define KSEQ_INITED
KSEQ_INIT(gzFile, gzread)
#endif

void LocalAssembler::ReadContigs(const char *fastx_file_name) {
	gzFile fp = strcmp(fastx_file_name, "-") == 0 ? gzdopen(fileno(stdin), "r") : gzopen(fastx_file_name, "r");
	assert(fp != NULL);
	if (contigs_ == NULL) {
		contigs_ = new SequencePackage();
	}

    kseq_t *seq = kseq_init(fp); // kseq to read files
    int ctg_len = 0;
    while ((ctg_len = kseq_read(seq)) >= 0) {
    	if (ctg_len < (int)min_contig_len_) {
    		continue;
    	}
    	contigs_->AppendSeq(seq->seq.s, seq->seq.l);
	}

	kseq_destroy(seq);
	gzclose(fp);
}

void LocalAssembler::AddReadLib(const char *file_name, int file_type, bool is_paired) {
	SequencePackage *read_lib = new SequencePackage();
	read_libs_.push_back(read_lib);
	if (is_paired) {
		insert_sizes_.push_back(tlen_t(0, 0));
	} else {
		insert_sizes_.push_back(tlen_t(-1, -1));
	}

	gzFile fp = strcmp(file_name, "-") == 0 ? gzdopen(fileno(stdin), "r") : gzopen(file_name, "r");
	assert(fp != NULL);

	if (file_type == kBinary) {
		int seq_capacity = 16;
		uint32_t *seq = (uint32_t*) MallocAndCheck(sizeof(uint32_t) * seq_capacity, __FILE__, __LINE__);
		int len, num_words;

		while (gzread(fp, &len, sizeof(uint32_t)) != 0) {
			num_words = (len + SequencePackage::kCharsPerWord - 1) / SequencePackage::kCharsPerWord;
			if (num_words > seq_capacity) {
				seq_capacity = num_words;
				seq = (uint32_t*) ReAllocAndCheck(seq, sizeof(uint32_t) * seq_capacity, __FILE__, __LINE__);
			}
			assert(gzread(fp, &seq[0], sizeof(uint32_t) * num_words) == num_words * (int)sizeof(uint32_t));

			read_lib->AppendSeq(seq, len);
		}

		free(seq);
	} else {
		kseq_t *seq = kseq_init(fp); // kseq to read files
	    while (kseq_read(seq) >= 0) {
	    	read_lib->AppendSeq(seq->seq.s, seq->seq.l);
		}
		kseq_destroy(seq);
	}

	gzclose(fp);
}

inline uint64_t EncodeContigOffset(unsigned contig_id, unsigned contig_offset, bool strand) {
	return (uint64_t(contig_id) << 32) | (contig_offset << 1) | strand;
}

inline void DecodeContigOffset(uint64_t code, uint32_t &contig_id, uint32_t &contig_offset, bool &strand) {
	contig_id = code >> 32;
	contig_offset = (code & 0xFFFFFFFFULL) >> 1;
	strand = code & 1ULL;
}

void LocalAssembler::AddToHashMapper_(mapper_t &mapper, unsigned contig_id, int sparcity) {
	kmer_plus_t kp;
	kmer_t key;
	kp.ann = ~0ULL;	// special marker
	for (int i = 0, len = contigs_->length(contig_id); i + seed_kmer_ <= len; i += sparcity) {
		uint64_t full_offset = contigs_->start_idx[contig_id] + i;
		key.init(&contigs_->packed_seq[full_offset / SequencePackage::kCharsPerWord], full_offset % SequencePackage::kCharsPerWord, seed_kmer_);
		kp.kmer = key.unique_format(seed_kmer_);
		kmer_plus_t &kp_in_mapper = mapper.find_or_insert_with_lock(kp);

		if (kp_in_mapper.ann != ~0ULL) {
			kp_in_mapper.ann |= 1ULL << 63; // mark the highest bit as unused
		} else {
			kp_in_mapper.ann = EncodeContigOffset(contig_id, i, key != kp_in_mapper.kmer);
		}

		mapper.unlock(kp);
	}
}

inline uint32_t GetWord(uint32_t *first_word, uint32_t first_shift, int from, int len, bool strand) {
	int from_word_idx = (first_shift + from) / 16;
	int from_word_shift = (first_shift + from) % 16;
	uint32_t ret = *(first_word + from_word_idx) << from_word_shift * 2;
	assert(len <= 16);

	if (16 - from_word_shift < len) {
		ret |= *(first_word + from_word_idx + 1) >> (16 - from_word_shift) * 2;
	}	

	if (len < 16) {
		ret >>= (16 - len) * 2;
		ret <<= (16 - len) * 2;
	}

	if (strand == 1) {
		bit_operation::ReverseComplement(ret);
		ret <<= (16 - len) * 2;
	}

	return ret;
}

inline int Mismatch(uint32_t x, uint32_t y) {
	x ^= y;
	x |= x >> 1;
	x &= 0x55555555U;
	return __builtin_popcount(x);
}

bool LocalAssembler::Match_(SequencePackage *read_lib, size_t read_id, int query_from, int query_to,
							size_t contig_id, int ref_from, int ref_to, bool strand) {
	uint32_t *query_first_word = &read_lib->packed_seq[read_lib->start_idx[read_id] / 16];
	int query_shift = read_lib->start_idx[read_id] % 16;
	uint32_t *ref_first_word = &contigs_->packed_seq[contigs_->start_idx[contig_id] / 16]; 
	int ref_shift = contigs_->start_idx[contig_id] % 16;

	int match_len = query_to - query_from + 1;
	int threshold = similarity_ * match_len + 0.5;

	for (int i = query_from; i <= query_to; i += 16) {
		int len = std::min(16, query_to - i + 1);
		uint32_t qw = GetWord(query_first_word, query_shift, i, len, 0);
		int ref_i = strand == 0 ? ref_from + i - query_from : ref_to - (i + len - 1 - query_from);
		uint32_t rw = GetWord(ref_first_word, ref_shift, ref_i, len, strand);

		match_len -= Mismatch(qw, rw);
		if (match_len < threshold) {
			return false;
		}
	}

	return true;
}

bool LocalAssembler::MapToHashMapper_(const mapper_t &mapper, SequencePackage *read_lib, size_t read_id, MappingRecord &rec) {
	int len = read_lib->length(read_id);
	if (len < seed_kmer_) return false;

	int tested = 0;
	MappingRecord tested_rec[3];

	uint32_t *packed_seq = &read_lib->packed_seq[0];
	kmer_t kmer_f(packed_seq + read_lib->start_idx[read_id] / SequencePackage::kCharsPerWord,
				  read_lib->start_idx[read_id] % SequencePackage::kCharsPerWord, seed_kmer_);
	kmer_t kmer_r = kmer_f;
	kmer_r.ReverseComplement(seed_kmer_);
	int num_mapped = 0;
	uint32_t contig_id, contig_offset;
	bool contig_strand;

	for (int i = seed_kmer_ - 1; i < len; ++i) {
		if (i >= seed_kmer_) {
			uint8_t ch = read_lib->get_base(read_id, i);
			kmer_f.ShiftAppend(ch, seed_kmer_);
			kmer_r.ShiftPreappend(3 - ch, seed_kmer_);
		}

		bool query_strand = kmer_f.cmp(kmer_r, seed_kmer_) <= 0 ? 0 : 1;

		auto iter = mapper.find(query_strand == 0 ? kmer_f : kmer_r);
		if (iter == mapper.end() || (iter->ann >> 63) != 0) {
			continue;
		}

		DecodeContigOffset(iter->ann, contig_id, contig_offset, contig_strand);
		assert(contig_id < contigs_->size());
		assert(contig_offset < contigs_->length(contig_id));
		
		bool mapping_strand = contig_strand ^ query_strand;
		int contig_from = mapping_strand == 0 ? contig_offset - (i - seed_kmer_ + 1) : contig_offset - (len - 1 - i);
		int contig_to = mapping_strand == 0 ? contig_offset + seed_kmer_ - 1 + len - 1 - i : contig_offset + i;
		contig_from = std::max(contig_from, 0);
		contig_to = std::min((int)contigs_->length(contig_id) - 1, contig_to);

		if (contig_to - contig_from + 1 < min_mapped_len_) { continue; }

		int query_from = mapping_strand == 0 ? i - (seed_kmer_ - 1) - (contig_offset - contig_from) : i - (contig_to - contig_offset);
		int query_to = mapping_strand == 0 ? i - (seed_kmer_ - 1) + (contig_to - contig_offset) : i + (contig_offset - contig_from);

		bool has_tested = false;
		for (int j = 0; j < tested; ++j) {
			if (contig_id == tested_rec[j].contig_id &&
				contig_from == tested_rec[j].contig_from &&
				contig_to == tested_rec[j].contig_to &&
				query_from == tested_rec[j].query_from &&
				query_to == tested_rec[j].query_to &&
				mapping_strand == tested_rec[j].strand) {
				has_tested = true;
				break;
			}
		}

		if (has_tested) { continue; }
		else {
			if (tested >= 3) {
				tested--;
			}
			tested_rec[tested].contig_id = contig_id;
			tested_rec[tested].query_from = query_from;
			tested_rec[tested].query_to = query_to;
			tested_rec[tested].contig_from = contig_from;
			tested_rec[tested].contig_to = contig_to;
			tested_rec[tested].strand = mapping_strand;
			++tested;
		}

		if (Match_(read_lib, read_id, query_from, query_to, contig_id, contig_from, contig_to, mapping_strand)) {
			if (num_mapped > 0) {
				return false;
			} else {
				rec = tested_rec[tested-1];
				num_mapped = 1;
			}
		}
	}

	return num_mapped == 1;
}

void LocalAssembler::EstimateInsertSize(bool show_stat) {
	// first select the longest reads totaling ~1Mbp
	Histgram<int> contig_hist;
	for (size_t i = 0, sz = contigs_->size(); i < sz; ++i) {
		contig_hist.insert(contigs_->length(i));
	}

	int cutoff_len = 0; //contig_hist.Nx(kSizeToEstimateInsert);
	mapper_t mapper;

	if (show_stat) {
		fprintf(stderr, "Total contig size %lubp\n", contigs_->base_size());
    	fprintf(stderr, "Building hash mapper from contigs >= %dbp\n", cutoff_len);
    }

    size_t sz = contigs_->size();
#pragma omp parallel for
    for (size_t i = 0; i < sz; ++i) {
    	if (contigs_->length(i) >= (unsigned)cutoff_len) {
    		AddToHashMapper_(mapper, i, sparcity_);
    	}
    }

    if (show_stat) {
    	fprintf(stderr, "Mapper size :%lu\n", mapper.size());
    }

    for (unsigned lib_id = 0; lib_id < read_libs_.size(); ++lib_id) {
    	size_t sz = read_libs_[lib_id]->size();
    	MappingRecord rec1, rec2;
    	Histgram<int> insert_hist;

#pragma omp parallel for private(rec1, rec2)
    	for (size_t i = 0; i < std::min(sz, size_t(1 << 20)); i += 2) {
    		if (MapToHashMapper_(mapper, read_libs_[lib_id], i, rec1) &&
    			MapToHashMapper_(mapper, read_libs_[lib_id], i^1, rec2)) {
    			if (rec1.contig_id == rec2.contig_id && rec1.strand != rec2.strand) {
    				int insert_size = -1;
    				if (rec1.strand == 0) {
    					insert_size = rec2.contig_to + read_libs_[lib_id]->length(i^1) - rec2.query_to - (rec1.contig_from - rec1.query_from);
    				} else {
    					insert_size = rec1.contig_to + read_libs_[lib_id]->length(i) - rec1.query_to - (rec2.contig_from - rec2.query_from);
    				}

    				if (insert_size >= (int)read_libs_[lib_id]->length(i) &&
    					insert_size >= (int)read_libs_[lib_id]->length(i^1)) {
    					insert_hist.insert(insert_size);
    				}
    			}
    		}
    	}
    	insert_hist.Trim(0.01);

    	if (show_stat) {
	    	fprintf(stderr, "Lib %d, mapped pairs: %u\n", lib_id, insert_hist.size());
	    	fprintf(stderr, "insert size: %.2lf sd: %.2lf\n", insert_hist.mean(), insert_hist.sd());
	    }
    }
}