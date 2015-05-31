#include "local_assembler.h"

#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

#include <zlib.h>
#include <pthread.h>
#include "lib_idba/sequence.h"
#include "lib_idba/hash_graph.h"
#include "lib_idba/contig_graph.h"

#include "utils.h"
#include "mem_file_checker-inl.h"
#include "kseq.h"
#include "histgram.h"
#include "bit_operation.h"
#include "atomic_bit_vector.h"

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

void LocalAssembler::BuildHashMapper(bool show_stat) {
    size_t sz = contigs_->size();
	size_t estimate_num_kmer = 0;

#pragma omp parallel for reduction(+: estimate_num_kmer)
	for (size_t i = 0; i < sz; ++i) {
		estimate_num_kmer += (contigs_->length(i) - seed_kmer_ + sparcity_) / sparcity_;
	}
	mapper_.reserve(estimate_num_kmer);

#pragma omp parallel for
    for (size_t i = 0; i < sz; ++i) {
    	AddToHashMapper_(mapper_, i, sparcity_);
    }

    if (show_stat) {
    	fprintf(stderr, "Mapper size :%lu\n", mapper_.size());
    }
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

int LocalAssembler::Match_(SequencePackage *read_lib, size_t read_id, int query_from, int query_to,
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
			return 0;
		}
	}

	return match_len;
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

		int match_bases = Match_(read_lib, read_id, query_from, query_to, contig_id, contig_from, contig_to, mapping_strand);
		if (match_bases > 0) {
			if (num_mapped > 0) {
				return false;
			} else {
				rec = tested_rec[tested-1];
				rec.mismatch = query_to - query_from + 1 - match_bases;
				num_mapped = 1;
			}
		}
	}

	return num_mapped == 1;
}

void LocalAssembler::EstimateInsertSize(bool show_stat) {
    for (unsigned lib_id = 0; lib_id < read_libs_.size(); ++lib_id) {
    	if (insert_sizes_[lib_id].first < -0.5) { continue; }
    	MappingRecord rec1, rec2;
    	Histgram<int> insert_hist;
    	size_t start_read_id = 0, end_read_id = 0;

    	while (insert_hist.size() < (1 << 20) && end_read_id < read_libs_[lib_id]->size()) {
    		start_read_id = end_read_id;
    		end_read_id = std::min(read_libs_[lib_id]->size(), start_read_id + size_t(2 << 20));
#pragma omp parallel for private(rec1, rec2)
	    	for (size_t i = start_read_id; i < end_read_id; i += 2) {
	    		if (MapToHashMapper_(mapper_, read_libs_[lib_id], i, rec1) &&
	    			MapToHashMapper_(mapper_, read_libs_[lib_id], i^1, rec2)) {
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
    	}
    	insert_hist.Trim(0.01);
    	insert_sizes_[lib_id] = tlen_t(insert_hist.mean(), insert_hist.sd());

    	if (show_stat) {
	    	fprintf(stderr, "Lib %d, mapped pairs: %u\n", lib_id, insert_hist.size());
	    	fprintf(stderr, "insert size: %.2lf sd: %.2lf\n", insert_hist.mean(), insert_hist.sd());
	    }
    }
}

int LocalAssembler::LocalRange_(int lib_id) {
	int local_range = read_libs_[lib_id]->max_read_len() - 1;
	if (insert_sizes_[lib_id].first >= read_libs_[lib_id]->max_read_len()) {
		local_range = std::min(2 * insert_sizes_[lib_id].first,
			                   insert_sizes_[lib_id].first + 3 * insert_sizes_[lib_id].second);
	}

	return local_range;
}

inline uint64_t PackMappingResult(uint64_t contig_offset, uint64_t is_mate, uint64_t mismatch,
								  uint64_t strand, uint64_t read_id, uint64_t lib_id, uint64_t num_libs) {
	assert(contig_offset < (1ULL << 14));
	return (contig_offset << 50) | (is_mate << 49) | (std::min(uint64_t(15), mismatch) << 45) |
		   (strand << 44) | (read_id * num_libs + lib_id);
}

bool LocalAssembler::AddToMappingDeque_(int lib_id, size_t read_id, const MappingRecord &rec, int local_range) {
	assert(read_id < read_libs_[lib_id]->size());
	assert(rec.contig_id < contigs_->size());

	int contig_len = contigs_->length(rec.contig_id);
	int read_len = read_libs_[lib_id]->length(read_id);
	if (rec.contig_to < local_range && rec.query_from != 0) {
		uint64_t res = PackMappingResult(rec.contig_to, 0, rec.mismatch, rec.strand, read_id, lib_id, read_libs_.size());
		while (!locks_.lock(rec.contig_id)) {
			continue;
		}
		mapped_f_[rec.contig_id].push_back(res);
		locks_.unset(rec.contig_id);
		return true;
	} else if (rec.contig_from + local_range >= contig_len && rec.query_to < read_len - 1) {
		uint64_t res = PackMappingResult(contig_len - 1 - rec.contig_from, 0, rec.mismatch, rec.strand, read_id, lib_id, read_libs_.size());
		while (!locks_.lock(rec.contig_id)) {
			continue;
		}
		mapped_r_[rec.contig_id].push_back(res);
		locks_.unset(rec.contig_id);
		return true;
	}

	return false;
}

bool LocalAssembler::AddMateToMappingDeque_(int lib_id, size_t read_id, const MappingRecord &rec1, const MappingRecord &rec2, bool mapped2, int local_range) {
	// assert(read_id < read_libs_[lib_id]->size());
	// assert((read_id ^ 1) < read_libs_[lib_id]->size());
	// assert(rec1.contig_id < contigs_->size());
	// assert(!mapped2 || rec2.contig_id < contigs_->size());

	if (mapped2 && rec2.contig_id == rec1.contig_id)
		return false;

	int contig_len = contigs_->length(rec1.contig_id);
	int read_len = read_libs_[lib_id]->length(read_id);

	if (rec1.contig_to < local_range && rec1.strand == 1) {
		uint64_t res = PackMappingResult(rec1.contig_to, 1, rec1.mismatch, rec1.strand, read_id^1, lib_id, read_libs_.size());
		while (!locks_.lock(rec1.contig_id)) {
			continue;
		}
		mapped_f_[rec1.contig_id].push_back(res);
		locks_.unset(rec1.contig_id);	
		return true;
	} else if (rec1.contig_from + local_range >= contig_len && rec1.strand == 0) {
		uint64_t res = PackMappingResult(contig_len - 1 - rec1.contig_from, 1, rec1.mismatch, rec1.strand, read_id^1, lib_id, read_libs_.size());
		while (!locks_.lock(rec1.contig_id)) {
			continue;
		}
		mapped_r_[rec1.contig_id].push_back(res);
		locks_.unset(rec1.contig_id);
		return true;
	}

	return false;
}

void LocalAssembler::MapToContigs() {
	mapped_f_.resize(contigs_->size());
	mapped_r_.resize(contigs_->size());
	locks_.reset(contigs_->size());

	max_read_len_ = 1;
	local_range_ = 0;

	for (unsigned lib_id = 0; lib_id < read_libs_.size(); ++lib_id) {
		int local_range = LocalRange_(lib_id);
		bool is_paired = insert_sizes_[lib_id].first >= read_libs_[lib_id]->max_read_len();

		local_range_ = std::max(local_range, local_range_);
		max_read_len_ = std::max(max_read_len_, (int)read_libs_[lib_id]->max_read_len());

    	size_t sz = read_libs_[lib_id]->size();
    	MappingRecord rec1, rec2;
    	size_t num_added = 0;

#pragma omp parallel for private(rec1, rec2) reduction(+: num_added)
    	for (size_t i = 0; i < sz; i += 2) {
    		bool map1 = MapToHashMapper_(mapper_, read_libs_[lib_id], i, rec1);
    		bool map2 = (i^1) < sz ? MapToHashMapper_(mapper_, read_libs_[lib_id], i^1, rec2) : false;

    		if (map1) {
    			num_added += AddToMappingDeque_(lib_id, i, rec1, local_range);
    			if (is_paired) {
    				num_added += AddMateToMappingDeque_(lib_id, i, rec1, rec2, map2, local_range);
    			}
    		}
    		if (map2) {
    			num_added += AddToMappingDeque_(lib_id, i^1, rec2, local_range);
    			if (is_paired) {
    				num_added += AddMateToMappingDeque_(lib_id, i^1, rec2, rec1, map1, local_range);
    			}
    		}
    	}
	    fprintf(stderr, "Lib %d: total %lu reads, added %lu reads for local assembly\n", lib_id, read_libs_[lib_id]->size(), num_added);
    }
}

void LocalAssembler::LocalAssemble() {
	omp_set_num_threads(1);

	std::vector<pthread_t> threads(num_threads_);
	std::vector<AssembleTask> tasks(num_threads_);
	for (int tid = 0; tid < num_threads_; ++tid) {
		tasks[tid].tid = tid;
		tasks[tid].local_assembler = this;
        pthread_create(&threads[tid], NULL, LocalAssembleThread_, (void *)&tasks[tid]);
	}

	for (int tid = 0; tid < num_threads_; ++tid) {
		pthread_join(threads[tid], NULL);
	}

	omp_set_num_threads(num_threads_);	
}

inline void LaunchIDBA(std::deque<Sequence> &reads, Sequence &contig_end,
					   std::deque<Sequence> &out_contigs,
					   int mink, int maxk, int step) {
	int local_range = contig_end.size();
	HashGraph hash_graph;
	hash_graph.reserve(4 * local_range);

	ContigGraph contig_graph;
	std::deque<ContigInfo> contig_infos;
	out_contigs.clear();

	int max_read_len = 0;
	for (unsigned i = 0; i < reads.size(); ++i) {
		max_read_len = std::max(max_read_len, (int)reads[i].size());
	}

    for (int kmer_size = mink; kmer_size <= std::min(maxk, max_read_len); kmer_size = std::min(kmer_size + step, max_read_len))
    {
        int64_t sum = 0;
        hash_graph.set_kmer_size(kmer_size);
        hash_graph.clear();
        for (int64_t i = 0; i < (int64_t)reads.size(); ++i)
        {
            if ((int)reads[i].size() < kmer_size)
                continue;

            Sequence seq(reads[i]);
            hash_graph.InsertKmers(seq);
            sum += seq.size() - kmer_size + 1;
        }

        Histgram<int> histgram = hash_graph.coverage_histgram();
        double mean = histgram.percentile(1 - 1.0 * local_range / hash_graph.num_vertices());
        double threshold = mean;

        hash_graph.InsertKmers(contig_end);
        for (int64_t i = 0; i < (int64_t)out_contigs.size(); ++i)
            hash_graph.InsertUncountKmers(out_contigs[i]);

        hash_graph.Assemble(out_contigs, contig_infos);
        contig_graph.set_kmer_size(kmer_size);
        contig_graph.Initialize(out_contigs, contig_infos);
        contig_graph.RemoveDeadEnd(kmer_size*2);
        
        contig_graph.RemoveBubble();
        contig_graph.IterateCoverage(kmer_size*2, 1, threshold);
        contig_graph.Assemble(out_contigs, contig_infos);

        if (out_contigs.size() == 1)
            break;
    }
}

void* LocalAssembler::LocalAssembleThread_(void *data) {
	omp_set_num_threads(1);

	AssembleTask *task = (AssembleTask*)data;
	LocalAssembler *la = task->local_assembler;
	int min_num_reads = la->local_range_ / la->max_read_len_;

	Sequence seq, contig_end;
	std::deque<Sequence> reads;
	std::deque<Sequence> out_contigs;

	for (size_t i = task->tid, csz = la->contigs_->size(); i < csz; i += la->num_threads_) {
		for (int strand = 0; strand < 2; ++strand) {
			std::deque<uint64_t> &mapped_reads = strand == 0 ? la->mapped_f_[i] : la->mapped_r_[i];
			if ((int)mapped_reads.size() <= min_num_reads) {
				continue;
			}

			// collect local reads, convert them into Sequence
			reads.clear();

			std::sort(mapped_reads.begin(), mapped_reads.end());
			int last_mapping_pos = -1;
			int pos_count = 0;

			for (size_t j = 0; j < mapped_reads.size(); ++j) {
				int pos = mapped_reads[j] >> 49;
				assert((pos >> 1) < la->contigs_->length(i));
				pos_count = pos == last_mapping_pos ? pos_count + 1 : 1;
				last_mapping_pos = pos;

				if (pos_count <= 3) {
					seq.clear();
					int lib_id = (mapped_reads[j] & ((1ULL << 44) - 1)) % la->read_libs_.size();
					uint64_t read_id = (mapped_reads[j] & ((1ULL << 44) - 1)) / la->read_libs_.size();
					assert(read_id < la->read_libs_[lib_id]->size());
					for (unsigned ri = 0, rsz = la->read_libs_[lib_id]->length(read_id); ri < rsz; ++ri) {
						seq.Append(la->read_libs_[lib_id]->get_base(read_id, ri));
					}
					reads.push_back(seq);
					// if (i == 0 && strand == 0) {
					// 	while (!la->locks_.lock(1)) {}
					// 	WriteFasta(std::cerr, seq, FormatString("read_%d", read_id));
					// 	la->locks_.unset(1);
					// }
				}
			}

			contig_end.clear();
			int cl = la->contigs_->length(i);
			if (strand == 0) {
				for (int j = 0, e = std::min(la->local_range_, cl); j < e; ++j) {
					contig_end.Append(la->contigs_->get_base(i, j));
				}
			} else {
				for (int j = std::max(0,  cl - la->local_range_); j < cl; ++j) {
					contig_end.Append(la->contigs_->get_base(i, j));
				}
			}

			out_contigs.clear();
			LaunchIDBA(reads, contig_end, out_contigs, la->local_kmin_, la->local_kmax_, la->local_step_);

			for (size_t j = 0; j < out_contigs.size(); ++j) {
				if (out_contigs[j].size() > la->min_contig_len_) {
					while (!la->locks_.lock(1)) {}
					WriteFasta(std::cout, out_contigs[j], FormatString("localcontig_%llu_strand_%d_id_%lu", i, strand, j));
					la->locks_.unset(1);
				}
			}
		}
	}

	return NULL;
}