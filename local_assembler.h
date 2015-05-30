#ifndef LOCAL_ASSEMBLER_H__
#define LOCAL_ASSEMBLER_H__

#include <stdint.h>
#include <vector>
#include <deque>

#include "atomic_bit_vector.h"
#include "hash_table.h"
#include "sequence_package.h"
#include "kmer_plus.h"

struct MappingRecord {
	uint32_t contig_id;
	int32_t contig_from : 28;
	int32_t contig_to : 28;
	int16_t query_from : 15;
	int16_t query_to : 15;
	bool strand : 1;
	int mismatch: 9;
};

struct LocalAssembler {
	typedef std::pair<double, double> tlen_t;
	typedef KmerPlus<2, uint32_t, uint64_t> kmer_plus_t;
	typedef Kmer<2, uint32_t> kmer_t;
	typedef HashTable<kmer_plus_t, kmer_t> mapper_t;

	enum ReadFormat {
		kFastx,
		kBinary,
	} read_format;

	int min_contig_len_;	// only align reads to these contigs
	int seed_kmer_;			// kmer size for seeding
	int similarity_;		// similarity threshold for alignment
	int sparcity_;			// sparcity of hash mapper
	int min_mapped_len_;

	SequencePackage *contigs_;
	mapper_t mapper_;
	std::vector<SequencePackage*> read_libs_;
	std::vector<tlen_t> insert_sizes_;
	AtomicBitVector locks_;

	std::vector<std::deque<uint64_t> > mapped_f, mapped_r;

	LocalAssembler(int min_contig_len, int seed_kmer)
		:min_contig_len_(min_contig_len), seed_kmer_(seed_kmer),
		 contigs_(NULL) {
		similarity_ = 0.95;
		sparcity_ = 8;
		min_mapped_len_ = 75;
	}

	~LocalAssembler() {
		delete contigs_;
		for (auto it = read_libs_.begin(); it != read_libs_.end(); ++it) {
			delete *it;
		}
	}

	void ReadContigs(const char *fastx_file_name);
	void BuildHashMapper(bool show_stat = true);
	void AddReadLib(const char *file_name, int file_type, bool is_paired);
	void EstimateInsertSize(bool show_stat = true);
	void MapToContigs();
	void LocalAssemble();

	void AddToHashMapper_(mapper_t &mapper, unsigned contig_id, int sparcity);
	int Match_(SequencePackage *read_lib, size_t read_id, int query_from, int query_to, size_t contig_id, int ref_from, int ref_to, bool strand);
	int LocalRange_(int lib_id);
	bool AddToMappingDeque_(int lib_id, size_t read_id, const MappingRecord &rec, int local_range);
	bool AddMateToMappingDeque_(int lib_id, size_t read_id, const MappingRecord &rec1, const MappingRecord &rec2, bool mapped2, int local_range); 
	bool MapToHashMapper_(const mapper_t &mapper, SequencePackage *read_lib, size_t read_id, MappingRecord &rec);
};

#endif