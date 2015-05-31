#ifndef LOCAL_ASSEMBLER_H__
#define LOCAL_ASSEMBLER_H__

#include <stdint.h>
#include <vector>
#include <string>
#include <deque>

#include "atomic_bit_vector.h"
#include "hash_table.h"
#include "sequence_package.h"
#include "kmer_plus.h"

struct LocalAssembler {
	typedef std::pair<double, double> tlen_t;
	typedef KmerPlus<2, uint32_t, uint64_t> kmer_plus_t;
	typedef Kmer<2, uint32_t> kmer_t;
	typedef HashTable<kmer_plus_t, kmer_t> mapper_t;

	enum ReadFormat {
		kFastx,
		kBinary,
	} read_format;

	struct MappingRecord {
		uint32_t contig_id;
		int32_t contig_from : 28;
		int32_t contig_to : 28;
		int16_t query_from : 15;
		int16_t query_to : 15;
		bool strand : 1;
		int mismatch: 9;
	};

	struct AssembleTask {
		int tid;
		LocalAssembler *local_assembler;
	};

	int min_contig_len_;	// only align reads to these contigs
	int seed_kmer_;			// kmer size for seeding
	int similarity_;		// similarity threshold for alignment
	int sparcity_;			// sparcity of hash mapper
	int min_mapped_len_;
	std::string local_filename_;
	int local_kmin_, local_kmax_, local_step_;
	int num_threads_;

	// auto calculated after calling EstimateInsertSize() and MapToContigs()
	int local_range_;
	int max_read_len_;

	SequencePackage *contigs_;
	mapper_t mapper_;
	std::vector<SequencePackage*> read_libs_;
	std::vector<tlen_t> insert_sizes_;
	AtomicBitVector locks_;

	std::vector<std::deque<uint64_t> > mapped_f_, mapped_r_;

	LocalAssembler(int min_contig_len, int seed_kmer, int sparcity)
		:min_contig_len_(min_contig_len), seed_kmer_(seed_kmer), sparcity_(sparcity),
		 contigs_(NULL) {
		similarity_ = 0.95;
		min_mapped_len_ = 75;
		num_threads_ = omp_get_max_threads();

		local_kmin_ = 21;
		local_kmax_ = 41;
		local_step_ = 6;
	}

	~LocalAssembler() {
		delete contigs_;
		for (auto it = read_libs_.begin(); it != read_libs_.end(); ++it) {
			delete *it;
		}
	}

	void set_kmer(int kmin, int kmax, int step) {
		local_kmin_ = kmin;
		local_kmax_ = kmax;
		local_step_ = step;
	}

	void set_mapping_threshold(double similarity, int mapping_len) {
		similarity_ = similarity;
		min_mapped_len_ = mapping_len;
	}

	void set_num_threads(int num_threads) {
		num_threads_ = num_threads;
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
	static void* LocalAssembleThread_(void *task);
};

#endif