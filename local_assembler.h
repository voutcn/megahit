#ifndef LOCAL_ASSEMBLER_H__
#define LOCAL_ASSEMBLER_H__

#include <stdint.h>
#include <vector>
#include <deque>
#include "sequence_package.h"

struct LocalAssembler {
	typedef std::pair<double, double> tlen_t;

	enum ReadFormat {
		kFastx,
		kBinary,
	} read_format;

	unsigned min_contig_len_;	// only align reads to these contigs
	unsigned seed_kmer_;			// kmer size for seeding

	SequencePackage *contigs_;
	std::vector<SequencePackage*> read_libs_;
	std::vector<tlen_t> insert_sizes_;

	std::vector<std::deque<uint64_t> > mapped_f, mapped_r;

	LocalAssembler(unsigned min_contig_len, unsigned seed_kmer)
		:min_contig_len_(min_contig_len), seed_kmer_(seed_kmer),
		 contigs_(NULL) {}

	~LocalAssembler() {
		delete contigs_;
		for (auto it = read_libs_.begin(); it != read_libs_.end(); ++it) {
			delete *it;
		}
	}

	void ReadContigs(const char *fastx_file_name);
	void AddReadLib(const char *file_name, int file_type, bool is_paired);
	void EstimateInsertSize();
	void MapToContigs();
	void LocalAssemble();
};

#endif