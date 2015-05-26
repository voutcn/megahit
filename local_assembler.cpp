#include "local_assembler.h"

#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <zlib.h>
#include <vector>
#include <string>

#include "utils.h"
#include "mem_file_checker-inl.h"
#include "kseq.h"

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
		uint32_t seq_capacity = 16;
		uint32_t *seq = (uint32_t*) MallocAndCheck(sizeof(uint32_t) * seq_capacity, __FILE__, __LINE__);
		uint32_t len, num_words;

		while (gzread(fp, &len, sizeof(uint32_t)) != 0) {
			num_words = (len + SequencePackage::kCharsPerWord - 1) / SequencePackage::kCharsPerWord;
			if (num_words > seq_capacity) {
				seq_capacity = num_words;
				seq = (uint32_t*) ReAllocAndCheck(seq, sizeof(uint32_t) * seq_capacity, __FILE__, __LINE__);
			}
			assert(gzread(fp, &seq[0], sizeof(uint32_t) * num_words) == num_words * sizeof(uint32_t));

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