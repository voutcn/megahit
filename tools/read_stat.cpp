#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <algorithm>
#include "../kseq.h"

KSEQ_INIT(gzFile, gzread)

int main(int argc, char **argv) {
	if (argc > 1) {
		fprintf(stderr, "Usage: cat *.fq | %s\n", argv[0]);
		exit(1);
	}

	gzFile fp = gzdopen(fileno(stdin), "r");
    kseq_t *seq = kseq_init(fp); // kseq to read files
    int max_len = 0;
    int min_len = 999999999;
    int avg_len = 0;
    long long total_len = 0;
    long long num_reads = 0;

    while (kseq_read(seq) >= 0) {
    	++num_reads;
    	total_len += seq->seq.l;
    	max_len = std::max(seq->seq.l, (size_t)max_len);
    	min_len = std::min(seq->seq.l, (size_t)min_len);
    }

    printf("number: %lld\nlongest: %d\nshortest: %d\navg: %f\n", num_reads, max_len, min_len, total_len * 1.0 / num_reads);

    kseq_destroy(seq);
    gzclose(fp);
	return 0;
}