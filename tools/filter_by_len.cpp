#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <algorithm>
#include "kseq.h"


#ifndef KSEQ_INITED
#define KSEQ_INITED
KSEQ_INIT(gzFile, gzread)
#endif

int main_filter_by_len(int argc, char **argv) {
	if (argc < 2) {
		fprintf(stderr, "Usage: cat contigs.fa | %s <min_len>\n", argv[0]);
		exit(1);
	}

	gzFile fp = gzdopen(fileno(stdin), "r");
    kseq_t *seq = kseq_init(fp); // kseq to read files
    unsigned min_len = atoi(argv[1]);

    while (kseq_read(seq) >= 0) {
        if (seq->seq.l >= min_len) {
            printf(">%s %s\n%s\n", seq->name.s, seq->comment.s, seq->seq.s);
        }
    }

    kseq_destroy(seq);
    gzclose(fp);
	return 0;
}