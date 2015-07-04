#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <algorithm>
#include "histgram.h"
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

    Histgram<long long> hist;

    while (kseq_read(seq) >= 0) {
        if (seq->seq.l >= min_len) {
            hist.insert(seq->seq.l);
            printf(">%s %s\n%s\n", seq->name.s, seq->comment.s, seq->seq.s);
        }
    }

    long long total_bases = hist.sum();

    fprintf(stderr, "%d contigs, total %lld bp, min %lld bp, max %lld bp, avg %d bp, N50 %lld bp\n",
                    (int)hist.size(), total_bases, hist.minimum(), hist.maximum(), int(hist.mean() + 0.5), hist.Nx(total_bases * 0.5));

    kseq_destroy(seq);
    gzclose(fp);
	return 0;
}