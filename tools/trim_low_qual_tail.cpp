#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "../kseq.h"

KSEQ_INIT(gzFile, gzread)

int main(int argc, char **argv) {
	if (argc < 4) {
		fprintf(stderr, "Usage: cat *.fq | %s <format=64|33> <qual_threshold> <min_len>\n", argv[0]);
		exit(1);
	}

	gzFile fp = gzdopen(fileno(stdin), "r");
    kseq_t *seq = kseq_init(fp); // kseq to read files
    int qb = atoi(argv[1]);
    int qt = atoi(argv[2]);
    int min_len = atoi(argv[3]);

    while (kseq_read(seq) >= 0) {
    	int i;
    	for (i = 0; i < seq->seq.l; ++i) {
    		if (seq->qual.s[i] - qb <= qt) {
    			break;
    		}
    	}
		if (i >= min_len) {
			seq->qual.s[i] = 0;
			seq->seq.s[i] = 0;
			printf("@%s", seq->name.s);
			if (seq->comment.s != NULL) {
				printf("%s", seq->comment.s);
			}
			printf("\n%s\n+\n%s\n", seq->seq.s, seq->qual.s);
		}
    }

    kseq_destroy(seq);
    gzclose(fp);
	return 0;
}