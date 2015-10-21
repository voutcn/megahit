/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong & L3 Bioinformatics Limited
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* contact: Dinghua Li <dhli@cs.hku.hk> */

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "kseq.h"

#ifndef KSEQ_INITED
    #define KSEQ_INITED
    KSEQ_INIT(gzFile, gzread)
#endif

int main_trim_lowq_tail(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage: cat *.fq | %s <format=64|33> <qual_threshold> <min_len_remain>\n", argv[0]);
        exit(1);
    }

    gzFile fp = gzdopen(fileno(stdin), "r");
    kseq_t *seq = kseq_init(fp); // kseq to read files
    int qb = atoi(argv[1]);
    int qt = atoi(argv[2]);
    unsigned min_len = atoi(argv[3]);

    while (kseq_read(seq) >= 0) {
        unsigned i;

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