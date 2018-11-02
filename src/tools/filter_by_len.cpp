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