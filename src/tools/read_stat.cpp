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
#include "kseq.h"


#ifndef KSEQ_INITED
    #define KSEQ_INITED
    KSEQ_INIT(gzFile, gzread)
#endif

int main_read_stat(int argc, char **argv) {
    if (argc > 1) {
        fprintf(stderr, "Usage: cat *.fq | %s\n", argv[0]);
        exit(1);
    }

    gzFile fp = gzdopen(fileno(stdin), "r");
    kseq_t *seq = kseq_init(fp); // kseq to read files
    int max_len = 0;
    int min_len = 999999999;
    long long total_len = 0;
    long long num_reads = 0;

    while (kseq_read(seq) >= 0) {
        ++num_reads;
        total_len += seq->seq.l;
        max_len = std::max(seq->seq.l, (size_t)max_len);
        min_len = std::min(seq->seq.l, (size_t)min_len);
    }

    double avg_len = total_len * 1.0 / num_reads;

    printf("number reads: %lld\ntotal size: %lld\nlongest: %d\nshortest: %d\navg: %lf\n", num_reads, total_len, max_len, min_len, avg_len);

    kseq_destroy(seq);
    gzclose(fp);
    return 0;
}