/*
 *  MEGAHIT
 *  Copyright (C) 2014 The University of Hong Kong
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

#ifndef ITERATE_EDGE_H_
#define ITERATE_EDGE_H_

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <zlib.h>
#include "definitions.h"
#include "kmer_plus.h"
#include "hash_table.h"

struct IterateGlobalData {
    char dna_map[256];

    gzFile contigs_file;
    gzFile contigs_multi_file;
    gzFile addi_contig_file;
    gzFile addi_multi_file;
    gzFile read_file;
    FILE *output_edge_file;
    FILE *output_read_file;

    enum ReadFormats {
        kFasta,
        kFastq,
        kBinary,
    } read_format;

    int kmer_k;
    int step;
    int max_read_len;
    int num_cpu_threads;

    // stat
    int64_t num_of_reads;
    int64_t num_of_contigs;
    int64_t num_of_iterative_edges;
    int64_t num_of_remaining_reads;
};

#endif // ITERATE_EDGE_H_