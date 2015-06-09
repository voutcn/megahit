/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong
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

#ifndef ITERATE_EDGE_H_
#define ITERATE_EDGE_H_

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <zlib.h>
#include "definitions.h"

struct IterateGlobalData {
    char dna_map[256];

    std::string contig_file;
    std::string read_file;
    std::string read_format;
    std::string output_prefix;

    int kmer_k;
    int step;
    int next_k1; // = next_k + 1
    int num_cpu_threads;

    // stat
    int64_t num_of_reads;
    int64_t num_of_contigs;
    int64_t num_of_iterative_edges;
    int64_t num_of_remaining_reads;
};

#endif // ITERATE_EDGE_H_