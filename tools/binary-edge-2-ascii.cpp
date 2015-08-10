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

#include <cstdio>
#include <string>
#include <cstdlib>
#include <cassert>

#include "kseq.h"
#include "utils.h"

using namespace std;

int main_binary_edge_to_ascii(int argc, char const *argv[]) {
    if (argc < 5) {
        fprintf(stderr, "Usage: %s <edge_file_prefix> <num_files> <has_header=0|1> <kmer_k>\n", argv[0]);
        exit(1);
    }

    const char *file_prefix = argv[1];
    int num_files = atoi(argv[2]);
    int has_header = atoi(argv[3]);
    int kmer_k = atoi(argv[4]);
    int words_per_edge = ((kmer_k + 1) * 2 + 16 + 31) / 32;
    unsigned int buf[1024];

    for (int i = 0; i < num_files; ++i) {
        FILE *cur_file = fopen(FormatString("%s.%d", file_prefix, i), "rb");
        assert(cur_file != NULL);

        if (has_header && i == 0) {
            int k, w;
            fread(&k, sizeof(unsigned int), 1, cur_file);
            assert(k == kmer_k);
            fread(&w, sizeof(unsigned int), 1, cur_file);
            assert(w == words_per_edge);
        }

        while (fread(buf, sizeof(unsigned int), words_per_edge, cur_file) == words_per_edge) {
            for (int k = 0; k < kmer_k + 1; ++k) {
                printf("%c", "ACGT"[3 - ((buf[k / 16] >> (15 - k % 16) * 2) & 3)]);
            }

            puts("");

            for (int k = kmer_k; k >= 0; --k) {
                printf("%c", "ACGT"[((buf[k / 16] >> (15 - k % 16) * 2) & 3)]);
            }

            puts("");
        }
    }

    return 0;
}