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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>
#include <assert.h>

using namespace std;

#define READLEN 100
#define WORDS_PER_READ 7
#define MOD_WPR 127
#define DIV_WPR 7

#define BUFFER 8192

typedef long long int int64;

int main(int argc, char **argv) {

    char charmap[256];
    charmap['A'] = 0;
    charmap['C'] = 1;
    charmap['G'] = 2;
    charmap['T'] = 3;
    charmap['N'] = 2;

    char buffer[BUFFER+10];
    unsigned int bytes_read;
    unsigned int p = 0;
    char r[READLEN+20];
    unsigned int rp = 0;
    unsigned int out_buffer[7000+100];
    unsigned int op = 0;

#define FILL_BUFFER do { if (p>=bytes_read) { bytes_read = gzread(infile, buffer, BUFFER); p=0; } } while (0);

    unsigned int read_cnt = 0;

    for (int f=1; f<argc; ++f) {
        gzFile infile = gzopen(argv[f], "r");

        fprintf(stderr, "Processing file %s ...", argv[f]);
        fflush(stderr);

        FILL_BUFFER;
        while (1) {
            // read @.....\n
            if (buffer[p] != '>') goto finish;
            while (buffer[p]!='\n') {
                if (bytes_read==0) goto finish;
                p++;
                FILL_BUFFER;
            }
            p++;
            FILL_BUFFER;
            read_cnt++;
            // ACGT
            rp = 0;
            while (buffer[p]!='\n') {
                r[rp] = buffer[p];
                p++;
                rp++;
                FILL_BUFFER;
            }
            p++;
            FILL_BUFFER;
            if (buffer[p] == '+') {
                // read +
                while (buffer[p]!='\n') {
                    if (bytes_read==0) goto finish;
                    p++;
                    FILL_BUFFER;
                }
                p++;
                FILL_BUFFER;
                // read quality
                while (buffer[p]!='\n') {
                    if (bytes_read==0) goto finish;
                    p++;
                    FILL_BUFFER;
                }
                p++;
                FILL_BUFFER;
            }

            // binary-ize the read
            unsigned int bits = 0;
            for (int j = 0; j < READLEN; ++j) {
                if (!(j & 15)) {
                    if (j) out_buffer[op++] = bits;
                    bits = 0;
                }
                bits |= (charmap[r[j]] << (2 * (15-(j & 15))));
            }
            out_buffer[op++] = bits;
            if (op==7000) {
                fwrite(out_buffer, sizeof(unsigned int), op, stdout);
                op = 0;
            }

        }
finish:
        ;
        if (op) {
            fwrite(out_buffer, sizeof(unsigned int), op, stdout);
        }

        gzclose(infile);

        fprintf(stderr, "done. %u reads.\n", read_cnt);
        fflush(stderr);
    }

    fprintf(stderr, "read count = %u\n", read_cnt);

    return 0;

}

