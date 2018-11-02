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
#include <cstdint>
#include <cstdlib>
#include <unordered_map>
#include <map>
#include <string>
#include <cstring>
using namespace std;

unordered_map<uint64_t, int> count;

int dna_char_map[256];
char rc_map[256];
uint64_t mask = 0x3FFFFFFFFFFFFFFFULL;

void init(int k) {
    dna_char_map['A'] = dna_char_map['a'] = 0;
    dna_char_map['C'] = dna_char_map['c'] = 1;
    dna_char_map['G'] = dna_char_map['g'] = 2;
    dna_char_map['T'] = dna_char_map['t'] = 3;
    dna_char_map['N'] = dna_char_map['n'] = 2;

    rc_map['A'] = rc_map['a'] = 'T';
    rc_map['C'] = rc_map['c'] = 'G';
    rc_map['G'] = rc_map['g'] = 'C';
    rc_map['T'] = rc_map['t'] = 'A';
    rc_map['N'] = rc_map['n'] = 'G';

    if (k < 32) {
        mask = (1ULL << k * 2) - 1;
    }
}

void rc(char *buf, int len) {
    int i, j;

    for (i = 0, j = len - 1; i < j; ++i) {
        swap(buf[i], buf[j]);
        buf[i] = rc_map[buf[i]];
        buf[j] = rc_map[buf[j]];
    }

    if (i == j) {
        buf[i] = rc_map[buf[i]];
    }
}

void gao(char *buf, int len, int k) {
    uint64_t key = 0;

    for (int i = 0; i < k; ++i) {
        key = key << 2;
        key += dna_char_map[buf[i]];
    }

    count[key]++;

    for (int i = k; i < len; ++i) {
        key &= mask;
        key <<= 2;
        key |= dna_char_map[buf[i]];
        count[key]++;
    }
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage %s [k_mer (<=32)]\n", argv[0]);
        fprintf(stderr, "Reads (fasta or fastq) should be inputed in stdin.\n");
        exit(1);
    }

    int k = atoi(argv[1]);
    long long numReads = 1;

    init(k);

    char buf[1024];
    // first read
    gets(buf);
    bool is_fq = buf[0] == '@';
    gets(buf);
    int len = strlen(buf);
    printf("Format: %s, Read Length: %d\n", is_fq ? "Fastq" : "Fasta", len);
    gao(buf, len, k);
    rc(buf, len);
    gao(buf, len, k);

    if (is_fq) {
        gets(buf);
        gets(buf);
    }

    while (gets(buf)) {
        gets(buf);

        if (strlen(buf) == len) {
            gao(buf, len, k);
            rc(buf, len);
            gao(buf, len, k);
            ++numReads;
        }

        if (is_fq) {
            gets(buf);
            gets(buf);
        }
    }

    printf("The number of reads: %lld\n", numReads);
    printf("The number of unique %d-mers: %lld\n", k, count.size());

    map<int, uint64_t> dis;

    for (auto it = count.begin(); it != count.end(); ++it) {
        dis[it->second]++;
    }

    for (auto it = dis.begin(); it != dis.end(); ++it) {
        fprintf(stderr, "%d %llu\n", it->first, it->second);
    }

    return 0;
}
