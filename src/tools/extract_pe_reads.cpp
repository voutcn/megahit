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
#include <assert.h>
#include <algorithm>
#include <unordered_map>
#include <string>
#include "kseq.h"

using namespace std;

#ifndef KSEQ_INITED
    #define KSEQ_INITED
    KSEQ_INIT(gzFile, gzread)
#endif

struct Seq {
    string name;
    string comment;
    string seq;
    string qual;

    Seq(const string &name = "", const string &comment = "", const string &seq = "", const string &qual = ""):
        name(name), comment(comment), seq(seq), qual(qual) {}
};

void output_seq(kseq_t *seq, FILE *file) {
    const char *comm = seq->comment.s ? seq->comment.s : "";
    const char *qual = seq->qual.s ? seq->qual.s : "";
    if (seq->qual.l != 0) {
        fprintf(file, "@%s %s\n%s\n+%s\n%s\n", seq->name.s, comm, seq->seq.s, seq->name.s, qual);
    } else {
        fprintf(file, ">%s %s\n%s\n", seq->name.s, comm, seq->seq.s);
    }
}

void output_seq(const Seq &seq, FILE *file) {
    if (seq.qual.length() != 0) {
        fprintf(file, "@%s %s\n%s\n+%s\n%s\n", seq.name.c_str(), seq.comment.c_str(), seq.seq.c_str(), seq.name.c_str(), seq.qual.c_str());
    } else {
        fprintf(file, ">%s %s\n%s\n", seq.name.c_str(), seq.comment.c_str(), seq.seq.c_str());
    }
}

int main_extract_pe(int argc, char **argv) {
    if (argc < 2 || ((string(argv[1]) == "-" ) && argc < 3)) {
        fprintf(stderr, "%s <in.[fa|fq][.gz]> [out_prefix]\n", argv[0]);
        exit(1);
    }

    gzFile fp = string(argv[1]) == "-" ? gzdopen(fileno(stdin), "r") : gzopen(argv[1], "r");
    FILE *out_se = fopen((string(argc >= 3 ? argv[2] : argv[1])+".se").c_str(), "w");
    FILE *out_pe = fopen((string(argc >= 3 ? argv[2] : argv[1])+".pe").c_str(), "w");
    kseq_t *seq = kseq_init(fp); // kseq to read files

    unordered_map<string, Seq> reads; // header -> seq, qual, which strand

    while (kseq_read(seq) >= 0) {
        bool has12 = false;
        if (seq->name.l > 2 && seq->name.s[seq->name.l-2] == '/' && (seq->name.s[seq->name.l-1] == '1' || seq->name.s[seq->name.l-1] != '2')) {
            has12 = true;
        }

        string header = string(seq->name.s, seq->name.l - (has12 ? 2 : 0));
        auto it = reads.find(header);
        if (it != reads.end()) {
            int cmp = strcmp(it->second.name.c_str(), seq->name.s);
            if (cmp == 0) {
                cmp = strcmp(it->second.comment.c_str(), seq->comment.s);
            }
            if (cmp < 0) {
                output_seq(it->second, out_pe);
                output_seq(seq, out_pe);
            } else {
                output_seq(seq, out_pe);
                output_seq(it->second, out_pe);
            }

            reads.erase(it);
        } else {
            // output_seq(seq, stdout);
            const char *comm = seq->comment.s ? seq->comment.s : "";
            const char *qual = seq->qual.s ? seq->qual.s : "";
            reads[header] = Seq(seq->name.s, comm, seq->seq.s, qual);
        }
    }

    for (auto it = reads.begin(); it != reads.end(); ++it) {
        output_seq(it->second, out_se);
    }

    kseq_destroy(seq);
    gzclose(fp);
    fclose(out_se);
    fclose(out_pe);
    return 0;
}