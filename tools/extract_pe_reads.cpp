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
    if (argc < 2) {
        fprintf(stderr, "%s <in.fa[fq]>\n", argv[0]);
        exit(1);
    }

    gzFile fp = gzopen(argv[1], "r");
    FILE *out_se = fopen((string(argv[1])+".se").c_str(), "w");
    FILE *out_pe = fopen((string(argv[1])+".pe").c_str(), "w");
    kseq_t *seq = kseq_init(fp); // kseq to read files

    unordered_map<string, Seq> reads; // header -> seq, qual, which strand

    while (kseq_read(seq) >= 0) {
        if (seq->name.l < 2 || seq->name.s[seq->name.l-2] != '/' || (seq->name.s[seq->name.l-1] != '1' && seq->name.s[seq->name.l-1] != '2')) {
            output_seq(seq, out_se);
        }

        string header = string(seq->name.s, seq->name.l - 2);
        auto it = reads.find(header);
        if (it != reads.end()) {
            int cmp = strcmp(it->second.name.c_str(), seq->name.s);
            if (cmp < 0) {
                output_seq(it->second, out_pe);
                output_seq(seq, out_pe);
            } else if (cmp > 0) {
                output_seq(seq, out_pe);
                output_seq(it->second, out_pe);
            } else {
                assert(false);
            }

            reads.erase(it);
        } else {
            // output_seq(seq, stdout);
            const char *comm = seq->comment.s ? seq->comment.s : "";
            const char *qual = seq->qual.s ? seq->qual.s : "";
            reads[header] = Seq(seq->name.s, comm, seq->seq.s, qual);
        }
    }

    for (auto &item : reads) {
        output_seq(item.second, out_se);
    }

    kseq_destroy(seq);
    gzclose(fp);
    fclose(out_se);
    fclose(out_pe);
    return 0;
}