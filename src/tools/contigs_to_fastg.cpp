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

#include <string>
#include <map>
#include <vector>
#include <cstdio>
#include <cassert>
#include <zlib.h>

#include "kseq.h"

using namespace std;

#ifndef KSEQ_INITED
    #define KSEQ_INITED
    KSEQ_INIT(gzFile, gzread)
#endif

char Comp(char c) {
    switch (c) {
    case 'A':
    case 'a':
        return 'T';

    case 'C':
    case 'c':
        return 'G';

    case 'G':
    case 'g':
        return 'C';

    case 'T':
    case 't':
        return 'A';

    default:
        assert(false);
    }
}

string RevComp(const string &s) {
    string ret;

    for (unsigned i = 0; i < s.length(); ++i) {
        ret.push_back(Comp(s[s.length() - 1 - i]));
    }

    return ret;
}

string NodeName(int i, int len, double mul, bool is_rc) {
    static char buf[10240];
    sprintf(buf, "NODE_%d_length_%d_cov_%.4f_ID_%d", i, len, mul, i * 2 - 1);

    if (is_rc) return string(buf) + "'";
    else return string(buf);
}

int main_contig2fastg(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <kmer_size> <k_{kmer_size}.contigs.fa>\n", argv[0]);
        exit(1);
    }

    unsigned k = atoi(argv[1]);
    gzFile fp = gzopen(argv[2], "r");
    assert(fp != NULL);
    kseq_t *seq = kseq_init(fp); // kseq to read files

    vector<string> ctgs;
    vector<double> muls;
    vector<string> node_names;
    vector<string> rev_node_names;
    map<string, vector<int> > start_kmer_to_id;

    while (kseq_read(seq) >= 0) {
        if (seq->seq.l < k + 1) {
            continue;
        }

        double mul;
        assert(sscanf(seq->comment.s + 7, "multi=%lf", &mul) == 1);

        muls.push_back(mul);
        ctgs.push_back(string(seq->seq.s));
    }

    for (int i = 0; i < (int)ctgs.size(); ++i) {
        start_kmer_to_id[ctgs[i].substr(0, k)].push_back(i + 1);
        start_kmer_to_id[RevComp(ctgs[i].substr(ctgs[i].length() - k))].push_back(-i - 1);
    }

    for (int i = 0; i < (int)ctgs.size(); ++i) {
        node_names.push_back(NodeName(i + 1, ctgs[i].length(), muls[i], false));
        rev_node_names.push_back(NodeName(i + 1, ctgs[i].length(), muls[i], true));
    }

    for (int i = 0; i < (int)ctgs.size(); ++i) {
        for (int dir = 0; dir < 2; ++dir) {
            string header = dir == 0 ? node_names[i] : rev_node_names[i];
            header = ">" + header;
            string s = dir == 0 ? ctgs[i] : RevComp(ctgs[i]);
            auto mit = start_kmer_to_id.find(s.substr(s.length() - k));

            if (mit != start_kmer_to_id.end()) {
                for (unsigned j = 0; j < mit->second.size(); ++j) {
                    if (j == 0) {
                        header += ":";
                    }
                    else {
                        header += ",";
                    }

                    if (mit->second[j] > 0) {
                        header += node_names[mit->second[j] - 1];
                    }
                    else {
                        header += rev_node_names[-mit->second[j] - 1];
                    }
                }
            }

            header += ";";
            printf("%s\n%s\n", header.c_str(), s.c_str());
        }
    }

    return 0;
}