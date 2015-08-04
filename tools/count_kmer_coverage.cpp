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

#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include "../container/hash_map.h"

using namespace std;

unsigned int BKDRHash(const string &str) {
    unsigned int seed = 131; // 31 131 1313 13131 131313 etc..
    unsigned int hash = 0;

    for (unsigned i = 0; i < str.length(); ++i) {
        hash = hash * seed + str[i];
    }

    return (hash & 0x7FFFFFFF);
}

struct StringHash: public std::unary_function<string, uint64_t> {
    uint64_t operator ()(const string &value) const {
        return BKDRHash(value);
    }
};

static inline char Complement(char c) {
    switch (c) {
    case 'A': {
        return 'T';
    }

    case 'C': {
        return 'G';
    }

    case 'G': {
        return 'C';
    }

    case 'T': {
        return 'A';
    }

    default: {
        return 'N';
    }
    }
}

static inline void ReverseComplement(string &s) {
    int i, j;

    for (i = 0, j = s.length() - 1; i < j; ++i, --j) {
        std::swap(s[i], s[j]);
        s[i] = Complement(s[i]);
        s[j] = Complement(s[j]);
    }

    if (i == j) {
        s[i] = Complement(s[i]);
    }
}

inline void AddKmer(string &read, HashMap<string, int, StringHash > &cnt, int k_value) {
    if (read.length() >= k_value) {
        for (int i = 0; i + k_value <= read.length(); ++i) {
            string s = read.substr(i, k_value);
            string t = s;
            ReverseComplement(t);
            auto it = cnt.find(s);

            if (it != cnt.end()) ++it->second;

            it = cnt.find(t);

            if (it != cnt.end()) ++it->second;
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " [k_value] [kmer_file] [reads.fa]" << endl;
        return 1;
    }

    ifstream kmer_file(argv[2]);
    ifstream read_file(argv[3]);

    string buf;
    string read;
    int k_value = atoi(argv[1]);
    HashMap<string, int, StringHash > cnt;

    while (getline(kmer_file, buf)) {
        cnt[buf] = 0;
    }

    cout << "Read kmer file done!" << endl;

    vector<string> vs;

    while (getline(read_file, buf)) {
        if (buf[0] == '>') {
            if (read.length() >= k_value) {
                vs.push_back(read);

                if (vs.size() >= 1024 * 1024) {
                    #pragma omp parallel for

                    for (unsigned i = 0; i < vs.size(); ++i) {
                        AddKmer(vs[i], cnt, k_value);
                    }

                    vs.clear();
                }
            }

            read.clear();
            continue;
        }

        read += buf;
    }

    if (read.length() >= k_value) {
        if (read.length() >= k_value) {
            vs.push_back(read);
        }
    }

    #pragma omp parallel for

    for (unsigned i = 0; i < vs.size(); ++i) {
        AddKmer(vs[i], cnt, k_value);
    }

    vs.clear();

    HashMap<int, int> histo;

    for (auto it = cnt.begin(); it != cnt.end(); ++it) {
        ++histo[it->second];
        cerr << it->first << ' ' << it->second << "\n";
    }

    for (auto it = histo.begin(); it != histo.end(); ++it) {
        cout << it->first << ' ' << it->second << endl;
    }

    return 0;
}