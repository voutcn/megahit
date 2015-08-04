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
#include <algorithm>
#include <vector>
#include <cstring>
#include <parallel/algorithm>

using namespace std;

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

inline void AddKmer(string &read, std::vector<string> &vs, int k_value) {
    if (read.length() >= k_value) {
        for (int i = 0; i + k_value <= read.length(); ++i) {
            string s = read.substr(i, k_value);
            string t = s;
            ReverseComplement(t);

            if (t < s) {
                vs.push_back(t);
            }
            else {
                vs.push_back(s);
            }
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " [k_value] [1.fa] [2.fa]" << endl;
        return 1;
    }

    string buf;
    int k_value = atoi(argv[1]);
    ifstream fa[2];
    ofstream uniq_out[2];
    vector<string> v_kmer[2];
    string read;

    for (int x = 0; x < 2; ++x) {
        fa[x].open(argv[2 + x]);
        uniq_out[x].open((string(argv[2 + x]) + ".uniq.kmer").c_str());

        while (getline(fa[x], buf)) {
            if (buf[0] == '>') {
                AddKmer(read, v_kmer[x], k_value);
                read.clear();
                continue;
            }

            read += buf;
        }

        if (read.length() >= k_value) {
            AddKmer(read, v_kmer[x], k_value);
            read.clear();
        }

        __gnu_parallel::sort(v_kmer[x].begin(), v_kmer[x].end());
        unsigned new_size = unique(v_kmer[x].begin(), v_kmer[x].end()) - v_kmer[x].begin();
        cout << argv[2 + x] << " contains " << new_size << " kmers, and " << v_kmer[x].size() - new_size << " redundant kmers." << endl;
        v_kmer[x].resize(new_size);
    }

    unsigned p1 = 0, p2 = 0;
    int num_uniq1 = 0, num_uniq2 = 0;
    int num_overlap = 0;

    while (p1 < v_kmer[0].size() && p2 < v_kmer[1].size()) {
        int cmp = strcmp(v_kmer[0][p1].c_str(), v_kmer[1][p2].c_str());

        if (cmp == 0) {
            ++p1;
            ++p2;
            ++num_overlap;
        }
        else if (cmp < 0) {
            uniq_out[0] << v_kmer[0][p1] << "\n";
            ++p1;
            ++num_uniq1;
        }
        else {
            uniq_out[1] << v_kmer[1][p2] << "\n";
            ++p2;
            ++num_uniq2;
        }
    }

    while (p1 < v_kmer[0].size()) {
        uniq_out[0] << v_kmer[0][p1] << "\n";
        ++p1;
        ++num_uniq1;
    }

    while (p2 < v_kmer[1].size()) {
        uniq_out[1] << v_kmer[1][p2] << "\n";
        ++p2;
        ++num_uniq2;
    }

    cout << "overlap: " << num_overlap << endl;
    cout << argv[2] << " unique kmer: " << num_uniq1 << endl;
    cout << argv[3] << " unique kmer: " << num_uniq2 << endl;

    uniq_out[0].close();
    uniq_out[1].close();

    return 0;
}