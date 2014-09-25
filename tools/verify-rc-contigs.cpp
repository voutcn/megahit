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

/**
 * @file verify-rc-contigs.cpp
 * @ brief To verify whether the reverse-complement of a contig is also output,
 *         when we do not mark the rc-nodes when traversing the graph
 * @author Dinghua Li (dhli@cs.hku.hk)
 * @date 2014-05-09
 */

#include <stdlib.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <unordered_set>

using namespace std;

inline char Complement(char c) {
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
        assert(false);
    }
    }
}

void ReverseComplement(string &s) {
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

int main(int argc, char **argv) {
    if (argc < 2) {
        cerr << "Usage :" << argv[0] << " [kmer_k]\n";
        cerr << "contigs should be input by stdin.\n";
        exit(1);
    }

    int kmer_k = atoi(argv[1]);
    long long id;
    string contig;
    unordered_set<string> st;
    long long num_palindrome = 0;
    long long num_verified = 0;
    long long num_not_verified = 0;
    long long num_loop = 0;
    long long num_not_verified_but_loop = 0;

    while (cin >> id >> contig) {
        if (contig.substr(0, kmer_k - 1) == contig.substr(contig.length() - kmer_k + 1)) {
            num_loop++;
        }
        st.insert(contig);
    }

    for (auto it = st.begin(); it != st.end(); ++it) {
        string s = *it;
        ReverseComplement(s);
        if (s == *it) {
            ++num_palindrome;
            ++num_verified;
        } else {
            if (st.count(s)) {
                ++num_verified;
            } else {
                ++num_not_verified;
                if (it->substr(0, kmer_k - 1) == it->substr(it->length() - kmer_k + 1)) {
                    num_not_verified_but_loop++;
                }
                cerr << *it << endl;
            }
        }
    }

    cout << "Number verified: " << num_verified << endl;
    cout << "Number palindrome: " << num_palindrome << endl;
    cout << "Number loop: " << num_loop << endl;
    cout << "Number not verified: " << num_not_verified << endl;
    cout << "Number not verified but loop: " << num_not_verified_but_loop << endl;

    return 0;
}