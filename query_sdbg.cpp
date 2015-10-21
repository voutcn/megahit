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

#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <omp.h>
#include "succinct_dbg.h"
#include "assembly_algorithms.h"

using namespace std;

void CheckAndPrintUsage(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s [DBG_name] [mode: 1 for query many kmers, 2 for query kmers manually]\n", argv[0]);
        exit(1);
    }
}

void InitCharMap(uint8_t *char_map) {
    fill(char_map, char_map + 256, 0);
    char_map['A'] = 1;
    char_map['C'] = 2;
    char_map['G'] = 3;
    char_map['T'] = 4;
}

bool Query(SuccinctDBG &dbg, string &kmer, uint8_t *char_map, bool print_message) {
    static char acgt[] = "$ACGT";
    uint8_t seq[dbg.kMaxKmerK];

    for (unsigned i = 0; i < kmer.length(); ++i) {
        seq[i] = char_map[(int)kmer[i]];
    }

    int64_t x = dbg.IndexBinarySearch(seq);

    if (x == -1) {
        if (print_message) cout << "Can not index this kmer." << endl;

        return false;
    }

    if (print_message) {
        int64_t edge[4];
        int edge_mul[4];
        cout << "Node index: " << x << endl;
        cout << "Indegree: " << dbg.Indegree(x) << endl;
        int out_degree = dbg.Outgoings(x, edge, edge_mul);
        cout << "Outdegree: " << out_degree << endl;

        for (int i = 0; i < out_degree; ++i) {
            cout << acgt[dbg.GetNodeLastChar(edge[i])] << ": " << edge_mul[i] << endl;
        }

        cout << "NodeMultiplicity: " << dbg.NodeMultiplicity(x) << endl;

        {
            // construct its unitigs
            string unitig;
            int64_t curr_node = x;
            int64_t prev, next;

            while ((prev = assembly_algorithms::PrevSimplePathNode(dbg, curr_node)) != -1) {
                unitig.push_back(acgt[dbg.GetW(prev)]);
                curr_node = prev;
            }

            dbg.Label(curr_node, seq);

            for (int k = dbg.kmer_k - 1; k >= 0; --k) {
                unitig.push_back(acgt[seq[k]]);
            }

            reverse(unitig.begin(), unitig.end());
            curr_node = x;

            while ((next = assembly_algorithms::NextSimplePathNode(dbg, curr_node)) != -1) {
                unitig.push_back(acgt[dbg.GetW(curr_node)]);
                curr_node = next;
            }

            cout << unitig << endl;
        }

        int64_t y = dbg.ReverseComplement(x);
        cout << "Reverse Complement:" << endl;
        cout << "Node index: " << y << endl;
        cout << "Indegree: " << dbg.Indegree(y) << endl;
        out_degree = dbg.Outgoings(y, edge, edge_mul);
        cout << "Outdegree: " << out_degree << endl;

        for (int i = 0; i < out_degree; ++i) {
            cout << acgt[dbg.GetNodeLastChar(edge[i])] << ": " << edge_mul[i] << endl;
        }

        cout << "NodeMultiplicity: " << dbg.NodeMultiplicity(y) << endl;
    }

    return true;
}

int main(int argc, char **argv) {
    CheckAndPrintUsage(argc, argv);
    const char *dbg_name = argv[1];
    int mode = atoi(argv[2]);
    SuccinctDBG dbg;

    {
        // graph loading
        printf("Loading succinct de Bruijn graph: %s\n", dbg_name);
        dbg.LoadFromFile(dbg_name);
        printf("Number of Edges: %ld\n", dbg.size);;
        printf("K value: %d\n", dbg.kmer_k);
    }

    uint8_t char_map[256];
    InitCharMap(char_map);

    int64_t bad_list[] = { 58418278, 602975976, 1033673088, 172278868, 775254816, 3 };

    for (int i = 0; i < 6; ++i) {
        printf("Node: %ld, multi: %d\n", bad_list[i], dbg.NodeMultiplicity(bad_list[i]));
        printf("Rev_node: %ld, multi: %d\n", dbg.ReverseComplement(bad_list[i]), dbg.NodeMultiplicity(dbg.ReverseComplement(bad_list[i])));
        char acgt[] = "$ACGT";
        uint8_t seq[dbg.kMaxKmerK];
        dbg.Label(bad_list[i], seq);

        for (int j = 0; j < dbg.kmer_k; ++j) {
            printf("%c", acgt[seq[j]]);
        }

        puts("");
    }

    if (false) {
        omp_lock_t lock;
        omp_init_lock(&lock);
        #pragma omp parallel for

        for (int64_t i = 0; i < dbg.size; ++i) {
            if (dbg.IsValidNode(i) && dbg.IsLast(i)) {
                if (dbg.NodeMultiplicity(i) != dbg.NodeMultiplicity(dbg.ReverseComplement(i))) {
                    omp_set_lock(&lock);
                    printf("Node: %ld, multi: %d\n", i, dbg.NodeMultiplicity(i));
                    printf("Rev_node: %ld, multi: %d\n", dbg.ReverseComplement(i), dbg.NodeMultiplicity(dbg.ReverseComplement(i)));
                    omp_unset_lock(&lock);
                }
            }
        }

        omp_destroy_lock(&lock);
        puts("Evaluate NodeMultiplicity done!");
    }

    string kmer;
    bool print_message = mode != 1;
    int num_found = 0;
    int num_total = 0;
    vector<string> vs;

    if (mode == 2) {
        while (cin >> kmer) {
            if (kmer == "Remove") {
                int len;
                cin >> len;
                assembly_algorithms::Trim(dbg, len, 200);
                continue;
            }

            if ((int)kmer.length() != dbg.kmer_k) {
                cout << "k value not match!" << endl;
                continue;
            }
            else {
                ++num_total;
                num_found += Query(dbg, kmer, char_map, print_message);
            }
        }
    }
    else {
        while (cin >> kmer) {
            if ((int)kmer.length() != dbg.kmer_k) {
                cout << "k value not match!" << endl;
                continue;
            }

            vs.push_back(kmer);
        }

        num_total = vs.size();
        omp_lock_t lock;
        omp_init_lock(&lock);
        #pragma omp parallel for

        for (unsigned i = 0; i < vs.size(); ++i) {
            if (Query(dbg, vs[i], char_map, print_message)) {
                #pragma omp atomic
                ++num_found;
            }
            else {
                omp_set_lock(&lock);
                cerr << vs[i] << endl;
                omp_unset_lock(&lock);
            }
        }

        omp_destroy_lock(&lock);
    }

    cout << "Total kmer: " << num_total << endl;
    cout << "Found kmer: " << num_found << endl;

    return 0;
}