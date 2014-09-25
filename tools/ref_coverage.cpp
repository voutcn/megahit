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

#include <map>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

using namespace std;

int main(int argc, char** argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s [ref-id_start_end_file] [ref-len_file]\n", argv[0]);
        exit(1);
    }

    ifstream align_file(argv[1]);
    ifstream len_file(argv[2]);
    map<string, vector<pair<long long, long long> > > cov_map;
    map<string, long long> len_map;
    string id;
    long long start, end, len;
    long long total_align_len = 0;

    while (align_file >> id) {
        align_file >> start >> end;
        if (start > end) {
            swap(start, end);
        }
        cov_map[id].push_back(make_pair(start, end));
        total_align_len += end - start + 1;
    }

    cout << "Total number of segments: " << cov_map.size() << endl;
    cout << "Total aligned length: " << total_align_len << endl;

    while (len_file >> id) {
        len_file >> len;
        len_map[id] = len;
    }

    cout << "Total number of ref_seqs: " << len_map.size() << endl;

    long long total_cov_len = 0;
    vector<int> cov_level_count(100);
    for (auto it = cov_map.begin(); it != cov_map.end(); ++it) {
        sort(it->second.begin(), it->second.end());
        long long cur_end = 0;
        long long cov_len = 0;
        long long rep_len = 0;

        for (auto jt = it->second.begin(); jt != it->second.end(); ++jt) {
            if (cur_end >= jt->first) {
                if (cur_end <= jt->second) {
                    rep_len += cur_end - jt->first + 1;
                } else {
                    rep_len += jt->second - jt->first + 1;
                }
                jt->first = cur_end + 1;
            }
            if (jt->second >= jt->first) {
                cov_len += jt->second - jt->first + 1;
            }

            cur_end = max(cur_end, jt->second);
        }
        int cov_level = cov_len * 1.0 / len_map[it->first] * 100;
        total_cov_len += cov_len;
        cov_level_count[cov_level]++;
        if (cov_level >= 5) {
            cout << "ref_id " << it->first << " coverage >= 50%, ref length " << len_map[it->first] << endl;
        }
        if (true || rep_len > cov_len) {
            // cerr << "ref_id " << it->first << " coverage by many redundant contigs (" <<  rep_len << "/" << cov_len << "/" << len_map[it->first] << ")" << endl;
        }
    }

    for (int i = 0; i <= 100; ++i) {
        cout << "Refs with cov no less than " << i * 1.0 / 100 << " " << cov_level_count[i] << endl;
    }

    cout << "Total covered length: " << total_cov_len << endl;

    return 0;
}