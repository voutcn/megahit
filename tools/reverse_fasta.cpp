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
#include <string>
#include <algorithm>
#include <vector>
using namespace std;

int main() {
    string buf;
    vector<string> a[2];
    while (getline(cin, buf)) {
        a[0].push_back(buf);
        getline(cin, buf);
        a[1].push_back(buf);
        reverse(a[1].back().begin(), a[1].back().end());
    }

    for (int i = a[0].size() - 1; i >= 0; --i) {
        cout << a[0][i] << '\n';
        cout << a[1][i] << '\n';
    }
    return 0;
}
