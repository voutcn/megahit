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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>
#include <assert.h>

using namespace std;

#define READLEN 100
#define WORDS_PER_READ 7
#define MOD_WPR 127
#define DIV_WPR 7

#define BUFFER 8192

typedef long long int int64;

int main(int argc, char **argv) {


    char buff[1000];
    while (fgets(buff, 1000, stdin) != NULL) {
        printf(">XX\n");
        fgets(buff, 1000, stdin);
        buff[100] = 0;
        printf("%s", buff);
        fgets(buff, 1000, stdin);
        fgets(buff, 1000, stdin);
        printf("%s", buff);
    }

    return 0;

}

