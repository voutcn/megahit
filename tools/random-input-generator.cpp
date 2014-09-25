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
#include <time.h>

int main(int argc, char **argv) {

    srand(time(0));
    long long int n;
    int len;

    if (argc != 3) {
        fprintf(stderr, "Usage: %s num_reads read_length\n", argv[0]);
        return 0;
    }
    n = atoll(argv[1]);
    len = atoi(argv[2]);

    char *buff = (char*)malloc(sizeof(char) * (len*2+50));

    buff[0] = '@'; // >
    buff[1] = 'X';
    buff[2] = '\n';
    buff[len+3] = '\n';
    buff[len+4] = '+';
    buff[len+5] = '\n';
    for (int i=len+6; i<len+6+len; ++i)
        buff[i] = 'g';
    buff[len+6+len] = '\n';
    buff[len+7+len] = 0;
    //  buff[len+4] = 0;

    char map[4] = {'A','C','G','T'};

    while (n--) {
        for (int i=3; i<len+3; ++i)
            buff[i] = map[rand()%4];
        printf("%s", buff);
    }

    free(buff);

    return 0;
}
