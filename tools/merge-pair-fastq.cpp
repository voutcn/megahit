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

using namespace std;

int main(int argc, char **argv) {

    FILE *f1, *f2;

    f1 = fopen( argv[1], "r" );
    f2 = fopen( argv[2], "r" );

    char buff1[1000], buff2[1000];

    while ( fgets( buff1, 1000, f1 ) ) {
        printf("%s", buff1);
        fgets( buff1, 1000, f1 );
        printf("%s", buff1);
        fgets( buff1, 1000, f1 );
        printf("%s", buff1);
        fgets( buff1, 1000, f1 );
        printf("%s", buff1);

        fgets( buff2, 1000, f2 );
        printf("%s", buff2);
        fgets( buff2, 1000, f2 );
        printf("%s", buff2);
        fgets( buff2, 1000, f2 );
        printf("%s", buff2);
        fgets( buff2, 1000, f2 );
        printf("%s", buff2);

    }

    fclose(f1);
    fclose(f2);

    return 0;
}
