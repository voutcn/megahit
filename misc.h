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

#ifndef DBJ_MISC_H_
#define DBJ_MISC_H_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <fcntl.h>
#include <unistd.h>
#include <string>

typedef long long int int64;
typedef unsigned int word;
typedef unsigned long long int uint64;

struct xtimer_t {
    struct timeval tv1, tv2;
    int64 time_elapsed;
    void reset() {
        time_elapsed = 0;
    }
    void start() {
        gettimeofday(&tv1, NULL);
    }
    void stop() {
        gettimeofday(&tv2, NULL);
        time_elapsed += (int64)(tv2.tv_sec - tv1.tv_sec) * 1000000 + tv2.tv_usec - tv1.tv_usec;
    }
    double elapsed() {
        return time_elapsed / 1000000.0;
    }
};

#endif