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

#ifndef MISC_TIMER_H_
#define MISC_TIMER_H_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <fcntl.h>
#include <unistd.h>

struct xtimer_t {
    struct timeval tv1, tv2;
    long long time_elapsed;

    xtimer_t() { reset(); }
    void reset() { time_elapsed = 0; }
    void start() { gettimeofday(&tv1, NULL); }
    void stop() {
        gettimeofday(&tv2, NULL);
        time_elapsed += (long long)(tv2.tv_sec - tv1.tv_sec) * 1000000 + tv2.tv_usec - tv1.tv_usec;
    }
    double elapsed() { return time_elapsed / 1000000.0; }
};

struct AutoMaxRssRecorder {
    struct timeval tv1, tv2;

    AutoMaxRssRecorder() {
        gettimeofday(&tv1, NULL);
    }

    ~AutoMaxRssRecorder() {
        watch();
    }

    void watch() {
#define TURN_ON_MAX_RSS_LOG
#ifdef TURN_ON_MAX_RSS_LOG
        gettimeofday(&tv2, NULL);
        struct rusage usage;
        if (getrusage(RUSAGE_SELF, &usage)) {
            fprintf(stderr, "Fail to getrusage()\n");
            return;
        }
        double utime = 1e-6 * usage.ru_utime.tv_usec + usage.ru_utime.tv_sec;
        double stime = 1e-6 * usage.ru_stime.tv_usec + usage.ru_stime.tv_sec;

        long long real_time = (long long)(tv2.tv_sec - tv1.tv_sec) * 1000000 + tv2.tv_usec - tv1.tv_usec;
        fprintf(stderr, "Real: %.4lf", real_time / 1000000.0);
        fprintf(stderr, "\tuser: %.4lf\tsys: %.4lf\tmaxrss: %ld\n", utime, stime, usage.ru_maxrss);
#endif     
    }
};

#endif