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

#ifndef MEGAHIT_UTILS_H__
#define MEGAHIT_UTILS_H__

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <fcntl.h>
#include <unistd.h>

static const int kSignBitMask = 0x80000000; // the MSB of 32-bit

template<typename T1, typename T2>
inline T1 DivCeiling(T1 a, T2 b) {
    return (a + b - 1) / b;
}

inline void megahit_log__(const char *format, ...) {
    static volatile int lock_ = 0;
    while (__sync_lock_test_and_set(&lock_, 1)) while (lock_);
    va_list args;
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    __sync_lock_release(&lock_);
}

#define xlog_ext(str, args...) megahit_log__(str, ##args);
#define xlog(str, args...) megahit_log__("    [%-26s:%4d]     ", __FILE__, __LINE__); megahit_log__(str, ##args);
#define xerr(str, args...) megahit_log__("    [ERROR] [%-25s:%4d]: ", __FILE__, __LINE__); megahit_log__(str, ##args);
#define xwarning(str, args...) megahit_log__("    [WARNING] [%-25s:%4d] ", __FILE__, __LINE__); megahit_log__(str, ##args);
#define xerr_and_exit(str, args...) megahit_log__("    [ERROR] [%-25s:%4d]: ", __FILE__, __LINE__); megahit_log__(str, ##args); exit(1);

#ifdef __GNUC__
    #define LIKELY(x) __builtin_expect((x),1)
    #define UNLIKELY(x) __builtin_expect((x),0)
#else
    #define LIKELY(x) (x)
    #define UNLIKELY(x) (x)
#endif

inline char *FormatString(const char *fmt, ...) {
    static char buffer[1 << 20];
    va_list args;
    va_start(args, fmt);
    vsprintf(buffer, fmt, args);
    va_end(args);
    return buffer;
}

inline double time_elapsed(struct timeval &tv_start) {
    struct timeval tv_now;
    gettimeofday(&tv_now, NULL);
    return (long long)(tv_now.tv_sec - tv_start.tv_sec) * 1000000 + tv_now.tv_usec - tv_start.tv_usec;
}

struct xtimer_t {
    struct timeval tv1, tv2;
    long long time_elapsed;

    xtimer_t() {
        reset();
    }
    void reset() {
        time_elapsed = 0;
    }
    void start() {
        gettimeofday(&tv1, NULL);
    }
    void stop() {
        gettimeofday(&tv2, NULL);
        time_elapsed += (long long)(tv2.tv_sec - tv1.tv_sec) * 1000000 + tv2.tv_usec - tv1.tv_usec;
    }
    double elapsed() {
        return time_elapsed / 1000000.0;
    }
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
            xwarning("Fail to getrusage()\n");
        }

        double utime = 1e-6 * usage.ru_utime.tv_usec + usage.ru_utime.tv_sec;
        double stime = 1e-6 * usage.ru_stime.tv_usec + usage.ru_stime.tv_sec;

        long long real_time = (long long)(tv2.tv_sec - tv1.tv_sec) * 1000000 + tv2.tv_usec - tv1.tv_usec;
        xlog("Real: %.4lf", real_time / 1000000.0);
        xlog_ext("\tuser: %.4lf\tsys: %.4lf\tmaxrss: %ld\n", utime, stime, usage.ru_maxrss);
#endif
    }
};

#endif // MEGAHIT_UTILS_H__