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
#include <vector>
#include <cstdlib>
#include "rank_and_select.h"
#include "misc.h"

using namespace std;

template<int ALPHA_SIZE = 16, int BIT_PER_CHAR = 4, int CHAR_PER_WORD = 16>
void testRankAndSelect4Bits(unsigned long long NUM_CHAR, unsigned long long REPEAT) {

    xtimer_t timer;
    unsigned long long *v = (unsigned long long*) malloc(sizeof(unsigned long long) * NUM_CHAR / CHAR_PER_WORD);
    unsigned long long a = 0;

    printf("[START] Generating %llu characters, ALPHA_SIZE: %d.\n", NUM_CHAR, ALPHA_SIZE);
    timer.reset();
    timer.start();
    for (unsigned long long i = 0; i < NUM_CHAR / CHAR_PER_WORD; ++i) {
        v[i] = rand();
        v[i] |= (unsigned long long)rand() << 32;
    }
    timer.stop();
    printf("[DONE] Generated. Time :%lf\n", timer.elapsed());

    vector<char> random_c(REPEAT);
    vector<unsigned long long> random_ll(REPEAT);

    for (unsigned long long i = 0; i < REPEAT; ++i) {
        random_c[i] = rand() % ALPHA_SIZE;
        random_ll[i] = (unsigned long long)rand() * rand() % NUM_CHAR;
    }

    printf("[START] Building rank&select\n");
    timer.reset();
    timer.start();
    RankAndSelect4Bits rsGeneral4Bits;
    rsGeneral4Bits.Build(v, NUM_CHAR);
    timer.stop();
    printf("[DONE] Built rand&select. Time: %lf\n", timer.elapsed());

    printf("[START] Testing rank, number of calling: %llu\n", REPEAT);
    timer.reset();
    timer.start();
    a = 0;
    for (int t = 0; t < REPEAT; ++t) {
        a += rsGeneral4Bits.Rank(random_c[t], random_ll[t]);
    }
    timer.stop();
    printf("[DONE] End testing rank. Time: %lf\n", timer.elapsed());
    printf("%llu\n", a);

    for (unsigned long long i = 0; i < REPEAT; ++i) {
        random_ll[i] = random_ll[i] % rsGeneral4Bits.char_frequency[random_c[i]] + 1;
    }

    printf("[START] Testing select, number of calling: %llu\n", REPEAT);
    timer.reset();
    timer.start();
    a = 0;
    for (unsigned long long t = 0; t < REPEAT; ++t) {
        a += rsGeneral4Bits.Select(random_c[t], random_ll[t]);
    }
    timer.stop();
    printf("[DONE] End testing select. Time: %lf\n", timer.elapsed());
    printf("%llu\n", a);

    free(v);
}

void testRankAndSelect1Bit(long long NUM_CHAR, long long REPEAT) {
    xtimer_t timer;
    unsigned long long *v = (unsigned long long*) malloc(sizeof(unsigned long long) * NUM_CHAR / 64);
    unsigned long long a = 0;

    printf("[START] Generating %llu characters.\n", NUM_CHAR);
    timer.reset();
    timer.start();
    for (unsigned long long i = 0; i < NUM_CHAR / 64; ++i) {
        v[i] = rand();
        v[i] |= (unsigned long long)rand() << 32;
    }
    timer.stop();
    printf("[DONE] Generated. Time :%lf\n", timer.elapsed());

    vector<unsigned long long> random_ll(REPEAT);
    for (int i = 0; i < REPEAT; ++i) {
        random_ll[i] = (unsigned long long)rand() * rand() % NUM_CHAR;
    }

    printf("[START] Building rank&select_1bit\n");
    timer.reset();
    timer.start();
    RankAndSelect1Bit rs01;
    rs01.Build(v, NUM_CHAR);
    timer.stop();
    printf("[DONE] Built rand&select_1bit. Time: %lf\n", timer.elapsed());

    a = 0;
    printf("[START] Testing rank_1bit, number of calling: %llu\n", REPEAT);
    timer.reset();
    timer.start();
    for (int t = 0; t < REPEAT; ++t) {
        a += rs01.Rank(random_ll[t]);
    }
    timer.stop();
    printf("[DONE] End testing rank_1bit. Time: %lf\n", timer.elapsed());
    printf("%llu\n", a);

    for (int i = 0; i < REPEAT; ++i) {
        random_ll[i] = random_ll[i] % rs01.total_num_ones + 1;
    }

    a = 0;
    printf("[START] Testing select_1bit, number of calling: %llu\n", REPEAT);
    timer.reset();
    timer.start();
    for (int t = 0; t < REPEAT; ++t) {
        a += rs01.Select(random_ll[t]);
    }
    timer.stop();
    printf("[DONE] End testing select_1bit. Time: %lf\n", timer.elapsed());
    printf("%llu\n", a);

    free(v);
}

int main(int argc, char **argv) {
    int numOfG = 3;
    if (argc > 1) {
        numOfG = atoi(argv[1]);
    }
    printf("%d\n", RAND_MAX);
    testRankAndSelect4Bits(1073741824LL * numOfG, 1000000000);
    testRankAndSelect1Bit(1073741824LL * numOfG, 1000000000);
    return 0;
}