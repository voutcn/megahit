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

#ifndef DBG_RANK_AND_SELECT_H_
#define DBG_RANK_AND_SELECT_H_

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

const int kBitsPerByte = 8;
const int kBitsPerULL = sizeof(unsigned long long) * kBitsPerByte;
typedef uint32_t interval_t;

class RankAndSelect4Bits {
  public:
    // constants
    static const int kCharPerInterval = 256;	// tunable
    static const int kSelectSampleSize = 256;	// tunable
    static const int kAlphabetSize = 9;
    static const int kBitsPerChar = 4;
    static const int kCharPerWord = sizeof(unsigned long long) * kBitsPerByte / kBitsPerChar;
    static const int kCharPerIntervalMajor = 65536;
    static const int kMinorPerMajor = kCharPerIntervalMajor / kCharPerInterval;
    // public data, can call directly
    int64_t length;
    int64_t char_frequency[kAlphabetSize];

    RankAndSelect4Bits();
    ~RankAndSelect4Bits();

    void Build(unsigned long long *packed_text, int64_t length);	// initialize
    int64_t Rank(uint8_t c, int64_t pos);	// the number of c's in [0...pos]
    int64_t Select(uint8_t c, int64_t ranking); // return the pos of the ranking_th c (0-based)
    int64_t Pred(uint8_t c, int64_t pos);	// the last c in [0...pos]
    int64_t PredLimitedStep(uint8_t c, int64_t pos, int step); // the last c in [pos-step, pos], return pos-step-1 if not exist
    int64_t Succ(uint8_t c, int64_t pos);	// the first c in [pos...length]
    int64_t SuccLimitedStep(uint8_t c, int64_t pos, int step); // the first c in [pos, pos+step], return pos+step+1 if not exist

  private:
    int CountCharInWord_(uint8_t c, unsigned long long x, unsigned long long mask) {
        x ^= popcount_char_xorer_[c];
        x &= x >> 2;
        x &= x >> 1;
        return __builtin_popcountll(x & popcount_mask_ & mask);
    }

    int CountCharInWord_(uint8_t c, unsigned long long x) {
        x ^= popcount_char_xorer_[c];
        x &= x >> 2;
        x &= x >> 1;
        return __builtin_popcountll(x & popcount_mask_);
    }

    int SelectInWord_(uint8_t c, int num_c, unsigned long long x) {
        int tailing_zero = 0;
        x ^= popcount_char_xorer_[c];
        x &= x >> 2;
        x &= x >> 1;
        x &= popcount_mask_;

        while (num_c > 0) {
            tailing_zero = __builtin_ctzll(x);
            x ^= 1ULL << tailing_zero;
            --num_c;
        }
        return tailing_zero / kBitsPerChar; // 0-based
    }

    void PrefectchOccValue_(uint8_t c, int64_t i) {
        __builtin_prefetch(occ_value_explicit_major_[c] + i / kMinorPerMajor, 0);
        __builtin_prefetch(occ_value_explicit_minor_[c] + i, 0);
    }

    int64_t OccValue_(uint8_t c, int64_t i) {
        return occ_value_explicit_major_[c][i / kMinorPerMajor] +
               occ_value_explicit_minor_[c][i];
    }

  private:
    // main memory for the structure
    unsigned long long *packed_text_;

    // sampled structure for rank and select
    // two level sampling for rank (occ value)
    // call the function OccValue_(c, i) to get the number of c's in packed_text_[0...i*kCharPerInterval-1]
    int64_t *occ_value_explicit_major_[kAlphabetSize];
    uint16_t *occ_value_explicit_minor_[kAlphabetSize];

    // sampling for select
    // rank_to_interval_explicit_[c][i]=j: the jth interval (0 based) contains the (i*kSelectSampleSize)th (0 based) c
    // i.e. OccValue_(c, j)<=i*kSelectSampleSize and OccValue_(c, j+1)>i*kSelectSampleSize
    interval_t *rank_to_interval_explicit_[kAlphabetSize];

    // popcount masks
    unsigned long long popcount_char_xorer_[kAlphabetSize]; // e.g. if c = 0110(2), popcount_char_xorer_[kAlphabetSize] = 1001 1001 1001 1001...(2), to make all c's in a word 1111
    unsigned long long popcount_mask_; // e.g. for kBitsPerChar=4, popcount_mask_ = 0x1111111111111111ULL
};

class RankAndSelect1Bit {
  public:
    static const int kBitsPerInterval = 256;	// tunable
    static const int kSelectSampleSize = 256;	// tunable
    static const int kBitsPerWord = sizeof(unsigned long long) * kBitsPerByte;
    static const int kBitsPerMajorInterval = 65536;
    static const int kMinorPerMajor = kBitsPerMajorInterval / kBitsPerInterval;
    int64_t length;
    int64_t total_num_ones;

    RankAndSelect1Bit();
    ~RankAndSelect1Bit();

    void Build(unsigned long long *packed_text, int64_t length);
    int64_t Rank(int64_t pos);
    int64_t Select(int64_t ranking);
    int64_t Pred(int64_t pos);
    int64_t Succ(int64_t pos);

  private:
    void PrefectchOccValue_(int64_t i) {
        __builtin_prefetch(occ_value_explicit_major_ + i / kMinorPerMajor, 0);
        __builtin_prefetch(occ_value_explicit_minor_ + i, 0);
    }

    int64_t OccValue_(int64_t i) {
        return occ_value_explicit_major_[i / kMinorPerMajor] +
               occ_value_explicit_minor_[i];
    }

    int SelectInWord_(int num, unsigned long long x) {
        int tailing_zero = 0;
        while (num > 0) {
            tailing_zero = __builtin_ctzll(x);
            x ^= 1LL << tailing_zero;
            --num;
        }
        return tailing_zero; // 0-based
    }

  private:
    unsigned long long *packed_text_;
    int64_t *occ_value_explicit_major_;
    uint16_t *occ_value_explicit_minor_;
    uint32_t *rank_to_interval_explicit_;

};

#endif  // DBG_RANK_AND_SELECT_H_