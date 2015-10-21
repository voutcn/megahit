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

#ifndef DBG_RANK_AND_SELECT_H_
#define DBG_RANK_AND_SELECT_H_

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

const int kBitsPerByte = 8;
const int kBitsPerULL = sizeof(unsigned long long) * kBitsPerByte;
typedef uint32_t interval_t;
#define DIFF_TO_DO_BINARY_SEARCH 2

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

    RankAndSelect4Bits() {
        for (int i = 0; i < kAlphabetSize; ++i) {
            occ_value_explicit_major_[i] = NULL;
            occ_value_explicit_minor_[i] = NULL;
            rank_to_interval_explicit_[i] = NULL;
            popcount_char_xorer_[i] = 0;

            for (int j = 0; j < kCharPerWord; ++j) {
                popcount_char_xorer_[i] |= (unsigned long long)i << (kBitsPerChar * j);
            }

            popcount_char_xorer_[i] = ~popcount_char_xorer_[i];
        }

        popcount_mask_ = ~popcount_char_xorer_[1];
    }

    ~RankAndSelect4Bits() {
        for (int i = 0; i < kAlphabetSize; ++i) {
            if (occ_value_explicit_major_[i] != NULL) {
                free(occ_value_explicit_major_[i]);
            }

            if (occ_value_explicit_minor_[i] != NULL) {
                free(occ_value_explicit_minor_[i]);
            }

            if (rank_to_interval_explicit_[i] != NULL) {
                free(rank_to_interval_explicit_[i]);
            }
        }
    }

    void Build(unsigned long long *packed_text, int64_t length) {
        int64_t count[kAlphabetSize];
        int64_t num_intervals = (length + kCharPerInterval - 1) / kCharPerInterval + 1;
        int64_t num_intervals_major = (length + kCharPerIntervalMajor - 1) / kCharPerIntervalMajor + 1;
        unsigned long long *cur_word = packed_text;

        // build rank
        for (int i = 0; i < kAlphabetSize; ++i) {
            count[i] = 0;
            occ_value_explicit_major_[i] = (int64_t *) malloc(sizeof(int64_t) * num_intervals_major);

            if (occ_value_explicit_major_[i] == NULL) {
                fprintf(stderr, "Malloc Failed: %s: %d\n", __FILE__, __LINE__);
                exit(1);
            }

            occ_value_explicit_minor_[i] = (uint16_t *) malloc(sizeof(uint16_t) * num_intervals);

            if (occ_value_explicit_minor_[i] == NULL) {
                fprintf(stderr, "Malloc Failed: %s: %d\n", __FILE__, __LINE__);
                exit(1);
            }
        }

        for (int64_t i = 0; i < length; i += kCharPerWord, ++cur_word) {
            if (i % kCharPerInterval == 0) {
                for (int j = 0; j < kAlphabetSize; ++j) {
                    if (i % kCharPerIntervalMajor == 0) {
                        occ_value_explicit_major_[j][i / kCharPerIntervalMajor] = count[j];
                    }

                    occ_value_explicit_minor_[j][i / kCharPerInterval] = count[j] - occ_value_explicit_major_[j][i / kCharPerIntervalMajor];
                }
            }

            for (int j = 0; j < kAlphabetSize; ++j) {
                count[j] += CountCharInWord_(j, *cur_word);
            }
        }

        for (int j = 0; j < kAlphabetSize; ++j) {
            occ_value_explicit_major_[j][num_intervals_major - 1] = count[j];
            occ_value_explicit_minor_[j][num_intervals - 1] = count[j] - occ_value_explicit_major_[j][(num_intervals - 1) / kMinorPerMajor];
            char_frequency[j] = count[j];
        }

        // build select look up table
        // rank_to_interval_explicit_[c][i]=j: the jth interval (0 based) contains the (i*kSelectSampleSize)th (1 based) c
        // i.e. OccValue_(c,j)<=i*kSelectSampleSize and OccValue_(c,j+1)>i*kSelectSampleSize

        for (int j = 0; j < kAlphabetSize; ++j) {
            interval_t s_table_size = (count[j] + kSelectSampleSize - 1) / kSelectSampleSize + 1;
            rank_to_interval_explicit_[j] = (interval_t *) malloc(sizeof(interval_t) * s_table_size);

            if (rank_to_interval_explicit_[j] == NULL) {
                fprintf(stderr, "Malloc Failed: %s: %d\n", __FILE__, __LINE__);
                exit(1);
            }

            for (int64_t i = 0, s_table_idx = 0; i < num_intervals; ++i) {
                while (s_table_idx * kSelectSampleSize < OccValue_(j, i)) {
                    rank_to_interval_explicit_[j][s_table_idx] = i - 1;
                    ++s_table_idx;
                }
            }

            rank_to_interval_explicit_[j][s_table_size - 1] = num_intervals - 1;
        }

        packed_text_ = packed_text;
        this->length = length;
    }

    int64_t Rank(uint8_t c, int64_t pos) {
        // the number of c's in [0...pos]
        if (pos >= length - 1) {
            return char_frequency[c];
        }

        unsigned long long mask;
        unsigned long long *first_word;
        int total_chars_to_count;
        int words_to_count;
        int chars_to_count;
        int count_c = 0;

        ++pos;
        int64_t which_interval = (pos + kCharPerInterval / 2 - 1) / kCharPerInterval;
        int64_t sampled_index = which_interval * kCharPerInterval;

        if (sampled_index >= length) {
            sampled_index -= kCharPerInterval;
            which_interval--;
        }

        PrefectchOccValue_(c, which_interval);

        if (sampled_index > pos) {
            total_chars_to_count = sampled_index - pos;
            words_to_count = total_chars_to_count / kCharPerWord;
            first_word = packed_text_ + sampled_index / kCharPerWord - words_to_count - 1;
            __builtin_prefetch(first_word);
            chars_to_count = total_chars_to_count % kCharPerWord;

            if (chars_to_count > 0) {
                mask = -(1ULL << kBitsPerChar * (kCharPerWord - chars_to_count));
                count_c += CountCharInWord_(c, first_word[0], mask);
            }

            for (int i = 1; i <= words_to_count; ++i) {
                count_c += CountCharInWord_(c, first_word[i]);
            }

            return OccValue_(c, which_interval) - count_c;

        }
        else if (sampled_index < pos) {
            first_word = packed_text_ + sampled_index / kCharPerWord;
            __builtin_prefetch(first_word);
            total_chars_to_count = pos - sampled_index;
            words_to_count = total_chars_to_count / kCharPerWord;
            chars_to_count = total_chars_to_count % kCharPerWord;

            for (int i = 0; i < words_to_count; ++i) {
                count_c += CountCharInWord_(c, first_word[i]);
            }

            if (chars_to_count > 0) {
                mask = (1ULL << kBitsPerChar * chars_to_count) - 1;
                count_c += CountCharInWord_(c, first_word[words_to_count], mask);
            }

            return OccValue_(c, which_interval) + count_c;

        }
        else {
            return OccValue_(c, which_interval);
        }
    }

    int64_t Select(uint8_t c, int64_t ranking) {
        // return the pos of the ranking_th c (0-based)
        if (ranking >= char_frequency[c]) {
            return length;
        }
        else if (ranking < 0) {
            return -1;
        }

        // first locate which interval Select(c, ranking) falls
        interval_t interval_l = rank_to_interval_explicit_[c][ranking / kSelectSampleSize];
        interval_t interval_r = rank_to_interval_explicit_[c][(ranking + kSelectSampleSize - 1) / kSelectSampleSize];
        interval_t interval_m;

        while (interval_r > interval_l + DIFF_TO_DO_BINARY_SEARCH) {
            interval_m = (interval_r + interval_l + 1) / 2;

            // PrefectchOccValue_(c, (interval_m + interval_l) / 2);
            // PrefectchOccValue_(c, (interval_m + interval_r + 1) / 2);
            if (OccValue_(c, interval_m) > ranking) {
                interval_r = interval_m - 1;
            }
            else {
                interval_l = interval_m;
            }
        }

#if DIFF_TO_DO_BINARY_SEARCH > 0
        PrefectchOccValue_(c, interval_l);

        if (interval_r > interval_l) {
            while (OccValue_(c, interval_l + 1) <= ranking) {
                ++interval_l;
            }
        }

#endif

        int64_t pos = (int64_t)interval_l * kCharPerInterval;
        unsigned long long *cur_word = packed_text_ + pos / kCharPerWord;
        int pos_in_word = 0;
        __builtin_prefetch(cur_word);

        int remaining_c = ranking + 1 - OccValue_(c, interval_l);
        int popcnt;

        for (; ; pos_in_word += kCharPerWord) {
            popcnt = CountCharInWord_(c, *cur_word);

            if (popcnt >= remaining_c) {
                break;
            }
            else {
                remaining_c -= popcnt;
            }

            ++cur_word;
        }

        return pos + pos_in_word + SelectInWord_(c, remaining_c, *cur_word);
    }

    int64_t Pred(uint8_t c, int64_t pos) {
        // the last c in [0...pos]
        if (((*(packed_text_ + pos / kCharPerWord) >> (pos % kCharPerWord * kBitsPerChar)) & ((1 << kBitsPerChar) - 1)) == c) {
            return pos;
        }

        return Select(c, Rank(c, pos) - 1);
    }

    int64_t PredLimitedStep(uint8_t c, int64_t pos, int step) {
        // the last c in [pos-step, pos], return pos-step-1 if not exist
        int64_t end = pos - step;

        if (end < 0) {
            end = 0;
        }

        while (pos >= end) {
            if (((*(packed_text_ + pos / kCharPerWord) >> (pos % kCharPerWord * kBitsPerChar)) & ((1 << kBitsPerChar) - 1)) == c) {
                return pos;
            }

            --pos;
        }

        return pos;
    }

    int64_t Succ(uint8_t c, int64_t pos) {
        // the first c in [pos...length]
        if (((*(packed_text_ + pos / kCharPerWord) >> (pos % kCharPerWord * kBitsPerChar)) & ((1 << kBitsPerChar) - 1)) == c) {
            return pos;
        }

        return Select(c, Rank(c, pos - 1));
    }

    int64_t SuccLimitedStep(uint8_t c, int64_t pos, int step) {
        // the first c in [pos, pos+step], return pos+step+1 if not exist
        int64_t end = pos - step;

        if (end >= length) {
            end = length;
        }

        while (pos <= end) {
            if (((*(packed_text_ + pos / kCharPerWord) >> (pos % kCharPerWord * kBitsPerChar)) & ((1 << kBitsPerChar) - 1)) == c) {
                return pos;
            }

            ++pos;
        }

        return pos;
    }

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

template <bool rank_only = false>
class RankAndSelect1Bit {
  public:
    static const int kBitsPerInterval = 256;	// tunable
    static const int kSelectSampleSize = 256;	// tunable
    static const int kBitsPerWord = sizeof(unsigned long long) * kBitsPerByte;
    static const int kBitsPerMajorInterval = 65536;
    static const int kMinorPerMajor = kBitsPerMajorInterval / kBitsPerInterval;
    int64_t length;
    int64_t total_num_ones;

    RankAndSelect1Bit() {
        occ_value_explicit_minor_ = NULL;
        occ_value_explicit_major_ = NULL;
        rank_to_interval_explicit_ = NULL;
    }

    ~RankAndSelect1Bit() {
        if (occ_value_explicit_major_ != NULL) {
            free(occ_value_explicit_major_);
        }

        if (occ_value_explicit_minor_ != NULL) {
            free(occ_value_explicit_minor_);
        }

        if (rank_to_interval_explicit_ != NULL) {
            free(rank_to_interval_explicit_);
        }
    }

    void Build(unsigned long long *packed_text, int64_t length) {
        int64_t count_ones = 0;
        int64_t num_intervals = (length + kBitsPerInterval - 1) / kBitsPerInterval + 1;
        int64_t num_intervals_major = (length + kBitsPerMajorInterval - 1) / kBitsPerMajorInterval + 1;
        unsigned long long *cur_word = packed_text;

        occ_value_explicit_major_ = (int64_t *) malloc(sizeof(int64_t) * num_intervals_major);

        if (occ_value_explicit_major_ == NULL) {
            fprintf(stderr, "Malloc Failed: %s: %d\n", __FILE__, __LINE__);
            exit(1);
        }

        occ_value_explicit_minor_ = (uint16_t *) malloc(sizeof(uint16_t) * num_intervals);

        if (occ_value_explicit_minor_ == NULL) {
            fprintf(stderr, "Malloc Failed: %s: %d\n", __FILE__, __LINE__);
            exit(1);
        }

        for (int64_t i = 0; i < length; i += kBitsPerWord, ++cur_word) {
            if (i % kBitsPerInterval == 0) {
                if (i % kBitsPerMajorInterval == 0) {
                    occ_value_explicit_major_[i / kBitsPerMajorInterval] = count_ones;
                }

                occ_value_explicit_minor_[i / kBitsPerInterval] = count_ones - occ_value_explicit_major_[i / kBitsPerMajorInterval];
            }

            count_ones += __builtin_popcountll(*cur_word);
        }

        occ_value_explicit_major_[num_intervals_major - 1] = count_ones;
        occ_value_explicit_minor_[num_intervals - 1] = count_ones - occ_value_explicit_major_[(num_intervals - 1) / kMinorPerMajor];
        total_num_ones = count_ones;

        // Build select look up table
        if (!rank_only) {
            uint32_t s_table_size = (count_ones + kSelectSampleSize - 1) / kSelectSampleSize + 1;
            rank_to_interval_explicit_ = (uint32_t *) malloc(sizeof(uint32_t) * s_table_size);

            if (rank_to_interval_explicit_ == NULL) {
                fprintf(stderr, "Malloc Failed: %s: %d\n", __FILE__, __LINE__);
                exit(1);
            }

            int64_t s_table_idx = 0;

            for (int64_t i = 0; i < num_intervals; ++i) {
                while (s_table_idx * kSelectSampleSize < OccValue_(i)) {
                    rank_to_interval_explicit_[s_table_idx] = i - 1;
                    ++s_table_idx;
                }
            }

            rank_to_interval_explicit_[s_table_size - 1] = num_intervals - 1;
        }

        packed_text_ = packed_text;
        this->length = length;
    }

    int64_t Rank(int64_t pos) {
        if (pos > length - 1) {
            return total_num_ones;
        }

        ++pos;
        unsigned long long mask;
        unsigned long long *first_word;
        int total_bits_to_count;
        int words_to_count;
        int bits_to_count;
        int count_ones = 0;

        int64_t which_interval = (pos + kBitsPerInterval / 2 - 1) / kBitsPerInterval;
        int64_t sampled_index = which_interval * kBitsPerInterval;

        if (sampled_index > length) {
            sampled_index -= kBitsPerInterval;
            which_interval--;
        }

        PrefectchOccValue_(which_interval);

        if (sampled_index > pos) {
            total_bits_to_count = sampled_index - pos;
            words_to_count = total_bits_to_count / kBitsPerWord;
            first_word = packed_text_ + sampled_index / kBitsPerWord - words_to_count - 1;
            __builtin_prefetch(first_word);

            bits_to_count = total_bits_to_count % kBitsPerWord;

            if (bits_to_count > 0) {
                mask = -(1ULL << (kBitsPerWord - bits_to_count));
                count_ones += __builtin_popcountll(first_word[0] & mask);
            }

            for (int i = 1; i <= words_to_count; ++i) {
                count_ones += __builtin_popcountll(first_word[i]);
            }

            return OccValue_(which_interval) - count_ones;

        }
        else if (sampled_index < pos) {
            first_word = packed_text_ + sampled_index / kBitsPerWord;
            __builtin_prefetch(first_word);

            total_bits_to_count = pos - sampled_index;
            words_to_count = total_bits_to_count / kBitsPerWord;
            bits_to_count = total_bits_to_count % kBitsPerWord;

            for (int i = 0; i < words_to_count; ++i) {
                count_ones += __builtin_popcountll(first_word[i]);
            }

            if (bits_to_count > 0) {
                mask = (1ULL << bits_to_count) - 1;
                count_ones += __builtin_popcountll(first_word[words_to_count] & mask);
            }

            return OccValue_(which_interval) + count_ones;

        }
        else {
            return OccValue_(which_interval);
        }
    }

    int64_t Select(int64_t ranking) {
        static_assert(rank_only == false, "cannot select in rank_only struct");

        if (ranking >= total_num_ones) {
            return length;
        }
        else if (ranking < 0) {
            return -1;
        }

        // first locate which interval Select(c, ranking) falls
        uint32_t interval_l = rank_to_interval_explicit_[ranking / kSelectSampleSize];
        uint32_t interval_r = rank_to_interval_explicit_[(ranking + kSelectSampleSize - 1) / kSelectSampleSize];
        uint32_t interval_m;

        while (interval_r > interval_l + DIFF_TO_DO_BINARY_SEARCH) {
            interval_m = (interval_r + interval_l + 1) / 2;

            // PrefectchOccValue_(c, (interval_m + interval_l) / 2);
            // PrefectchOccValue_(c, (interval_m + interval_r + 1) / 2);
            if (OccValue_(interval_m) > ranking) {
                interval_r = interval_m - 1;
            }
            else {
                interval_l = interval_m;
            }
        }

#if DIFF_TO_DO_BINARY_SEARCH > 0
        PrefectchOccValue_(interval_l);

        if (interval_r > interval_l) {
            while (OccValue_(interval_l + 1) <= ranking) {
                ++interval_l;
            }
        }

#endif

        int64_t pos = (int64_t)interval_l * kBitsPerInterval;
        unsigned long long *cur_word = packed_text_ + pos / kBitsPerWord;
        int pos_in_word = 0;
        __builtin_prefetch(cur_word);

        int remaining_ones = ranking + 1 - OccValue_(interval_l);
        int popcnt;

        for (; ; pos_in_word += kBitsPerWord) {
            popcnt = __builtin_popcountll(*cur_word);

            if (popcnt >= remaining_ones) {
                break;
            }
            else {
                remaining_ones -= popcnt;
            }

            ++cur_word;
        }

        return pos + pos_in_word + SelectInWord_(remaining_ones, *cur_word);
    }

    int64_t Pred(int64_t pos) {
        unsigned long long *word = packed_text_ + pos / kBitsPerWord;
        int idx_in_word = pos % kBitsPerWord;

        while (pos >= 0 && !((*word >> idx_in_word) & 1)) {
            --idx_in_word;
            --pos;

            if (idx_in_word < 0) {
                idx_in_word = kBitsPerWord - 1;
                --word;
            }
        }

        return pos;
    }

    int64_t Succ(int64_t pos) {
        unsigned long long *word = packed_text_ + pos / kBitsPerWord;
        int idx_in_word = pos % kBitsPerWord;

        while (pos < length && !((*word >> idx_in_word) & 1)) {
            ++idx_in_word;
            ++pos;

            if (idx_in_word == kBitsPerWord) {
                idx_in_word = 0;
                ++word;
            }
        }

        return pos;
    }

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