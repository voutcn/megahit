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

#include "rank_and_select.h"

#define DIFF_TO_DO_BINARY_SEARCH 2

RankAndSelect4Bits::RankAndSelect4Bits() {
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

RankAndSelect4Bits::~RankAndSelect4Bits() {
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

void RankAndSelect4Bits::Build(unsigned long long *packed_text, int64_t length) {
    int64_t count[kAlphabetSize];
    int64_t num_intervals = (length + kCharPerInterval - 1) / kCharPerInterval + 1;
    int64_t num_intervals_major = (length + kCharPerIntervalMajor - 1) / kCharPerIntervalMajor + 1;
    unsigned long long *cur_word = packed_text;

    // build rank
    for (int i = 0; i < kAlphabetSize; ++i) {
        count[i] = 0;
        occ_value_explicit_major_[i] = (int64_t*) malloc(sizeof(int64_t) * num_intervals_major);
        if (occ_value_explicit_major_[i] == NULL) {
            fprintf(stderr, "Malloc Failed: %s: %d\n", __FILE__, __LINE__);
            exit(1);
        }
        occ_value_explicit_minor_[i] = (uint16_t*) malloc(sizeof(uint16_t) * num_intervals);
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
        rank_to_interval_explicit_[j] = (interval_t*) malloc(sizeof(interval_t) * s_table_size);
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

int64_t RankAndSelect4Bits::Rank(uint8_t c, int64_t pos) {
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

    } else if (sampled_index < pos) {
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

    } else {
        return OccValue_(c, which_interval);
    }
}

int64_t RankAndSelect4Bits::Select(uint8_t c, int64_t ranking) {
    if (ranking >= char_frequency[c]) {
        return length;
    } else if (ranking < 0) {
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
        } else {
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
        } else {
            remaining_c -= popcnt;
        }
        ++cur_word;
    }

    return pos + pos_in_word + SelectInWord_(c, remaining_c, *cur_word);
}

int64_t RankAndSelect4Bits::Pred(uint8_t c, int64_t pos) {
    if (((*(packed_text_ + pos / kCharPerWord) >> (pos % kCharPerWord * kBitsPerChar)) & ((1 << kBitsPerChar) - 1)) == c) {
        return pos;
    }
    return Select(c, Rank(c, pos) - 1);
}

int64_t RankAndSelect4Bits::PredLimitedStep(uint8_t c, int64_t pos, int step) {
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

int64_t RankAndSelect4Bits::Succ(uint8_t c, int64_t pos) {
    if (((*(packed_text_ + pos / kCharPerWord) >> (pos % kCharPerWord * kBitsPerChar)) & ((1 << kBitsPerChar) - 1)) == c) {
        return pos;
    }
    return Select(c, Rank(c, pos - 1));
}

int64_t RankAndSelect4Bits::SuccLimitedStep(uint8_t c, int64_t pos, int step) {
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

RankAndSelect1Bit::RankAndSelect1Bit() {
    occ_value_explicit_minor_ = NULL;
    occ_value_explicit_major_ = NULL;
    rank_to_interval_explicit_ = NULL;
}

RankAndSelect1Bit::~RankAndSelect1Bit() {
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

void RankAndSelect1Bit::Build(unsigned long long *packed_text, int64_t length) {
    int64_t count_ones = 0;
    int64_t num_intervals = (length + kBitsPerInterval - 1) / kBitsPerInterval + 1;
    int64_t num_intervals_major = (length + kBitsPerMajorInterval - 1) / kBitsPerMajorInterval + 1;
    unsigned long long *cur_word = packed_text;

    occ_value_explicit_major_ = (int64_t*) malloc(sizeof(int64_t) * num_intervals_major);
    if (occ_value_explicit_major_ == NULL) {
        fprintf(stderr, "Malloc Failed: %s: %d\n", __FILE__, __LINE__);
        exit(1);
    }
    occ_value_explicit_minor_ = (uint16_t*) malloc(sizeof(uint16_t) * num_intervals);
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
    uint32_t s_table_size = (count_ones + kSelectSampleSize - 1) / kSelectSampleSize + 1;
    rank_to_interval_explicit_ = (uint32_t*) malloc(sizeof(uint32_t) * s_table_size);
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
    packed_text_ = packed_text;
    this->length = length;
}

int64_t RankAndSelect1Bit::Rank(int64_t pos) {
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

    } else if (sampled_index < pos) {
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

    } else {
        return OccValue_(which_interval);
    }
}

int64_t RankAndSelect1Bit::Select(int64_t ranking) {
    if (ranking >= total_num_ones) {
        return length;
    } else if (ranking < 0) {
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
        } else {
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
        } else {
            remaining_ones -= popcnt;
        }
        ++cur_word;
    }
    return pos + pos_in_word + SelectInWord_(remaining_ones, *cur_word);
}

int64_t RankAndSelect1Bit::Pred(int64_t pos) {
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

int64_t RankAndSelect1Bit::Succ(int64_t pos) {
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