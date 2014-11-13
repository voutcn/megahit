/*
 *  succinct_dbg.cpp
 *  This file is a part of MEGAHIT
 *  
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

#include "succinct_dbg.h"

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <string>
#include <algorithm>

#include "mem_file_checker-inl.h"

int SuccinctDBG::NodeMultiplicity(int64_t x) {
    int outgoing_multi = 0;
    int incoming_multi = 0;

    int64_t first_income = Backward(x);
    int8_t c = GetW(first_income);
    int count_ones = IsLastOrDollar(first_income);
    if (!IsDollarNode(first_income)) {
        incoming_multi += EdgeMultiplicity(first_income);
    }

    for (int64_t y = first_income + 1; count_ones < 5 && y < this->size; ++y) {
        count_ones += IsLastOrDollar(y);
        uint8_t cur_char = GetW(y);
        if (cur_char == c) { break; }
        else if (cur_char == c + 4 && !IsDollarNode(y)) {
            incoming_multi += EdgeMultiplicity(y);
        }
    }

    x = rs_last_.Succ(x);
    do {
        if (GetW(x) != 0) {
            outgoing_multi += EdgeMultiplicity(x);
        }
        --x;
    } while (x >= 0 && !IsLastOrDollar(x));

    return std::max(outgoing_multi, incoming_multi);
}

int SuccinctDBG::Outdegree(int64_t x) {
    int64_t outdegree = 0;
    x = rs_last_.Succ(x);
    do {
        if (GetW(x) != 0 && IsValidNode(Forward(x))) {
            ++outdegree;
        }
        --x;
    } while (x >= 0 && !IsLastOrDollar(x));
    return outdegree;
}

int SuccinctDBG::Indegree(int64_t x) {
    int64_t first_income = Backward(x);
    int8_t c = GetW(first_income);
    int count_ones = IsLastOrDollar(first_income);
    int64_t indegree = IsValidNode(first_income);

    for (int64_t y = first_income + 1; count_ones < 5 && y < this->size; ++y) {
        count_ones += IsLastOrDollar(y);
        uint8_t cur_char = GetW(y);
        if (cur_char == c) { break; }
        else if (cur_char == c + 4 && IsValidNode(y)) {
            ++indegree;
        }
    }
    return indegree;
}

int SuccinctDBG::Outgoings(int64_t x, int64_t *outgoings) {
    int64_t outdegree = 0;
    x = rs_last_.Succ(x);
    do {
        if (GetW(x) != 0 && IsValidNode(Forward(x))) {
            outgoings[outdegree] = Forward(x);
            ++outdegree;
        }
        --x;
    } while (x >= 0 && !IsLastOrDollar(x));
    return outdegree;
}

int SuccinctDBG::Outgoings(int64_t x, int64_t *outgoings, int *edge_countings) {
    int64_t outdegree = 0;
    x = rs_last_.Succ(x);
    do {
        if (GetW(x) != 0 && IsValidNode(Forward(x))) {
            outgoings[outdegree] = Forward(x);
            edge_countings[outdegree] = EdgeMultiplicity(x);
            ++outdegree;
        }
        --x;
    } while (x >= 0 && !IsLastOrDollar(x));
    return outdegree;
}

int SuccinctDBG::Incomings(int64_t x, int64_t *incomings) {
    int64_t first_income = Backward(x);
    int8_t c = GetW(first_income);
    int count_ones = IsLastOrDollar(first_income);
    int64_t indegree = IsValidNode(first_income);
    if (indegree > 0) { incomings[0] = GetLastIndex(first_income); }

    for (int64_t y = first_income + 1; count_ones < 5 && y < this->size; ++y) {
        count_ones += IsLastOrDollar(y);
        uint8_t cur_char = GetW(y);
        if (cur_char == c) { break; }
        else if (cur_char == c + 4 && IsValidNode(y)) {
            incomings[indegree] = GetLastIndex(y);
            ++indegree;
        }
    }
    return indegree;
}

int SuccinctDBG::Incomings(int64_t x, int64_t *incomings, int *edge_countings) {
    int64_t first_income = Backward(x);
    int8_t c = GetW(first_income);
    int count_ones = IsLastOrDollar(first_income);
    int64_t indegree = IsValidNode(first_income);
    if (indegree > 0) { incomings[0] = GetLastIndex(first_income); }

    for (int64_t y = first_income + 1; count_ones < 5 && y < this->size; ++y) {
        count_ones += IsLastOrDollar(y);
        uint8_t cur_char = GetW(y);
        if (cur_char == c) { break; }
        else if (cur_char == c + 4 && IsValidNode(y)) {
            incomings[indegree] = GetLastIndex(y);
            edge_countings[indegree] = EdgeMultiplicity(y);
            ++indegree;
        }
    }
    return indegree;
}

bool SuccinctDBG::OutdegreeZero(int64_t x) {
    x = rs_last_.Succ(x);
    do {
        if (GetW(x) != 0 && IsValidNode(Forward(x))) {
            return false;
        }
        --x;
    } while (x >= 0 && !IsLastOrDollar(x));
    return true;
}

bool SuccinctDBG::IndegreeZero(int64_t x) {
    int64_t first_income = Backward(x);
    int8_t c = GetW(first_income);
    int count_ones = IsLastOrDollar(first_income);
    if (IsValidNode(first_income) != 0) { return false; }

    for (int64_t y = first_income + 1; count_ones < 5 && y < this->size; ++y) {
        count_ones += IsLastOrDollar(y);
        uint8_t cur_char = GetW(y);
        if (cur_char == c) { break; }
        else if (cur_char == c + 4 && IsValidNode(y)) {
            return false;
        }
    }
    return true;
}

int64_t SuccinctDBG::UniqueOutgoing(int64_t x) {
    x = rs_last_.Succ(x);
    int64_t outgoing = -1;
    do {
        if (GetW(x) != 0) {
            int64_t y = Forward(x);
            if (IsValidNode(y)) {
                if (outgoing != -1) {
                    return -1;
                } else {
                    outgoing = y;
                }
            }
        }
        --x;
    } while (x >= 0 && !IsLastOrDollar(x));
    return outgoing;
}

int64_t SuccinctDBG::UniqueIncoming(int64_t x) {
    int64_t y = Backward(x);
    int64_t incoming = IsValidNode(y) ? y : -1;
    uint8_t c = GetW(y);
    int count_ones = IsLastOrDollar(y);

    for (++y; count_ones < 5 && y < this->size; ++y) {
        count_ones += IsLastOrDollar(y);
        uint8_t cur_char = GetW(y);
        if (cur_char == c) { break; }
        else if (cur_char == c + 4 && IsValidNode(y)) {
            if (incoming != -1) {
                return -1;
            } else {
                incoming = y;
            }
        }
    }
    return incoming;
}

int64_t SuccinctDBG::Index(uint8_t *seq) {
    int64_t l = f_[seq[0]];
    int64_t r = f_[seq[0] + 1] - 1;

    for (int i = 1; i < kmer_k; ++i) {
        PrefixRangeSearch_(seq[i], l, r);
        if (l == -1 || r == -1) {
            return IndexBinarySearch(seq);
        }
    }
    assert(l == r);
    return r;
}

int64_t SuccinctDBG::IndexBinarySearch(uint8_t *seq) {
    int64_t l = f_[seq[kmer_k - 1]];
    int64_t r = f_[seq[kmer_k - 1] + 1] - 1;

    while (l <= r) {
        int cmp = 0;
        int64_t mid = (l + r) / 2;
        int64_t y = mid;
        for (int i = kmer_k - 1; i >= 0; --i) {
            if (IsDollarNode(y)) {
                uint32_t *dollar_node_seq = dollar_node_seq_ + (size_t)uint32_per_dollar_nodes_ * (rs_is_dollar_.Rank(y) - 1);
                for (int j = 0; j < i; ++j) {
                    uint8_t c = (dollar_node_seq[j / kCharsPerUint32] >> (kCharsPerUint32 - 1 - j % kCharsPerUint32) * kBitsPerChar) & 3;
                    c++;
                    if (c < seq[i - j]) { cmp = -1; break; }
                    else if (c > seq[i - j]) { cmp = 1; break; }
                }

                if (cmp == 0) {
                    if(IsDollarNode(mid)) { cmp = -1; }
                    else {
                        uint8_t c = (dollar_node_seq[i / kCharsPerUint32] >> (kCharsPerUint32 - 1 - i % kCharsPerUint32) * kBitsPerChar) & 3;
                        c++;
                        if (c < seq[0]) { cmp = -1; break; }
                        else if (c > seq[0]) { cmp = 1; break; }
                    }
                }
                break;
            }

            y = Backward(y);
            uint8_t c = GetW(y);
            if (c < seq[i]) { cmp = -1; break; }
            else if (c > seq[i]) { cmp = 1; break; }
        }

        if (cmp == 0) {
            return GetLastIndex(mid);
        } else if (cmp > 0) {
            r = mid - 1;
        } else {
            l = mid + 1;
        }
    }
    return -1;
}

int SuccinctDBG::Label(int64_t x, uint8_t *seq) {
    for (int i = kmer_k - 1; i >= 0; --i) {
        if (IsDollarNode(x)) {
            uint32_t *dollar_node_seq = dollar_node_seq_ + (size_t)uint32_per_dollar_nodes_ * (rs_is_dollar_.Rank(x) - 1);
            for (int j = 0; j <= i; ++j) {
                seq[i - j] = (dollar_node_seq[j / kCharsPerUint32] >> (kCharsPerUint32 - 1 - j % kCharsPerUint32) * kBitsPerChar) & 3;
                seq[i - j]++;
            }
            break;
        }
        x = Backward(x);
        seq[i] = GetW(x);
        assert(seq[i] > 0);
        if (seq[i] > 4) { seq[i] -= 4; }
    }
    return kmer_k;
}

int64_t SuccinctDBG::ReverseComplement(int64_t x) {
    if (!IsValidNode(x)) {
        return -1;
    }
    uint8_t seq[kMaxKmerK];
    assert(kmer_k == Label(x, seq));

    int i, j;
    for (i = 0, j = kmer_k - 1; i < j; ++i, --j) {
        std::swap(seq[i], seq[j]);
        seq[i] = 5 - seq[i];
        seq[j] = 5 - seq[j];
    }
    if (i == j) {
        seq[i] = 5 - seq[i];
    }
    return IndexBinarySearch(seq);
}

void SuccinctDBG::LoadFromFile(const char *dbg_name) {
    FILE *w_file = OpenFileAndCheck((std::string(dbg_name) + ".w").c_str(), "rb");
    FILE *last_file = OpenFileAndCheck((std::string(dbg_name) + ".last").c_str(), "rb");
    FILE *f_file = OpenFileAndCheck((std::string(dbg_name) + ".f").c_str(), "r");
    FILE *is_dollar_file = OpenFileAndCheck((std::string(dbg_name) + ".isd").c_str(), "rb");
    FILE *dollar_node_seq_file = OpenFileAndCheck((std::string(dbg_name) + ".dn").c_str(), "rb");
    FILE *edge_multiplicity_file = OpenFileAndCheck((std::string(dbg_name) + ".mul").c_str(), "rb");
    assert(w_file != NULL);
    assert(last_file != NULL);
    assert(f_file != NULL);
    assert(dollar_node_seq_file != NULL);

    for (int i = 0; i < kAlphabetSize + 2; ++i) {
        assert(fscanf(f_file, "%lld", &f_[i]) == 1);
    }
    fscanf(f_file, "%d", &kmer_k);
    fscanf(f_file, "%u", &num_dollar_nodes_);
    size = f_[kAlphabetSize + 1];

    size_t word_needed_w = (size + kWCharsPerWord - 1) / kWCharsPerWord;
    size_t word_needed_last = (size + kBitsPerULL - 1) / kBitsPerULL;
    w_ = (unsigned long long*) MallocAndCheck(sizeof(unsigned long long) * word_needed_w, __FILE__, __LINE__);
    last_ = (unsigned long long*) MallocAndCheck(sizeof(unsigned long long) * word_needed_last, __FILE__, __LINE__);
    is_dollar_ = (unsigned long long*) MallocAndCheck(sizeof(unsigned long long) * word_needed_last, __FILE__, __LINE__);
    invalid_ = (unsigned long long*) MallocAndCheck(sizeof(unsigned long long) * word_needed_last, __FILE__, __LINE__);
    edge_multiplicities_ = (multi_t*) MallocAndCheck(sizeof(multi_t) * size, __FILE__, __LINE__);

    assert(w_ != NULL);
    assert(last_ != NULL);
    assert(is_dollar_ != NULL);
    assert(invalid_ != NULL);

    size_t word_read = fread(w_, sizeof(unsigned long long), word_needed_w, w_file);
    assert(word_read == word_needed_w);
    word_read = fread(last_, sizeof(unsigned long long), word_needed_last, last_file);
    assert(word_read == word_needed_last);
    word_read = fread(is_dollar_, sizeof(unsigned long long), word_needed_last, is_dollar_file);
    assert(word_read == word_needed_last);
    memcpy(invalid_, is_dollar_, sizeof(unsigned long long) * word_needed_last);
    if (edge_multiplicity_file != NULL) {
        word_read = fread(edge_multiplicities_, sizeof(multi_t), size, edge_multiplicity_file);
        assert(word_read == (size_t)size);
    } else {
        memset(edge_multiplicities_, 0xFF, sizeof(multi_t) * size);
    }

    // read dollar nodes sequences
    assert(fread(&uint32_per_dollar_nodes_, sizeof(uint32_t), 1, dollar_node_seq_file) == 1);
    dollar_node_seq_ = (uint32_t*) MallocAndCheck(sizeof(uint32_t) * num_dollar_nodes_ * uint32_per_dollar_nodes_, __FILE__, __LINE__);
    assert(dollar_node_seq_ != NULL);

    word_read = fread(dollar_node_seq_, sizeof(uint32_t), (size_t)num_dollar_nodes_ * uint32_per_dollar_nodes_, dollar_node_seq_file);
    assert(word_read == (size_t)num_dollar_nodes_ * uint32_per_dollar_nodes_);

    rs_is_dollar_.Build(is_dollar_, size);

    init(w_, last_, f_, size, kmer_k);
    need_to_free_ = true;

    fclose(w_file);
    fclose(last_file);
    fclose(f_file); 
    fclose(is_dollar_file);
    fclose(dollar_node_seq_file);
    if (edge_multiplicity_file != NULL) {
        fclose(edge_multiplicity_file);
    }
}

void SuccinctDBG::PrefixRangeSearch_(uint8_t c, int64_t &l, int64_t &r) {
    int64_t low = l - 1;
    unsigned long long *word_last = last_ + low / 64;
    unsigned long long *word_is_d = is_dollar_ + low / 64;
    unsigned long long word = *word_last | *word_is_d;

    int idx_in_word = low % 64;
    while (low >= 0 && !((word >> idx_in_word) & 1)) {
        --idx_in_word;
        --low;
        if (idx_in_word < 0) {
            idx_in_word = 64 - 1;
            --word_last;
            --word_is_d;
            word = *word_last | *word_is_d;
        }
    }
    ++low;

    int64_t high = r;
    word_last = last_ + high / 64;
    word_is_d = is_dollar_ + high / 64;
    word = *word_last | *word_is_d;

    idx_in_word = high % 64;
    while (high < size && !((word >> idx_in_word) & 1)) {
        ++idx_in_word;
        ++high;
        if (idx_in_word == 64) {
            idx_in_word = 0;
            ++word_last;
            ++word_is_d;
            word = *word_last | *word_is_d;
        }
    }

    // the last c/c- in [low, high]
    int64_t c_pos = std::max(rs_w_.Pred(c + 4, high), rs_w_.Pred(c, high));
    if (c_pos >= low) {
        r = Forward(c_pos);
    } else {
        r = -1;
        return;
    }
    // the first c/c- in [low, high]
    c_pos = std::min(rs_w_.Succ(c + 4, low), rs_w_.Succ(c, low));
    if (c_pos <= high) {
        l = Forward(c_pos);
    } else {
        l = -1;
        return;
    }
}