/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong
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

#ifndef SUCCINCT_DBG_H_
#define SUCCINCT_DBG_H_
#include <assert.h>
#include <vector>
#include "definitions.h"
#include "rank_and_select.h"
#include "khash.h"

using std::vector;
KHASH_MAP_INIT_INT64(k64v16, multi_t); // declare khash

class SuccinctDBG {
  public:
    // constants
    static const int kAlphabetSize = 4;
    static const int kWAlphabetSize = 9;
    static const int kWBitsPerChar = 4;
    static const int kWCharsPerWord = sizeof(unsigned long long) * kBitsPerByte / kWBitsPerChar;
    static const int kWCharMask = 0xF;
    static const int kMaxKmerK = 128;
    static const int kCharsPerUint32 = 16;
    static const int kBitsPerChar = 2;

    int64_t size;
    int kmer_k;

  public:
    SuccinctDBG(): need_to_free_(false) { }
    ~SuccinctDBG() {
        if (need_to_free_) {
            free(last_);
            free(w_);
            free(invalid_);
            free(is_dollar_);
            free(dollar_node_seq_);
        }

        if (need_to_free_mul_) {
            free(edge_multiplicities_);
            kh_destroy(k64v16, large_multi_h_);
        }
    }

    void LoadFromFile(const char *dbg_name);
    void init(unsigned long long *w, unsigned long long *last, long long *f, int64_t size, int kmer_k) {
        w_ = w;
        last_ = last;
        this->size = size;
        this->kmer_k = kmer_k;
        for (int i = 0; i < kAlphabetSize + 2; ++i) {
            f_[i] = f[i];
        }
        rs_w_.Build(w_, size);
        rs_last_.Build(last_, size);
    }

    uint8_t GetW(int64_t x) {
        return (*(w_ + x / kWCharsPerWord) >> (x % kWCharsPerWord * kWBitsPerChar)) & kWCharMask;
    }

    bool IsLast(int64_t x) {
        return (*(last_ + x / 64) >> (x % 64)) & 1;
    }

    bool IsLastOrDollar(int64_t x) {
        return ((last_[x / 64] | is_dollar_[x / 64]) >> (x % 64)) & 1;
    }

    int64_t GetLastIndex(int64_t x) {
        return rs_last_.Succ(x);
    }

    uint8_t GetNodeLastChar(int64_t x) {
        for (int i = 1; ; ++i) {
            if (f_[i] > x) {
                return i - 1;
            }
        }
    }

    bool IsValidNode(int64_t x) {
        return !((*(invalid_ + x / 64) >> (x % 64)) & 1);
    }

    bool IsDollarNode(int64_t x) {
        return (*(is_dollar_ + x / 64) >> (x % 64)) & 1;
    }

    void SetValid(int64_t x) {
        x = rs_last_.Succ(x);
        do {
            __sync_fetch_and_and(invalid_ + x / 64, ~(1ULL << (x % 64)));
            --x;
        } while (x >= 0 && IsLastOrDollar(x) == 0);
    }

    void SetInvalid(int64_t x) {
        x = rs_last_.Succ(x);
        do {
            __sync_fetch_and_or(invalid_ + x / 64, 1ULL << (x % 64));
            --x;
        } while (x >= 0 && IsLastOrDollar(x) == 0);
    }

    int EdgeMultiplicity(int64_t x) {
        if (__builtin_expect(edge_multiplicities_[x] != kMulti2Sp, 1)) {
            return edge_multiplicities_[x];
        } else {
            return kh_value(large_multi_h_, kh_get(k64v16, large_multi_h_, x));
        }
    }

    int NodeMultiplicity(int64_t x);

    int64_t Forward(int64_t x) { // the last node x points to
        uint8_t a = GetW(x);
        if (a > 4) {
            a -= 4;
        }
        int64_t count_a = rs_w_.Rank(a, x);
        return rs_last_.Select(rs_last_.Rank(f_[a] - 1) + count_a - 1);
    }

    int64_t Backward(int64_t x) { // the first node points to x
        uint8_t a = GetNodeLastChar(x);
        int64_t count_a = rs_last_.Rank(x - 1) - rs_last_.Rank(f_[a] - 1);
        return rs_w_.Select(a, count_a);
    }

    int Indegree(int64_t x);
    int Outdegree(int64_t x);
    int Incomings(int64_t x, int64_t *incomings);
    int Outgoings(int64_t x, int64_t *outgoings);
    int Incomings(int64_t x, int64_t *incomings, int *edge_countings);
    int Outgoings(int64_t x, int64_t *outgoings, int *edge_countings);
    bool IndegreeZero(int64_t x);
    bool OutdegreeZero(int64_t x);

    int64_t UniqueIncoming(int64_t x); // whether there is a unique valid node y s.t. y->x, return y if exist, else return -1
    int64_t UniqueOutgoing(int64_t x); // whether there is a unique valid node y s.t. x->y, return y if exist, else return -1

    int64_t Outgoing(uint8_t c, int64_t x); // not implement
    int64_t Incoming(uint8_t c, int64_t x); // not implement

    int64_t Index(uint8_t *seq);
    int64_t IndexBinarySearch(uint8_t *seq);
    int Label(int64_t x, uint8_t *seq);
    int64_t ReverseComplement(int64_t x);

    // WARNING: use this with cautions
    // After that NodeMultiplicity() and EdgeMultiplicty() are invalid
    void FreeMul() {
        if (need_to_free_mul_) {
            free(edge_multiplicities_);
            kh_destroy(k64v16, large_multi_h_);
        }
        need_to_free_mul_ = false;
    }

  private:
    // main memory
    unsigned long long *w_;
    unsigned long long *last_;
    unsigned long long *is_dollar_;
    unsigned long long *invalid_;
    uint32_t *dollar_node_seq_;
    multi2_t *edge_multiplicities_;
    khash_t(k64v16) *large_multi_h_;

    long long f_[kAlphabetSize + 2];

    unsigned int num_dollar_nodes_;
    int uint32_per_dollar_nodes_;

    // auxiliary memory
    RankAndSelect4Bits rs_w_;
    RankAndSelect1Bit<false> rs_last_;
    RankAndSelect1Bit<true> rs_is_dollar_;
    bool need_to_free_;
    bool need_to_free_mul_;

    void PrefixRangeSearch_(uint8_t c, int64_t &l, int64_t &r);
};

#endif // SUCCINCT_DBG_H_