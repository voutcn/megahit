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
    static const int kMaxKmerK = kMaxK + 1;
    static const int kCharsPerUint32 = 16;
    static const int kBitsPerChar = 2;

    int64_t size;
    int kmer_k;

  public:
    SuccinctDBG(): need_to_free_(false), need_to_free_mul_(false), edge_multi_(NULL), edge_large_multi_(NULL), is_multi_1_(NULL) { }
    ~SuccinctDBG() {
        if (need_to_free_) {
            free(last_);
            free(w_);
            free(invalid_);
            free(is_tip_);
            free(tip_node_seq_);
            free(prefix_lkt_);
        }

        FreeMul();
    }

    void LoadFromMultiFile(const char *dbg_name, bool need_multiplicity = true);
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

        for (int i = 1; i < kAlphabetSize + 2; ++i) {
            rank_f_[i] = rs_last_.Rank(f_[i] - 1);
        }

        for (int64_t i = 0; i < size; ++i) {
            if (GetW(i) == 0) {
                SetInvalidEdge(i);
            }
        }
    }

    uint8_t GetW(int64_t x) {
        return (*(w_ + x / kWCharsPerWord) >> (x % kWCharsPerWord * kWBitsPerChar)) & kWCharMask;
    }

    bool IsLast(int64_t x) {
        return (*(last_ + x / 64) >> (x % 64)) & 1;
    }

    bool IsLastOrTip(int64_t x) {
        return ((last_[x / 64] | is_tip_[x / 64]) >> (x % 64)) & 1;
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

    bool IsValidEdge(int64_t edge_id) {
        return !((*(invalid_ + edge_id / 64) >> (edge_id % 64)) & 1);
    }

    bool IsTip(int64_t edge_id) {
        return (*(is_tip_ + edge_id / 64) >> (edge_id % 64)) & 1;
    }

    void SetValidEdge(int64_t edge_id) {
        __sync_fetch_and_and(invalid_ + edge_id / 64, ~(1ULL << (edge_id % 64)));
    }

    void SetInvalidEdge(int64_t edge_id) {
        __sync_fetch_and_or(invalid_ + edge_id / 64, 1ULL << (edge_id % 64));
    }

    int EdgeMultiplicity(int64_t edge_id) {
        // if (is_multi_1_) {
        //     return ((is_multi_1_[edge_id / 64] >> (edge_id % 64)) & 1) + 1;
        // }

        if (edge_large_multi_) {
            return edge_large_multi_[edge_id];
        }

        if (__builtin_expect(edge_multi_[edge_id] != kMulti2Sp, 1)) {
            return edge_multi_[edge_id];
        }
        else {
            return kh_value(large_multi_h_, kh_get(k64v16, large_multi_h_, edge_id));
        }
    }

    bool IsMulti1(int64_t edge_id) {
        if (is_multi_1_)
            return (is_multi_1_[edge_id / 64] >> (edge_id % 64)) & 1;
        return EdgeMultiplicity(edge_id) <= 1;
    }

    int64_t Forward(int64_t edge_id) { // the last edge edge_id points to
        uint8_t a = GetW(edge_id);

        if (a > 4) {
            a -= 4;
        }

        int64_t count_a = rs_w_.Rank(a, edge_id);
        return rs_last_.Select(rank_f_[a] + count_a - 1);
    }

    int64_t Backward(int64_t edge_id) { // the first edge points to edge_id
        uint8_t a = GetNodeLastChar(edge_id);
        int64_t count_a = rs_last_.Rank(edge_id - 1) - rank_f_[a];
        return rs_w_.Select(a, count_a);
    }

    int64_t Index(uint8_t *seq);
    int64_t IndexBinarySearch(uint8_t *seq);
    int64_t IndexBinarySearchEdge(uint8_t *seq);
    int Label(int64_t edge_or_node_id, uint8_t *seq);

    int EdgeIndegree(int64_t edge_id);
    int EdgeOutdegree(int64_t edge_id);
    int IncomingEdges(int64_t edge_id, int64_t *incomings);
    int OutgoingEdges(int64_t edge_id, int64_t *outgoings);
    bool EdgeIndegreeZero(int64_t edge_id);
    bool EdgeOutdegreeZero(int64_t edge_id);
    int64_t UniqueNextEdge(int64_t edge_id);
    int64_t UniquePrevEdge(int64_t edge_id);
    int64_t PrevSimplePathEdge(int64_t edge_id);
    int64_t NextSimplePathEdge(int64_t edge_id);
    int64_t EdgeReverseComplement(int64_t edge_id);

    bool NodeOutdegreeZero(int64_t node_id);
    bool NodeIndegreeZero(int64_t node_id);
    int64_t UniquePrevNode(int64_t node_id);
    int64_t UniqueNextNode(int64_t node_id);
    void DeleteAllEdges(int64_t node_id);

    int NextNodes(int64_t node_id, int64_t next[]); // return out-degree, ids stored in next[]

    // WARNING: use this with cautions
    // After that EdgeMultiplicty() are invalid
    void FreeMul() {
        if (need_to_free_mul_) {
            if (edge_multi_ != NULL) {
                free(edge_multi_);
                edge_multi_ = NULL;
                kh_destroy(k64v16, large_multi_h_);
            } else if (edge_large_multi_ != NULL) {
                free(edge_large_multi_);
                edge_large_multi_ = NULL;
            }

            if (is_multi_1_ != NULL) {
                free(is_multi_1_);
            }
        }

        need_to_free_mul_ = false;
    }

  private:
    bool need_to_free_;
    bool need_to_free_mul_;

    // main memory
    unsigned long long *w_;
    unsigned long long *last_;
    unsigned long long *is_tip_;
    unsigned long long *invalid_;
    uint32_t *tip_node_seq_;
    int64_t *prefix_lkt_;
    multi2_t *edge_multi_;
    multi_t *edge_large_multi_;
    khash_t(k64v16) *large_multi_h_;
    unsigned long long *is_multi_1_;

    long long f_[kAlphabetSize + 2];
    long long rank_f_[kAlphabetSize + 2]; // = rs_last_.Rank(f_[i] - 1)

    int64_t num_tip_nodes_;
    int uint32_per_tip_nodes_;
    int prefix_lk_len_;

    // auxiliary memory
    RankAndSelect4Bits rs_w_;
    RankAndSelect1Bit<false> rs_last_;
    RankAndSelect1Bit<true> rs_is_tip_;

    void PrefixRangeSearch_(uint8_t c, int64_t &l, int64_t &r);
};

#endif // SUCCINCT_DBG_H_