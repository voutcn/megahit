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

#include "succinct_dbg.h"

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <string>
#include <algorithm>

#include "sdbg_multi_io.h"
#include "mem_file_checker-inl.h"
#include "utils.h"

int SuccinctDBG::EdgeOutdegree(int64_t edge_id) {
    if (!IsValidEdge(edge_id)) {
        return -1;
    }

    int64_t outdegree = 0;
    int64_t next_edge = Forward(edge_id);

    do {
        if (IsValidEdge(next_edge)) {
            ++outdegree;
        }

        --next_edge;
    }
    while (next_edge >= 0 && !IsLastOrTip(next_edge));

    return outdegree;
}

int SuccinctDBG::EdgeIndegree(int64_t edge_id) {
    if (!IsValidEdge(edge_id)) {
        return -1;
    }

    int64_t first_income = Backward(edge_id);
    int8_t c = GetW(first_income);
    int count_ones = IsLastOrTip(first_income);
    int64_t indegree = IsValidEdge(first_income);

    for (int64_t y = first_income + 1; count_ones < 5 && y < this->size; ++y) {
        count_ones += IsLastOrTip(y);
        uint8_t cur_char = GetW(y);

        if (cur_char == c) {
            break;
        }
        else if (cur_char == c + 4 && IsValidEdge(y)) {
            ++indegree;
        }
    }

    return indegree;
}

int SuccinctDBG::OutgoingEdges(int64_t edge_id, int64_t *outgoings) {
    if (!IsValidEdge(edge_id)) {
        return -1;
    }

    int64_t outdegree = 0;
    int64_t next_edge = Forward(edge_id);

    do {
        if (IsValidEdge(next_edge)) {
            outgoings[outdegree] = next_edge;
            ++outdegree;
        }

        --next_edge;
    }
    while (next_edge >= 0 && !IsLastOrTip(next_edge));

    return outdegree;
}

int SuccinctDBG::IncomingEdges(int64_t edge_id, int64_t *incomings) {
    if (!IsValidEdge(edge_id)) {
        return -1;
    }

    int64_t first_income = Backward(edge_id);
    int8_t c = GetW(first_income);
    int count_ones = IsLastOrTip(first_income);
    int64_t indegree = IsValidEdge(first_income);

    if (indegree > 0) {
        incomings[0] = first_income;
    }

    for (int64_t y = first_income + 1; count_ones < 5 && y < this->size; ++y) {
        count_ones += IsLastOrTip(y);
        uint8_t cur_char = GetW(y);

        if (cur_char == c) {
            break;
        }
        else if (cur_char == c + 4 && IsValidEdge(y)) {
            incomings[indegree] = y;
            ++indegree;
        }
    }

    return indegree;
}


bool SuccinctDBG::EdgeIndegreeZero(int64_t edge_id) {
    if (!IsValidEdge(edge_id)) {
        return false;
    }

    int64_t first_income = Backward(edge_id);

    if (IsValidEdge(first_income)) {
        return false;
    }

    int8_t c = GetW(first_income);
    int count_ones = IsLastOrTip(first_income);

    for (int64_t y = first_income + 1; count_ones < 5 && y < this->size; ++y) {
        uint8_t cur_char = GetW(y);

        if (cur_char == c) {
            break;
        }
        else if (cur_char == c + 4 && IsValidEdge(y)) {
            return false;
        }

        count_ones += IsLastOrTip(y);
    }

    return true;
}

bool SuccinctDBG::EdgeOutdegreeZero(int64_t edge_id) {
    if (!IsValidEdge(edge_id)) {
        return false;
    }

    int64_t next_edge = Forward(edge_id);

    do {
        if (IsValidEdge(next_edge)) {
            return false;
        }

        --next_edge;
    }
    while (next_edge >= 0 && !IsLastOrTip(next_edge));

    return true;
}

int64_t SuccinctDBG::UniqueNextEdge(int64_t edge_id) {
    if (!IsValidEdge(edge_id)) {
        return -1;
    }

    int64_t next_edge = Forward(edge_id);
    int64_t ret = -1;

    do {
        if (IsValidEdge(next_edge)) {
            if (ret != -1) {
                return -1;
            }

            ret = next_edge;
        }

        --next_edge;
    }
    while (next_edge >= 0 && !IsLastOrTip(next_edge));

    return ret;
}

int64_t SuccinctDBG::UniquePrevEdge(int64_t edge_id) {
    if (!IsValidEdge(edge_id)) {
        return -1;
    }

    int64_t prev_edge = Backward(edge_id);
    uint8_t c = GetW(prev_edge);
    int64_t ret = IsValidEdge(prev_edge) ? prev_edge : -1;
    int count_ones = IsLastOrTip(prev_edge);

    for (++prev_edge; count_ones < 5 && prev_edge < this->size; ++prev_edge) {
        count_ones += IsLastOrTip(prev_edge);
        uint8_t cur_char = GetW(prev_edge);

        if (cur_char == c) {
            break;
        }
        else if (cur_char == c + 4 && IsValidEdge(prev_edge)) {
            if (ret != -1) {
                return -1;
            }
            else {
                ret = prev_edge;
            }
        }
    }

    return ret;
}

int64_t SuccinctDBG::PrevSimplePathEdge(int64_t edge_id) {
    int64_t prev_edge = UniquePrevEdge(edge_id);

    if (prev_edge != -1 && UniqueNextEdge(prev_edge) != -1) {
        return prev_edge;
    }
    else {
        return -1;
    }
}

int64_t SuccinctDBG::NextSimplePathEdge(int64_t edge_id) {
    int64_t next_edge = UniqueNextEdge(edge_id);

    if (next_edge != -1 && UniquePrevEdge(next_edge) != -1) {
        return next_edge;
    }
    else {
        return -1;
    }
}

bool SuccinctDBG::NodeOutdegreeZero(int64_t node_id) {
    int64_t edge_id = GetLastIndex(node_id);

    do {
        if (IsValidEdge(edge_id)) {
            return false;
        }

        --edge_id;
    }
    while (edge_id >= 0 && !IsLastOrTip(edge_id));

    return true;
}

bool SuccinctDBG::NodeIndegreeZero(int64_t node_id) {
    int64_t prev_edge = Backward(node_id);

    if (IsValidEdge(prev_edge)) {
        return false;
    }

    uint8_t c = GetW(prev_edge);
    int count_ones = IsLastOrTip(prev_edge);

    for (++prev_edge; count_ones < 5 && prev_edge < this->size; ++prev_edge) {
        count_ones += IsLastOrTip(prev_edge);
        uint8_t cur_char = GetW(prev_edge);

        if (cur_char == c) {
            break;
        }
        else if (cur_char == c + 4 && IsValidEdge(prev_edge)) {
            return false;
        }
    }

    return true;
}

int64_t SuccinctDBG::UniquePrevNode(int64_t node_id) {
    int64_t prev_edge = Backward(node_id);
    int64_t ret = IsValidEdge(prev_edge) ? prev_edge : -1;
    uint8_t c = GetW(prev_edge);
    int count_ones = IsLastOrTip(prev_edge);

    for (++prev_edge; count_ones < 5 && prev_edge < this->size; ++prev_edge) {
        count_ones += IsLastOrTip(prev_edge);
        uint8_t cur_char = GetW(prev_edge);

        if (cur_char == c) {
            break;
        }
        else if (cur_char == c + 4 && IsValidEdge(prev_edge)) {
            if (ret != -1) {
                return -1;
            }
            else {
                ret = prev_edge;
            }
        }
    }

    return ret == -1 ? -1 : GetLastIndex(ret);
}

int64_t SuccinctDBG::UniqueNextNode(int64_t node_id) {
    int64_t edge_id = GetLastIndex(node_id);
    int64_t ret = -1;

    do {
        if (IsValidEdge(edge_id)) {
            if (ret != -1) {
                return -1;
            }
            else {
                ret = Forward(edge_id);
            }
        }

        --edge_id;
    }
    while (edge_id >= 0 && !IsLastOrTip(edge_id));

    return ret == -1 ? -1 : GetLastIndex(ret);
}

void SuccinctDBG::DeleteAllEdges(int64_t node_id) {
    int64_t edge_id = GetLastIndex(node_id);

    do {
        SetInvalidEdge(edge_id);
        --edge_id;
    }
    while (edge_id >= 0 && !IsLastOrTip(edge_id));

    edge_id = Backward(node_id);
    uint8_t c = GetW(edge_id);
    int count_ones = IsLastOrTip(edge_id);
    SetInvalidEdge(edge_id);

    for (++edge_id; count_ones < 5 && edge_id < this->size; ++edge_id) {
        count_ones += IsLastOrTip(edge_id);
        uint8_t cur_char = GetW(edge_id);

        if (cur_char == c) {
            break;
        }
        else if (cur_char == c + 4) {
            SetInvalidEdge(edge_id);
        }
    }
}

int SuccinctDBG::NextNodes(int64_t node_id, int64_t next[]) {
    int64_t next_edge = GetLastIndex(node_id);
    int outd = 0;

    do {
        if (IsValidEdge(next_edge)) {
            next[outd++] = GetLastIndex(Forward(next_edge));
        }

        --next_edge;
    }
    while (next_edge >= 0 && !IsLastOrTip(next_edge));

    return outd;
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
    // int64_t l = f_[seq[kmer_k - 1]];
    // int64_t r = f_[seq[kmer_k - 1] + 1] - 1;

    // only work if k > 8
    int pre = 0;
    for (int i = 0; i < prefix_lk_len_; ++i) {
        pre = pre * 4 + seq[kmer_k - 1 - i] - 1;
    }
    int64_t l = prefix_lkt_[pre * 2];
    int64_t r = prefix_lkt_[pre * 2 + 1];

    while (l <= r) {
        int cmp = 0;
        int64_t mid = (l + r) / 2;
        int64_t y = mid;

        for (int i = kmer_k - 1; i >= 0; --i) {
            if (IsTip(y)) {
                uint32_t *tip_node_seq = tip_node_seq_ + (size_t)uint32_per_tip_nodes_ * (rs_is_tip_.Rank(y) - 1);

                for (int j = 0; j < i; ++j) {
                    uint8_t c = (tip_node_seq[j / kCharsPerUint32] >> (kCharsPerUint32 - 1 - j % kCharsPerUint32) * kBitsPerChar) & 3;
                    c++;

                    if (c < seq[i - j]) {
                        cmp = -1;
                        break;
                    }
                    else if (c > seq[i - j]) {
                        cmp = 1;
                        break;
                    }
                }

                if (cmp == 0) {
                    if(IsTip(mid)) {
                        cmp = -1;
                    }
                    else {
                        uint8_t c = (tip_node_seq[i / kCharsPerUint32] >> (kCharsPerUint32 - 1 - i % kCharsPerUint32) * kBitsPerChar) & 3;
                        c++;

                        if (c < seq[0]) {
                            cmp = -1;
                            break;
                        }
                        else if (c > seq[0]) {
                            cmp = 1;
                            break;
                        }
                    }
                }

                break;
            }

            y = Backward(y);
            uint8_t c = GetW(y);

            if (c < seq[i]) {   
                cmp = -1;
                break;
            }
            else if (c > seq[i]) {
                cmp = 1;
                break;
            }
        }

        if (cmp == 0) {
            return GetLastIndex(mid);
        }
        else if (cmp > 0) {
            r = mid - 1;
        }
        else {
            l = mid + 1;
        }
    }

    return -1;
}

int64_t SuccinctDBG::IndexBinarySearchEdge(uint8_t *seq) {
    int64_t node = IndexBinarySearch(seq);
    if (node == -1) { return -1; }
    do {
        uint8_t edge_label = GetW(node);
        if (edge_label == seq[kmer_k] || edge_label - 4 == seq[kmer_k]) {
            return node;
        }

        --node;
    }
    while (node >= 0 && !IsLastOrTip(node));
    return -1;
}


int SuccinctDBG::Label(int64_t edge_or_node_id, uint8_t *seq) {
    int64_t x = edge_or_node_id;

    for (int i = kmer_k - 1; i >= 0; --i) {
        if (IsTip(x)) {
            uint32_t *tip_node_seq = tip_node_seq_ + (size_t)uint32_per_tip_nodes_ * (rs_is_tip_.Rank(x) - 1);

            for (int j = 0; j <= i; ++j) {
                seq[i - j] = (tip_node_seq[j / kCharsPerUint32] >> (kCharsPerUint32 - 1 - j % kCharsPerUint32) * kBitsPerChar) & 3;
                seq[i - j]++;
            }

            break;
        }

        x = Backward(x);
        seq[i] = GetW(x);
        assert(seq[i] > 0);

        if (seq[i] > 4) {
            seq[i] -= 4;
        }
    }

    return kmer_k;
}

int64_t SuccinctDBG::EdgeReverseComplement(int64_t edge_id) {
    if (!IsValidEdge(edge_id)) {
        return -1;
    }

    uint8_t seq[kMaxKmerK + 1];
    assert(kmer_k == Label(edge_id, seq));
    seq[kmer_k] = GetW(edge_id);

    if (seq[kmer_k] > 4) {
        seq[kmer_k] -= 4;
    }

    int i, j;

    for (i = 0, j = kmer_k; i < j; ++i, --j) {
        std::swap(seq[i], seq[j]);
        seq[i] = 5 - seq[i];
        seq[j] = 5 - seq[j];
    }

    if (i == j) {
        seq[i] = 5 - seq[i];
    }

    int64_t rev_node = IndexBinarySearch(seq);

    if (rev_node == -1) return - 1;

    do {
        uint8_t edge_label = GetW(rev_node);

        if (edge_label == seq[kmer_k] || edge_label - 4 == seq[kmer_k]) {
            return rev_node;
        }

        --rev_node;
    }
    while (rev_node >= 0 && !IsLastOrTip(rev_node));

    return -1;
}

void SuccinctDBG::LoadFromMultiFile(const char *dbg_name, bool need_multiplicity) {
    SdbgReader sdbg_reader;
    sdbg_reader.set_file_prefix(std::string(dbg_name));
    sdbg_reader.read_info();
    sdbg_reader.init_files();

    for (int i = 0; i < 6; ++i) {
        f_[i] = sdbg_reader.f()[i];
    }

    kmer_k = sdbg_reader.kmer_size();
    size = sdbg_reader.num_items();
    num_tip_nodes_ = sdbg_reader.num_tips();
    uint32_per_tip_nodes_ = sdbg_reader.words_per_tip_label();
    prefix_lk_len_ = sdbg_reader.prefix_lkt_len();

    size_t word_needed_w = (size + kWCharsPerWord - 1) / kWCharsPerWord;
    size_t word_needed_last = (size + kBitsPerULL - 1) / kBitsPerULL;

    w_ = (unsigned long long *) MallocAndCheck(sizeof(unsigned long long) * word_needed_w, __FILE__, __LINE__);
    last_ = (unsigned long long *) MallocAndCheck(sizeof(unsigned long long) * word_needed_last, __FILE__, __LINE__);
    is_tip_ = (unsigned long long *) MallocAndCheck(sizeof(unsigned long long) * word_needed_last, __FILE__, __LINE__);
    tip_node_seq_ = (uint32_t *) MallocAndCheck(sizeof(uint32_t) * num_tip_nodes_ * sdbg_reader.words_per_tip_label(), __FILE__, __LINE__);
    prefix_lkt_ = (int64_t *) MallocAndCheck(sizeof(int64_t) * sdbg_reader.prefix_lkt_size() * 2, __FILE__, __LINE__);

    if (need_multiplicity) {
        if (sdbg_reader.num_large_mul() > (1 << 30) ||
            sdbg_reader.num_large_mul() > size * 0.08) {
            edge_large_multi_ = (multi_t *) MallocAndCheck(sizeof(multi_t) * size, __FILE__, __LINE__);
        } else {
            edge_multi_ = (multi2_t *) MallocAndCheck(sizeof(multi2_t) * size, __FILE__, __LINE__);
            large_multi_h_ = kh_init(k64v16);
        }
        need_to_free_mul_ = true;
    }
    else {
        is_multi_1_ = (unsigned long long *) MallocAndCheck(sizeof(unsigned long long) * word_needed_last, __FILE__, __LINE__);
        memset(is_multi_1_, 0, sizeof(unsigned long long) * word_needed_last);
        need_to_free_mul_ = true;
    }

    for (int i = 0; i < sdbg_reader.prefix_lkt_size() * 2; ++i) {
        prefix_lkt_[i] = sdbg_reader.prefix_lkt(i);
    }

    unsigned long long packed_w = 0;
    unsigned long long packed_last = 0;
    unsigned long long packed_tip = 0;

    int64_t w_word_idx = 0;
    int64_t last_word_idx = 0;

    int w_word_offset = 0;
    int last_word_offset = 0;

    int64_t tip_label_offset = 0;
    uint16_t item = 0;

    for (long long i = 0; LIKELY(i < size); ++i) {
        assert(sdbg_reader.NextItem(item));

        packed_w |= (unsigned long long)(item & 0xF) << w_word_offset;
        w_word_offset += kWBitsPerChar;

        if (w_word_offset == kBitsPerULL) {
            w_[w_word_idx++] = packed_w;
            w_word_offset = 0;
            packed_w = 0;
        }

        packed_last |= (unsigned long long)((item >> 4) & 1) << last_word_offset;
        packed_tip |= (unsigned long long)((item >> 5) & 1) << last_word_offset;

        last_word_offset++;

        if (last_word_offset == kBitsPerULL) {
            last_[last_word_idx] = packed_last;
            is_tip_[last_word_idx++] = packed_tip;
            last_word_offset = 0;
            packed_tip = packed_last = 0;
        }

        if (need_multiplicity) {
            if (edge_multi_)
                edge_multi_[i] = item >> 8;
            else
                edge_large_multi_[i] = item >> 8;
        } else {
            is_multi_1_[i / 64] |= (unsigned long long)((item >> 8) <= 1) << (i % 64);
        }

        if (UNLIKELY((item >> 8) == kMulti2Sp)) {
            multi_t mul = sdbg_reader.NextLargeMul();
            assert(mul >= kMulti2Sp);
            if (need_multiplicity) {
                if (edge_multi_) {
                    int ret;
                    khint_t k = kh_put(k64v16, large_multi_h_, i, &ret);
                    kh_value(large_multi_h_, k) = mul;
                } else {
                    edge_large_multi_[i] = mul;
                }
            }
        }

        if (UNLIKELY((item >> 5) & 1)) {
            sdbg_reader.NextTipLabel(tip_node_seq_ + tip_label_offset);
            tip_label_offset += sdbg_reader.words_per_tip_label();
        }
    }

    if (w_word_offset != 0) {
        w_[w_word_idx++] = packed_w;
    }

    if (last_word_offset != 0) {
        last_[last_word_idx] = packed_last;
        is_tip_[last_word_idx++] = packed_tip;
    }

    assert(!sdbg_reader.NextItem(item));
    assert(tip_label_offset == num_tip_nodes_ * sdbg_reader.words_per_tip_label());

    invalid_ = (unsigned long long *) MallocAndCheck(sizeof(unsigned long long) * word_needed_last, __FILE__, __LINE__);
    memcpy(invalid_, is_tip_, sizeof(unsigned long long) * word_needed_last);
    rs_is_tip_.Build(is_tip_, size);

    init(w_, last_, f_, size, kmer_k);
    need_to_free_ = true;
}

void SuccinctDBG::PrefixRangeSearch_(uint8_t c, int64_t &l, int64_t &r) {
    int64_t low = l - 1;
    unsigned long long *word_last = last_ + low / 64;
    unsigned long long *word_is_d = is_tip_ + low / 64;
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
    word_is_d = is_tip_ + high / 64;
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
    }
    else {
        r = -1;
        return;
    }

    // the first c/c- in [low, high]
    c_pos = std::min(rs_w_.Succ(c + 4, low), rs_w_.Succ(c, low));

    if (c_pos <= high) {
        l = Forward(c_pos);
    }
    else {
        l = -1;
        return;
    }
}