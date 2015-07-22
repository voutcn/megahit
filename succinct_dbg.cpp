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
    if (!IsValidEdge(edge_id)) { return -1; }

    int64_t outdegree = 0;
    int64_t next_edge = Forward(edge_id);
    do {
        if (IsValidEdge(next_edge)) {
            ++outdegree;
        }
        --next_edge;
    } while (next_edge >= 0 && !IsLastOrTip(next_edge));
    return outdegree;
}

int SuccinctDBG::EdgeIndegree(int64_t edge_id) {
    if (!IsValidEdge(edge_id)) { return -1; }

    int64_t first_income = Backward(edge_id);
    int8_t c = GetW(first_income);
    int count_ones = IsLastOrTip(first_income);
    int64_t indegree = IsValidEdge(first_income);

    for (int64_t y = first_income + 1; count_ones < 5 && y < this->size; ++y) {
        count_ones += IsLastOrTip(y);
        uint8_t cur_char = GetW(y);
        if (cur_char == c) {
            break;
        } else if (cur_char == c + 4 && IsValidEdge(y)) {
            ++indegree;
        }
    }
    return indegree;
}

int SuccinctDBG::OutgoingEdges(int64_t edge_id, int64_t *outgoings) {
    if (!IsValidEdge(edge_id)) { return -1; }

    int64_t outdegree = 0;
    int64_t next_edge = Forward(edge_id);

    do {
        if (IsValidEdge(next_edge)) {
            outgoings[outdegree] = next_edge;
            ++outdegree;
        }
        --next_edge;
    } while (next_edge >= 0 && !IsLastOrTip(next_edge));
    return outdegree;
}

int SuccinctDBG::IncomingEdges(int64_t edge_id, int64_t *incomings) {
    if (!IsValidEdge(edge_id)) { return -1; }

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
        } else if (cur_char == c + 4 && IsValidEdge(y)) {
            incomings[indegree] = y;
            ++indegree;
        }
    }
    return indegree;
}


bool SuccinctDBG::EdgeIndegreeZero(int64_t edge_id) {
    if (!IsValidEdge(edge_id)) { return false; }

    int64_t first_income = Backward(edge_id);
    if (IsValidEdge(first_income)) { return false; }
    int8_t c = GetW(first_income);
    int count_ones = IsLastOrTip(first_income);

    for (int64_t y = first_income + 1; count_ones < 5 && y < this->size; ++y) {
        uint8_t cur_char = GetW(y);
        if (cur_char == c) {
            break;
        } else if (cur_char == c + 4 && IsValidEdge(y)) {
            return false;
        }
        count_ones += IsLastOrTip(y);
    }
    return true;
}

bool SuccinctDBG::EdgeOutdegreeZero(int64_t edge_id) {
    if (!IsValidEdge(edge_id)) { return false; }

    int64_t next_edge = Forward(edge_id);

    do {
        if (IsValidEdge(next_edge)) {
            return false;
        }
        --next_edge;
    } while (next_edge >= 0 && !IsLastOrTip(next_edge));

    return true;
}

int64_t SuccinctDBG::UniqueNextEdge(int64_t edge_id) {
    if (!IsValidEdge(edge_id)) { return -1; }

    int64_t next_edge = Forward(edge_id);
    int64_t ret = -1;

    do {
        if (IsValidEdge(next_edge)) {
            if (ret != -1) { return -1; }
            ret = next_edge;
        }
        --next_edge;
    } while (next_edge >= 0 && !IsLastOrTip(next_edge));

    return ret;
}

int64_t SuccinctDBG::UniquePrevEdge(int64_t edge_id) {
    if (!IsValidEdge(edge_id)) { return -1; }

    int64_t prev_edge = Backward(edge_id);
    uint8_t c = GetW(prev_edge);
    int64_t ret = IsValidEdge(prev_edge) ? prev_edge : -1;
    int count_ones = IsLastOrTip(prev_edge);

    for (++prev_edge; count_ones < 5 && prev_edge < this->size; ++prev_edge) {
        count_ones += IsLastOrTip(prev_edge);
        uint8_t cur_char = GetW(prev_edge);
        if (cur_char == c) {
            break;
        } else if (cur_char == c + 4 && IsValidEdge(prev_edge)) {
            if (ret != -1) {
                return -1;
            } else {
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
    } else {
        return -1;
    }
}

int64_t SuccinctDBG::NextSimplePathEdge(int64_t edge_id) {
    int64_t next_edge = UniqueNextEdge(edge_id);
    if (next_edge != -1 && UniquePrevEdge(next_edge) != -1) {
        return next_edge;
    } else {
        return -1;
    }
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
            if (IsTip(y)) {
                uint32_t *tip_node_seq = tip_node_seq_ + (size_t)uint32_per_tip_nodes_ * (rs_is_tip_.Rank(y) - 1);
                for (int j = 0; j < i; ++j) {
                    uint8_t c = (tip_node_seq[j / kCharsPerUint32] >> (kCharsPerUint32 - 1 - j % kCharsPerUint32) * kBitsPerChar) & 3;
                    c++;
                    if (c < seq[i - j]) {
                        cmp = -1;
                        break;
                    } else if (c > seq[i - j]) {
                        cmp = 1;
                        break;
                    }
                }

                if (cmp == 0) {
                    if(IsTip(mid)) {
                        cmp = -1;
                    } else {
                        uint8_t c = (tip_node_seq[i / kCharsPerUint32] >> (kCharsPerUint32 - 1 - i % kCharsPerUint32) * kBitsPerChar) & 3;
                        c++;
                        if (c < seq[0]) {
                            cmp = -1;
                            break;
                        } else if (c > seq[0]) {
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
            } else if (c > seq[i]) {
                cmp = 1;
                break;
            }
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
    if (seq[kmer_k] > 4) { seq[kmer_k] -= 4; }

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
    } while (rev_node >= 0 && !IsLastOrTip(rev_node));

    return -1;
}

void SuccinctDBG::LoadFromFile(const char *dbg_name) {
    FILE *w_file = OpenFileAndCheck((std::string(dbg_name) + ".w").c_str(), "rb");
    FILE *last_file = OpenFileAndCheck((std::string(dbg_name) + ".last").c_str(), "rb");
    FILE *f_file = OpenFileAndCheck((std::string(dbg_name) + ".f").c_str(), "r");
    FILE *is_tip = OpenFileAndCheck((std::string(dbg_name) + ".isd").c_str(), "rb");
    FILE *tip_node_seq_file = OpenFileAndCheck((std::string(dbg_name) + ".dn").c_str(), "rb");
    FILE *edge_multiplicity_file = OpenFileAndCheck((std::string(dbg_name) + ".mul").c_str(), "rb");
    FILE *edge_multiplicity_file2 = OpenFileAndCheck((std::string(dbg_name) + ".mul2").c_str(), "rb");

    for (int i = 0; i < kAlphabetSize + 2; ++i) {
        assert(fscanf(f_file, "%lld", &f_[i]) == 1);
    }
    assert(fscanf(f_file, "%d", &kmer_k) == 1);
    assert(fscanf(f_file, "%u", &num_tip_nodes_) == 1);
    size = f_[kAlphabetSize + 1];

    size_t word_needed_w = (size + kWCharsPerWord - 1) / kWCharsPerWord;
    size_t word_needed_last = (size + kBitsPerULL - 1) / kBitsPerULL;
    w_ = (unsigned long long*) MallocAndCheck(sizeof(unsigned long long) * word_needed_w, __FILE__, __LINE__);
    last_ = (unsigned long long*) MallocAndCheck(sizeof(unsigned long long) * word_needed_last, __FILE__, __LINE__);
    is_tip_ = (unsigned long long*) MallocAndCheck(sizeof(unsigned long long) * word_needed_last, __FILE__, __LINE__);
    invalid_ = (unsigned long long*) MallocAndCheck(sizeof(unsigned long long) * word_needed_last, __FILE__, __LINE__);
    edge_multiplicities_ = (multi2_t*) MallocAndCheck(sizeof(multi2_t) * size, __FILE__, __LINE__);
    large_multi_h_ = kh_init(k64v16);

    size_t word_read = fread(w_, sizeof(unsigned long long), word_needed_w, w_file);
    assert(word_read == word_needed_w);
    word_read = fread(last_, sizeof(unsigned long long), word_needed_last, last_file);
    assert(word_read == word_needed_last);
    word_read = fread(is_tip_, sizeof(unsigned long long), word_needed_last, is_tip);
    assert(word_read == word_needed_last);
    memcpy(invalid_, is_tip_, sizeof(unsigned long long) * word_needed_last);
    word_read = fread(edge_multiplicities_, sizeof(multi2_t), size, edge_multiplicity_file);
    assert(word_read == (size_t)size);

    // read large multiplicities
    int64_t *buf = (int64_t*) MallocAndCheck(sizeof(int64_t) * 4096, __FILE__, __LINE__);
    while ((word_read = fread(buf, sizeof(int64_t), 4096, edge_multiplicity_file2)) != 0) {
        for (unsigned i = 0; i < word_read; ++i) {
            int ret;
            khint_t k = kh_put(k64v16, large_multi_h_, buf[i] >> 16, &ret);
            kh_value(large_multi_h_, k) = buf[i] & ((1 << 16) - 1);
        }
    }
    free(buf);

    // read tip nodes sequences
    assert(fread(&uint32_per_tip_nodes_, sizeof(uint32_t), 1, tip_node_seq_file) == 1);
    tip_node_seq_ = (uint32_t*) MallocAndCheck(sizeof(uint32_t) * num_tip_nodes_ * uint32_per_tip_nodes_, __FILE__, __LINE__);

    word_read = fread(tip_node_seq_, sizeof(uint32_t), (size_t)num_tip_nodes_ * uint32_per_tip_nodes_, tip_node_seq_file);
    assert(word_read == (size_t)num_tip_nodes_ * uint32_per_tip_nodes_);

    rs_is_tip_.Build(is_tip_, size);

    init(w_, last_, f_, size, kmer_k);
    need_to_free_ = true;

    fclose(w_file);
    fclose(last_file);
    fclose(f_file);
    fclose(is_tip);
    fclose(tip_node_seq_file);
    if (edge_multiplicity_file != NULL) {
        fclose(edge_multiplicity_file);
        fclose(edge_multiplicity_file2);
    }
}

void SuccinctDBG::LoadFromMultiFile(const char *dbg_name) {
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

    size_t word_needed_w = (size + kWCharsPerWord - 1) / kWCharsPerWord;
    size_t word_needed_last = (size + kBitsPerULL - 1) / kBitsPerULL;

    w_ = (unsigned long long*) MallocAndCheck(sizeof(unsigned long long) * word_needed_w, __FILE__, __LINE__);
    last_ = (unsigned long long*) MallocAndCheck(sizeof(unsigned long long) * word_needed_last, __FILE__, __LINE__);
    is_tip_ = (unsigned long long*) MallocAndCheck(sizeof(unsigned long long) * word_needed_last, __FILE__, __LINE__);
    edge_multiplicities_ = (multi2_t*) MallocAndCheck(sizeof(multi2_t) * size, __FILE__, __LINE__);
    tip_node_seq_ = (uint32_t*) MallocAndCheck(sizeof(uint32_t) * num_tip_nodes_ * sdbg_reader.words_per_tip_label(), __FILE__, __LINE__);
    large_multi_h_ = kh_init(k64v16);

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

        edge_multiplicities_[i] = item >> 8;

        if (UNLIKELY((item >> 8) == kMulti2Sp)) {
            int ret;
            khint_t k = kh_put(k64v16, large_multi_h_, i, &ret);
            kh_value(large_multi_h_, k) = sdbg_reader.NextLargeMul();
        }

        if (UNLIKELY((item >> 5) & 1)) {
            assert(sdbg_reader.NextTipLabel(tip_node_seq_ + tip_label_offset));
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

    invalid_ = (unsigned long long*) MallocAndCheck(sizeof(unsigned long long) * word_needed_last, __FILE__, __LINE__);
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