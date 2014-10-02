/**
 * @file compact_sequence.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-02
 * @modified by Dinghua Li
 * @date 2014-10-02
 */

#include "compact_sequence.h"

#include <string>
#include <algorithm>
#include "bit_operation.h"

#include <bitset>
#include <iostream>

using namespace std;

const CompactSequence &CompactSequence::Append(const CompactSequence &compact_seq, int offset, size_t length)
{
    if (length == std::string::npos || length > compact_seq.size() - offset)
        length = compact_seq.size() - offset;

    uint32_t old_size = size();
    resize(old_size + length);

    for (unsigned i = 0; i < length; ++i)
        set_base(i + old_size, compact_seq[i + offset]);
    return *this;
}

const CompactSequence &CompactSequence::Append(const std::string &seq, int offset, size_t length)
{
    if (length == std::string::npos || length > seq.size() - offset)
        length = seq.size() - offset;

    uint32_t old_size = size();
    resize(old_size + length);

    for (unsigned i = 0; i < length; ++i)
        set_base(i + old_size, seq[i + offset]);
    return *this;
}

const CompactSequence &CompactSequence::Append(uint8_t ch)
{
    resize(size() + 1);
    set_base(size()-1, ch);
    return *this;
}

const CompactSequence &CompactSequence::ReverseComplement()
{
    uint8_t remain_bases = data_[data_.size()-1];
    reverse(data_.begin(), data_.begin() + data_.size() - 1);
    for (unsigned i = 0; i+1 < data_.size(); ++i)
        data_[i] = bit_operation::ReverseComplement(data_[i]);
    for (unsigned i = 0; i+2 < data_.size(); ++i)
        data_[i] = (uint8_t(data_[i]) >> (remain_bases<<1)) | (uint8_t(data_[i+1]) << ((4-remain_bases)<<1));
    data_[data_.size()-2] >>= (remain_bases << 1);
    return *this;
}

const CompactSequence &CompactSequence::Reverse()
{
    uint8_t remain_bases = data_[data_.size()-1];
    reverse(data_.begin(), data_.begin() + data_.size() - 1);
    for (unsigned i = 0; i+1 < data_.size(); ++i)
        data_[i] = bit_operation::Reverse(data_[i]);
    for (unsigned i = 0; i+2 < data_.size(); ++i)
        data_[i] = (uint8_t(data_[i]) >> (remain_bases<<1)) | (uint8_t(data_[i+1]) << ((4-remain_bases)<<1));
    data_[data_.size()-2] >>= (remain_bases << 1);
    return *this;
}

std::string CompactSequence::ToDNAString()
{
    static char acgt[] = "ACGT";
    std::string dna;
    dna.resize(this->size());
    for (unsigned i = 0; i < this->size(); ++i) {
        dna[i] = acgt[this->get_base(i)];
    }
    return dna;
}
