/**
 * @file compact_sequence.h
 * @brief Compact Sequence Class.  
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-02
 * @modified by Dinghua Li
 * @date 2014-10-02
 */

#ifndef __SEQUENCE_COMPACT_SEQUENCE_H_

#define __SEQUENCE_COMPACT_SEQUENCE_H_

#include <stdint.h>
#include <assert.h>
#include <cstring>
#include <string>


class Sequence;
class ShortSequence;

/**
 * @brief It is a compact format of DNA sequence. It has the same 
 * functionalities as Sequence Class but uses less memory and more 
 * computational cost.
 */
class CompactSequence
{
public:

    CompactSequence() { clear(); }
    CompactSequence(const CompactSequence &compact_seq)
    { data_ = compact_seq.data_; }
    CompactSequence(const CompactSequence &compact_seq, int offset, size_t length)
    { Assign(compact_seq, offset, length); }
    explicit CompactSequence(const std::string &seq, int offset = 0, size_t length = std::string::npos)
    { Assign(seq, offset, length); }

    const CompactSequence &operator =(const CompactSequence &compact_seq) { Assign(compact_seq); return *this; }
    const CompactSequence &operator =(const std::string &seq) { Assign(seq); return *this; }

    const CompactSequence &operator +=(const CompactSequence &compact_seq) { Append(compact_seq); return *this; }
    const CompactSequence &operator +=(const std::string &seq) { Append(seq); return *this; }
    const CompactSequence &operator +=(uint8_t ch) { Append(ch); return *this; }

    bool operator ==(const CompactSequence &seq) const { return data_ == seq.data_; }
    bool operator !=(const CompactSequence &seq) const { return data_ != seq.data_; }
    bool operator <(const CompactSequence &seq) const { return data_ < seq.data_; }
    bool operator >(const CompactSequence &seq) const { return data_ > seq.data_; }

    const CompactSequence &Assign(const CompactSequence &compact_seq, int offset = 0, size_t length = std::string::npos)
    { if (&compact_seq != this) { resize(0); Append(compact_seq, offset, length); } return *this; }
    const CompactSequence &Assign(const std::string &seq, int offset = 0, size_t length = std::string::npos)
    { resize(0); Append(seq, offset, length); return *this; }

    const CompactSequence &Append(const CompactSequence &compact_seq, int offset = 0, size_t length = std::string::npos);
    const CompactSequence &Append(const std::string &seq, int offset = 0, size_t length = std::string::npos);
    const CompactSequence &Append(uint8_t ch);

    const CompactSequence &ReverseComplement();
    const CompactSequence &Reverse();

    std::string ToDNAString();

    uint8_t operator [](uint32_t index) const 
    { return (data_[index>>2] >> ((index&3) << 1)) & 3; }
    uint8_t get_base(uint32_t index) const 
    { return (data_[index>>2] >> ((index&3) << 1)) & 3; }
    void set_base(uint32_t index, uint8_t ch)
    { data_[index>>2] = (data_[index>>2] & ~(3U << ((index&3) << 1))) | (ch&3) << ((index&3) << 1); }

    void swap(CompactSequence &compact_seq) 
    { if (this != &compact_seq) data_.swap(compact_seq.data_); }

    uint32_t size() const 
    { return ((data_.size()-1) << 2) - data_[data_.size()-1]; }
    uint32_t length() const
    { return size(); }
    void resize(int new_size) 
    { data_.resize((new_size + 7) >> 2); data_[data_.size()-1] = ((data_.size() - 1) << 2) - new_size; }
    bool empty() const { return size() == 0; }

    void clear() { resize(0); }

private:
    std::basic_string<uint8_t> data_;
};

namespace std
{
template <> inline void swap(CompactSequence &compact_seq1, CompactSequence &compact_seq2)
{ compact_seq1.swap(compact_seq2); }
}

#endif

