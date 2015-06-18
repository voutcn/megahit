/**
 * @file sequence.h
 * @brief Sequence Class. 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-02
 */

#ifndef __SEQUENCE_SEQUENCE_H_

#define __SEQUENCE_SEQUENCE_H_

#include <stdint.h>

#include <istream>
#include <ostream>
#include <string>

#include "lib_idba/kmer.h"

/**
 * @brief It represents a DNA sequence ({A, C, G, T, N}) as a digit sequence 
 * ({0, 1, 2, 3, 4}). 
 */
class Sequence
{
public:
    friend std::istream &operator >>(std::istream &stream, Sequence &seq);
    friend std::ostream &operator <<(std::ostream &stream, const Sequence &seq);

    Sequence() {}
    Sequence(const Sequence &seq, int offset = 0, size_t length = std::string::npos)
    { Assign(seq, offset, length); }
    explicit Sequence(const std::string &seq, int offset = 0, size_t length = std::string::npos)
    { Assign(seq, offset, length); }
    Sequence(uint32_t num, uint8_t ch)
    { Assign(num, ch); }
    explicit Sequence(const IdbaKmer &kmer)
    { Assign(kmer); }

    ~Sequence() {}

    const Sequence &operator =(const Sequence &seq) { Assign(seq); return *this; }
    const Sequence &operator =(const std::string &seq) { Assign(seq); return *this; }
    const Sequence &operator =(const IdbaKmer &kmer) { Assign(kmer); return *this; }

    const Sequence &operator +=(const Sequence &seq) { Append(seq); return *this; }
    const Sequence &operator +=(uint8_t ch) { Append(ch); return *this; }

    bool operator ==(const Sequence &seq) const { return bases_ == seq.bases_; }
    bool operator !=(const Sequence &seq) const { return bases_ != seq.bases_; }
    bool operator <(const Sequence &seq) const { return bases_ < seq.bases_; }
    bool operator >(const Sequence &seq) const { return bases_ > seq.bases_; }

    const Sequence &Assign(const Sequence &seq, int offset = 0, size_t length = std::string::npos)
    { if (&seq != this) bases_.assign(seq.bases_, offset, length); return *this; }
    const Sequence &Assign(const std::string &s, int offset = 0, size_t length = std::string::npos)
    { bases_.assign(s, offset, length); Encode(); return *this; }
    const Sequence &Assign(uint32_t num, uint8_t ch)
    { bases_.assign(num, ch); return *this; }
    const Sequence &Assign(const IdbaKmer &kmer);

    const Sequence &Append(const Sequence &seq, int offset = 0, size_t length = std::string::npos) 
    { bases_.append(seq.bases_, offset, length); return *this; }
    const Sequence &Append(const std::string &seq, int offset = 0, size_t length = std::string::npos) 
    { bases_.append(seq, offset, length); return *this; }
    const Sequence &Append(uint32_t num, uint8_t ch)
    { bases_.append(num, ch); return *this; }
    const Sequence &Append(uint8_t ch)
    { bases_.append(1, ch); return *this; }

    const Sequence &ReverseComplement();
    bool IsValid() const; 
    bool IsPalindrome() const;
    void TrimN();

    IdbaKmer GetIdbaKmer(uint32_t offset, uint32_t kmer_size) const;

    uint8_t &operator [](unsigned index) { return (uint8_t &)bases_[index]; }
    const uint8_t &operator [](unsigned index) const { return (uint8_t &)bases_[index]; }
    uint8_t get_base(uint32_t index) const { return (uint8_t)bases_[index]; }
    void set_base(uint32_t index, uint8_t ch) { bases_[index] = ch; }

    void swap(Sequence &seq) 
    { if (this != &seq) bases_.swap(seq.bases_); }

    uint32_t size() const { return bases_.size(); }
    void resize(int new_size) { bases_.resize(new_size); }
    bool empty() const { return bases_.size() == 0; }

    void clear() { bases_.clear(); }

protected:
    void Encode();
    void Decode();

    bool IsValid(char ch) const
    { return ch == 'A' || ch ==  'C' || ch == 'G' || ch == 'T' || ch == 0 || ch == 1 || ch == 2 || ch == 3; }

private:
    std::string bases_;
};

namespace std
{
template <> inline void swap(Sequence &seq1, Sequence &seq2)
{ seq1.swap(seq2); }
}

std::istream &ReadFasta(std::istream &is, Sequence &seq, std::string &comment);
std::ostream &WriteFasta(std::ostream &os, const Sequence &seq, const std::string &comment);

std::istream &ReadFastq(std::istream &is, Sequence &seq, std::string &comment, std::string &quality);
std::ostream &WriteFastq(std::ostream &os, const Sequence &seq, const std::string &comment, const std::string &quality);


#endif

