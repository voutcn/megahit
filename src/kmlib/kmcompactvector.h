//
// Created by vout on 11/4/18.
//

#ifndef KMLIB_COMPACTVECTOR_H
#define KMLIB_COMPACTVECTOR_H

#include <vector>

namespace kmlib {

template<unsigned BaseSize, typename WordType = unsigned long, bool ReverseOrder = false>
class CompactVector {
  static_assert(sizeof(WordType) * 8 % BaseSize == 0,
                "BaseSize must of power of 2 and no larger than WordType");
 public:
  using word_type = WordType;
  using size_type = typename std::vector<word_type>::size_type;
  static const unsigned kBaseSize = BaseSize;
  static const bool kReverseOrder = ReverseOrder;
 public:
  class Adapter {
   private:
    word_type *word_;
    unsigned bit_offset_;
   public:
    Adapter(word_type *word, unsigned base_offset)
        : word_(word), bit_offset_(bit_offset(base_offset)) {}
    Adapter &operator=(const word_type &val) {
      *word_ &= ~(kBaseMask << bit_offset_);
      *word_ |= val << bit_offset_;
      return *this;
    }
    Adapter &operator=(const Adapter &rhs) {
      *this = word_type(rhs);
      return *this;
    }
    operator word_type() const {
      return *word_ >> bit_offset_ & kBaseMask;
    }
  };
 private:
  static size_type size_to_word_count(size_type size) {
    return (size + kBasePerWord - 1) / kBasePerWord;
  }
  static unsigned bit_offset(unsigned base_offset) {
    return kReverseOrder ?
           sizeof(word_type) * 8 - (base_offset + 1) * kBaseSize:
           base_offset * kBaseSize;
  }
 public:
  static word_type at(const word_type *v, size_type i) {
    return v[i / kBasePerWord] >> bit_offset(i % kBasePerWord) & kBaseMask;
  }
 public:
  CompactVector(size_type size = 0) : size_(size), data_(size_to_word_count(size)) {}
  ~CompactVector() = default;
  size_type size() const {
    return size_;
  }
  size_type word_count() const {
    return data_.size();
  }
  size_type capacity() const {
    return data_.capacity() * kBasePerWord;
  }
  const word_type *data() const {
    return data_.data();
  }
  word_type *data() {
    return data_.data();
  }
  word_type operator[](size_type i) const {
    return at(data_.data(), i);
  }
  Adapter operator[](size_type i) {
    return Adapter(&data_[i / kBasePerWord], i % kBasePerWord);
  }
  void reserve(size_type size) {
    data_.reserve(size_to_word_count(size));
  }
  void resize(size_type size) {
    size_ = size;  // WARNING bug if size < item_count_
    data_.resize(size_to_word_count(size), 0);
  }
  void push_back(word_type val) {
    if (size_ == data_.size() * kBasePerWord) {
      data_.emplace_back(val);
    } else {
      data_.back() |= val << (size_ % kBasePerWord * kBaseSize);
    }
    ++size_;
  }
  void pop_back() {
    (*this)[size_ - 1] = 0;
    --size_;
  }
 private:
  static const word_type kBaseMask = (word_type{1} << kBaseSize) - 1;
  static const unsigned kBasePerWord = sizeof(word_type) * 8 / kBaseSize;
  size_type size_;
  std::vector<word_type> data_;
};

} // namespace kmlib

#endif //KMLIB_COMPACTVECTOR_H
