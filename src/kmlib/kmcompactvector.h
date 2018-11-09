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
  class Adapter {
   private:
    WordType *word_;
    unsigned bit_offset_;
   public:
    Adapter(WordType *word, unsigned base_offset)
        : word_(word), bit_offset_(bit_offset(base_offset)) {}
    Adapter &operator=(const WordType &val) {
      *word_ &= ~(kBaseMask << bit_offset_);
      *word_ |= val << bit_offset_;
      return *this;
    }
    Adapter &operator=(const Adapter &rhs) {
      *this = WordType(rhs);
      return *this;
    }
    operator WordType() const {
      return *word_ >> bit_offset_ & kBaseMask;
    }
  };
 private:
  static size_t size_to_word_count(size_t size) {
    return (size + kBasePerWord - 1) / kBasePerWord;
  }
  static unsigned bit_offset(unsigned base_offset) {
    return ReverseOrder ?
           sizeof(WordType) * 8 - (base_offset + 1) * BaseSize:
           base_offset * BaseSize;
  }
 public:
  static WordType at(const WordType *v, size_t i) {
    return v[i / kBasePerWord] >> bit_offset(i % kBasePerWord) & kBaseMask;
  }
 public:
  CompactVector(size_t size = 0) : size_(size), data_(size_to_word_count(size)) {}
  ~CompactVector() = default;
  size_t size() const {
    return size_;
  }
  size_t word_count() const {
    return data_.size();
  }
  size_t capacity() const {
    return data_.capacity() * kBasePerWord;
  }
  const WordType *data() const {
    return data_.data();
  }
  WordType *data() {
    return data_.data();
  }
  WordType operator[](size_t i) const {
    return at(data_.data(), i);
  }
  Adapter operator[](size_t i) {
    return Adapter(&data_[i / kBasePerWord], i % kBasePerWord);
  }
  void reserve(size_t size) {
    data_.reserve(size_to_word_count(size));
  }
  void resize(size_t size) {
    size_ = size;  // WARNING bug if size < item_count_
    data_.resize(size_to_word_count(size), 0);
  }
  void push_back(WordType val) {
    if (size_ == data_.size() * kBasePerWord) {
      data_.emplace_back(val);
    } else {
      data_.back() |= val << (size_ % kBasePerWord * BaseSize);
    }
    ++size_;
  }
  void pop_back() {
    (*this)[size_ - 1] = 0;
    --size_;
  }
 private:
  static const WordType kBaseMask = (WordType{1} << BaseSize) - 1;
  static const unsigned kBasePerWord = sizeof(WordType) * 8 / BaseSize;
  size_t size_;
  std::vector<WordType> data_;
};

} // namespace kmlib

#endif //KMLIB_COMPACTVECTOR_H
