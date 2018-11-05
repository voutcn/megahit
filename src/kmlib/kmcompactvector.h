//
// Created by vout on 11/4/18.
//

#ifndef KMLIB_COMPACTVECTOR_H
#define KMLIB_COMPACTVECTOR_H

#include <cstdint>
#include <vector>

namespace kmlib {

template<uint8_t BaseSize, typename WordType = uint64_t>
class CompactVector {
  static_assert(sizeof(WordType) * 8 % BaseSize == 0,
                "BaseSize must of power of 2 and no larger than WordType");
 public:
  class Adapter {
   private:
    WordType *word_;
    uint8_t bit_offset_;
   public:
    Adapter(WordType *word, uint8_t base_offset)
        : word_(word), bit_offset_(base_offset * BaseSize) {}
    Adapter& operator=(const WordType &val) {
      *word_ &= ~(kBaseMask << bit_offset_);
      *word_ |= val << bit_offset_;
      return *this;
    }
    Adapter& operator=(const Adapter &rhs) {
      *this = WordType(rhs);
      return *this;
    }
    operator WordType() const {
      return *word_ >> bit_offset_ & kBaseMask;
    }
  };
  CompactVector(size_t size = 0) : size_(size), data_(ComputeWordSize(size)) {}
  ~CompactVector() = default;
  size_t size() const {
    return size_;
  }
  size_t capacity() const {
    return data_.capacity() * kBasePerWord;
  }
  WordType const *data() const {
    return data_.data();
  }
  WordType *data() {
    return data_.data();
  }
  WordType operator[](size_t i) const {
    return data_[i / kBasePerWord] >> (i % kBasePerWord * BaseSize) & kBaseMask;
  }
  Adapter operator[](size_t i) {
    return Adapter(&data_[i / kBasePerWord], i % kBasePerWord);
  }
  void reserve(size_t size) {
    data_.reserve(ComputeWordSize(size));
  }
  void resize(size_t size) {
    size_ = size;  // WARNING bug if size < size_
    data_.resize(ComputeWordSize(size), 0);
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
  static const uint8_t kBasePerWord = sizeof(WordType) * 8 / BaseSize;
  size_t size_;
  std::vector<WordType> data_;
 private:
  static size_t ComputeWordSize(size_t size) {
    return (size + kBasePerWord - 1) / kBasePerWord;
  }
};

} // namespace kmlib

#endif //KMLIB_COMPACTVECTOR_H
