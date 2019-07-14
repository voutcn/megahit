//
// Created by vout on 11/4/18.
//

#ifndef KMLIB_COMPACTVECTOR_H
#define KMLIB_COMPACTVECTOR_H

#include <vector>

namespace kmlib {

enum { kBigEndian, kLittleEndian };

template <unsigned BaseSize, typename WordType = unsigned long,
          int Endian = kLittleEndian>
class CompactVector {
  static_assert(sizeof(WordType) * 8 % BaseSize == 0,
                "BaseSize must of power of 2 and no larger than TWord");

 public:
  using word_type = WordType;
  using size_type = typename std::vector<word_type>::size_type;
  static const unsigned kBaseSize = BaseSize;
  static const int kEndian = Endian;
  static const word_type kBaseMask = (word_type{1} << kBaseSize) - 1;
  static const unsigned kBasesPerWord = sizeof(word_type) * 8 / kBaseSize;

 public:
  class Adapter {
   private:
    word_type *word_;
    unsigned bit_offset_;

   public:
    Adapter(word_type *word, unsigned base_offset)
        : word_(word), bit_offset_(bit_shift(base_offset)) {}
    Adapter &operator=(const word_type &val) {
      *word_ &= ~(kBaseMask << bit_offset_);
      *word_ |= val << bit_offset_;
      return *this;
    }
    Adapter &operator=(const Adapter &rhs) {
      *this = word_type(rhs);
      return *this;
    }
    operator word_type() const { return *word_ >> bit_offset_ & kBaseMask; }
  };

 public:
  static size_type size_to_word_count(size_type size) {
    return (size + kBasesPerWord - 1) / kBasesPerWord;
  }
  static unsigned bit_shift(unsigned pos, unsigned len = 1) {
    return Endian == kBigEndian
               ? sizeof(word_type) * 8 - (pos + len) * kBaseSize
               : pos * kBaseSize;
  }
  static word_type at(const word_type *v, size_type i) {
    return v[i / kBasesPerWord] >> bit_shift(i % kBasesPerWord) & kBaseMask;
  }
  static word_type sub_word(word_type val, unsigned pos, unsigned len = 1) {
    if (len == kBasesPerWord) {
      return val >> bit_shift(pos, len);
    } else {
      return val >> bit_shift(pos, len) &
             ((word_type{1} << len * kBaseSize) - 1);
    }
  }

 public:
  explicit CompactVector(size_type size = 0)
      : size_(size),
        underlying_vec_(size_to_word_count(size)),
        vec_ptr_(&underlying_vec_) {}
  explicit CompactVector(std::vector<word_type> *v, size_type size = 0)
      : size_(size), underlying_vec_(0), vec_ptr_(v) {
    resize(size);
  }
  ~CompactVector() = default;
  size_type size() const { return size_; }
  size_type word_count() const { return vec_ptr_->size(); }
  void clear() {
    size_ = 0;
    vec_ptr_->clear();
  }
  size_type capacity() const { return vec_ptr_->capacity() * kBasesPerWord; }
  size_type word_capacity() const { return vec_ptr_->capacity(); }
  const word_type *data() const { return vec_ptr_->data(); }
  word_type *data() { return vec_ptr_->data(); }
  word_type operator[](size_type i) const { return at(vec_ptr_->data(), i); }
  Adapter operator[](size_type i) {
    return Adapter(vec_ptr_->data() + i / kBasesPerWord, i % kBasesPerWord);
  }
  void reserve(size_type size) { vec_ptr_->reserve(size_to_word_count(size)); }
  void resize(size_type size) {
    size_type old_size = size_;
    size_ = size;
    vec_ptr_->resize(size_to_word_count(size), 0);
    if (size_ < old_size && size_ > 0 && size_ % kBasesPerWord != 0) {
      vec_ptr_->back() = sub_word(vec_ptr_->back(), 0, size_ % kBasesPerWord)
                         << bit_shift(0, size_ % kBasesPerWord);
    }
  }
  void push_back(word_type val) {
    if (size_ % kBasesPerWord == 0) {
      vec_ptr_->emplace_back(val << bit_shift(0));
    } else {
      vec_ptr_->back() |= val << bit_shift(size_ % kBasesPerWord);
    }
    ++size_;
  }
  void push_word(word_type val, unsigned pos, unsigned len) {
    unsigned pos_in_back = size_ % kBasesPerWord;
    unsigned remaining = kBasesPerWord - pos_in_back;
    if (pos_in_back == 0) {
      vec_ptr_->emplace_back(sub_word(val, pos, len) << bit_shift(0, len));
    } else if (remaining < len) {
      vec_ptr_->back() |= sub_word(val, pos, remaining)
                          << bit_shift(pos_in_back, remaining);
      vec_ptr_->emplace_back(sub_word(val, pos + remaining, len - remaining)
                             << bit_shift(0, len - remaining));
    } else {
      vec_ptr_->back() |= sub_word(val, pos, len)
                          << bit_shift(pos_in_back, len);
    }
    size_ += len;
  }
  void push_word(word_type val, unsigned pos = 0) {
    push_word(val, pos, kBasesPerWord - pos);
  }
  void pop_back() {
    (*this)[size_ - 1] = 0;
    --size_;
  }

 private:
  size_type size_;
  std::vector<word_type> underlying_vec_;
  std::vector<word_type> *vec_ptr_;
};

}  // namespace kmlib

#endif  // KMLIB_COMPACTVECTOR_H
