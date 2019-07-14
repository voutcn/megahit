//
// Created by vout on 11/5/18.
//

#ifndef MEGAHIT_BUFFERED_READER_H
#define MEGAHIT_BUFFERED_READER_H

#include <cstring>
#include <fstream>

/**
 * A buffered wrapper to speed up ifstream::read
 */
class BufferedReader {
 public:
  static constexpr size_t kMaxBufferSize = 65536;
  explicit BufferedReader() = default;
  void reset(std::istream *is, size_t buffer_size = kMaxBufferSize) {
    is_ = is;
    head_ = tail_ = 0;
    buffer_size_ = std::min(buffer_size, kMaxBufferSize * 1);
  }

  template <typename T>
  size_t read(T *dst, size_t size = 1) {
    if (is_ == nullptr) {
      return 0;
    }
    size_t wanted = sizeof(T) * size;
    size_t remained = wanted;
    auto dst_ptr = reinterpret_cast<char *>(dst);
    while (remained > 0) {
      if (remained <= tail_ - head_) {
        memcpy(dst_ptr, buffer_ + head_, remained);
        head_ += remained;
        return wanted;
      } else {
        if (tail_ > head_) {
          memcpy(dst_ptr, buffer_ + head_, tail_ - head_);
          remained -= tail_ - head_;
          dst_ptr += tail_ - head_;
        }
        if (refill() == 0) {
          return wanted - remained;
        }
      }
    }
    return 0;
  }

 private:
  size_t refill() {
    head_ = 0;
    is_->read(buffer_, buffer_size_);
    tail_ = *is_ ? buffer_size_ : is_->gcount();
    return tail_;
  }

 private:
  std::istream *is_{};
  char buffer_[kMaxBufferSize]{};
  size_t buffer_size_{kMaxBufferSize};
  size_t head_{0};
  size_t tail_{0};
};

#endif  // MEGAHIT_BUFFERED_READER_H
