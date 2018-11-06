//
// Created by vout on 11/5/18.
//

#ifndef MEGAHIT_BUFFERED_READER_H
#define MEGAHIT_BUFFERED_READER_H

#include <fstream>
#include <cstring>

/**
 * A buffered wrapper to speed up ifstream::read
 */
class BufferedReader {
 public:
  explicit BufferedReader() = default;
  void reset(std::ifstream *is) {
    is_ = is;
    head_ = tail_ = 0;
  }
  template <typename T>
  size_t read(T* dst, size_t size = 1) {
    size_t wanted = sizeof(T) * size;
    size_t remained = wanted;
    auto dst_ptr = reinterpret_cast<char*>(dst);
    while (remained > 0) {
      if (remained <= tail_ - head_) {
        memcpy(dst_ptr, buffer_ + head_, remained);
        head_ += remained;
        return wanted;
      } else {
        memcpy(dst_ptr, buffer_ + head_, tail_ - head_);
        remained -= tail_ - head_;
        dst_ptr += tail_ - head_;
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
    is_->read(buffer_, kBufferSize);
    tail_ = is_ ? kBufferSize : is_->gcount();
    return tail_;
  }
 private:
  static const size_t kBufferSize = 65536;
  std::ifstream *is_{};
  char buffer_[kBufferSize];
  size_t head_{kBufferSize};
  size_t tail_{kBufferSize};
};

#endif //MEGAHIT_BUFFERED_READER_H
