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

#ifndef ATOMIC_BIT_VECTOR_H_
#define ATOMIC_BIT_VECTOR_H_

#include <assert.h>
#include <stdlib.h>
#include <memory.h>
#include <stdint.h>
#include <algorithm>
#include "mem_file_checker-inl.h"

class AtomicBitVector {
  public:
    typedef uint8_t word_t;

    AtomicBitVector(size_t size = 0): size_(size) {
        num_words_ = ((size + kBitsPerWord - 1) / kBitsPerWord);
        capacity_ = num_words_;

        if (num_words_ != 0) {
            data_ = (word_t *) MallocAndCheck(sizeof(word_t) * num_words_, __FILE__, __LINE__);
            memset(data_, 0, sizeof(word_t) * num_words_);
        }
        else {
            data_ = NULL;
        }
    }

    ~AtomicBitVector() {
        if (data_ != NULL) {
            free(data_);
        }
    }

    size_t size() {
        return size_;
    }

    bool get(size_t i) {
        return bool((data_[i / kBitsPerWord] >> i % kBitsPerWord) & 1);
    }

    bool try_lock(size_t i) {
        // assert(i / kBitsPerWord < num_words_);
        word_t *p = data_ + i / kBitsPerWord;

        while (!((*p >> i % kBitsPerWord) & 1)) {
            word_t old_value = *p;
            word_t new_value = old_value | (word_t(1) << (i % kBitsPerWord));

            if (__sync_bool_compare_and_swap(p, old_value, new_value)) {
                return true;
            }
        }

        return false;
    }

    void set(size_t i) {
        __sync_fetch_and_or(data_ + i / kBitsPerWord, word_t(1) << (i % kBitsPerWord));
    }

    void unset(size_t i) {
        // assert(i / kBitsPerWord < num_words_);
        word_t mask = ~(word_t(1) << (i % kBitsPerWord));
        __sync_fetch_and_and(data_ + i / kBitsPerWord, mask);
    }

    void reset(size_t size = 0, int reset_value = 0) {
        size_ = size;
        num_words_ = (size + kBitsPerWord - 1) / kBitsPerWord;

        if (capacity_ < num_words_) {
            word_t *new_data = (word_t *) ReAllocAndCheck(data_, sizeof(word_t) * num_words_, __FILE__, __LINE__);
            data_ = new_data;
            capacity_ = num_words_;
        }

        if (num_words_ != 0) {
            if (reset_value) {
                memset(data_, -1, sizeof(word_t) * num_words_);
            }
            else {
                memset(data_, 0, sizeof(word_t) * num_words_);
            }
        }
    }

    void swap(AtomicBitVector &rhs) {
        if (data_ != rhs.data_) {
            std::swap(data_, rhs.data_);
            std::swap(size_, rhs.size_);
            std::swap(num_words_, rhs.num_words_);
            std::swap(capacity_, rhs.capacity_);
        }
    }

  private:
    static const int kBitsPerByte = 8;
    static const int kBitsPerWord = sizeof(word_t) * kBitsPerByte;
    size_t size_;
    size_t num_words_;
    size_t capacity_;
    word_t *data_;
};

#endif