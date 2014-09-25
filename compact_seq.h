/*
 *  MEGAHIT
 *  Copyright (C) 2014 The University of Hong Kong
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

#ifndef COMPACT_SEQ_H_
#define COMPACT_SEQ_H_

#include <assert.h>
#include <stdint.h>
#include <string>
#include <cstring>

class CompactSequence
{
public:
    CompactSequence() : length_(0), capacity_(0), data_(NULL) { }
    CompactSequence(const CompactSequence &rhs) : length_(rhs.length_) {
        uint32_t num_uint8 = (length_ + kCharsPerUint8 - 1) / kCharsPerUint8;
        capacity_ = num_uint8 * kCharsPerUint8;
        data_ = (uint8_t*) malloc(sizeof(uint8_t) * num_uint8);
        assert(data_ != NULL);
        memcpy(data_, rhs.data_, sizeof(uint8_t) * num_uint8);
    }
    CompactSequence(CompactSequence&& rhs) : CompactSequence() {
        std::swap(length_, rhs.length_);
        std::swap(capacity_, rhs.capacity_);
        std::swap(data_, rhs.data_);
    }

    ~CompactSequence() {
        if (data_ != NULL) {
            free(data_);
        }
    }

    CompactSequence& operator= (const CompactSequence &rhs) {
        ChangeSize_(rhs.length_);
        uint32_t new_num_uint8 = (length_ + kCharsPerUint8 - 1) / kCharsPerUint8;
        memcpy(data_, rhs.data_, sizeof(uint8_t) * new_num_uint8);
        return *this;
    }

    CompactSequence& operator= (CompactSequence&& rhs) {
        std::swap(length_, rhs.length_);
        std::swap(capacity_, rhs.capacity_);
        std::swap(data_, rhs.data_);
        return *this;
    }

    void resize(uint32_t new_size) {
        if (capacity_ < new_size) {
            uint32_t new_num_uint8 = (new_size + kCharsPerUint8 - 1) / kCharsPerUint8;
            uint8_t *new_data = (uint8_t*) malloc(sizeof(uint8_t) * new_num_uint8);
            assert(new_data != NULL);
            if (data_ != NULL) {
                memcpy(new_data, data_, sizeof(uint8_t) * (length_ + kCharsPerUint8 - 1) / kCharsPerUint8);
                free(data_);
            }
            data_ = new_data;
            capacity_ = new_num_uint8 * kCharsPerUint8;
        }

        length_ = new_size;
    }

    CompactSequence& operator= (const std::string &s) {
        ChangeSize_(s.length());
        memset(data_, 0, sizeof(uint8_t) * (length_ + kCharsPerUint8 - 1) / kCharsPerUint8);
        for (uint32_t i = 0; i < length_; ++i) {
            data_[i / kCharsPerUint8] |= (s[i] << (i % kCharsPerUint8) * kBitsPerChar);
        }

        return *this;
    }

    void ToString(std::string &s) const {
        s.clear();
        s.reserve(length_);
        for (uint32_t i = 0; i < length_; ++i) {
            s.push_back((data_[i / kCharsPerUint8] >> (i % kCharsPerUint8) * kBitsPerChar) & 3);
        }
    }

    void AppendToString(std::string &s) const {
        s.reserve(s.length() + length_);
        for (uint32_t i = 0; i < length_; ++i) {
            s.push_back((data_[i / kCharsPerUint8] >> (i % kCharsPerUint8) * kBitsPerChar) & 3);
        }
    }

    void destory() {
        if (data_ != NULL) {
            free(data_);
            data_ = NULL;
        }
        length_ = 0;
        capacity_ = 0;
    }

    uint32_t length() { return length_; }

private:
    void ChangeSize_(uint32_t new_size) {
        if (capacity_ < new_size) {
            uint32_t new_num_uint8 = (new_size + kCharsPerUint8 - 1) / kCharsPerUint8;
            uint8_t *new_data = (uint8_t*) malloc(sizeof(uint8_t) * new_num_uint8);
            assert(new_data != NULL);
            if (data_ != NULL) {
                free(data_);
            }
            data_ = new_data;
            capacity_ = new_num_uint8 * kCharsPerUint8;
        }

        length_ = new_size;
    }

private:
    static const int kBitsPerUint8 = 8;
    static const int kCharsPerUint8 = 4;
    static const int kBitsPerChar = 2;
    uint32_t length_;
    uint32_t capacity_;
    uint8_t *data_;
};

#endif