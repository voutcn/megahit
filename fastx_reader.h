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

#ifndef FASTX_READER_H_
#define FASTX_READER_H_

#include <assert.h>
#include <string.h>
#include <string>
#include <zlib.h>

/**
 * @brief ascii reader for gz file
 *
 */
class BufferReader {
public:
    // constructor
    BufferReader(): input_(NULL) { }
    BufferReader(gzFile input): input_(input) {
        Refill();
    }

    void init(gzFile input) {
        input_ = input;
        Refill();
    }

    bool eof() {
        if (*p_) return false;
        Refill();
        return num_bytes_in_buffer_ == 0;
    }

    // go to the next line (skip until the next '\n'
    void SkipLine() {
        p_ = (char*) memchr(p_, '\n', num_bytes_in_buffer_ - (p_ - buffer_));
        while (!p_) { // if no end of line, Refill buffer
            Refill();
            if (!num_bytes_in_buffer_) return; // eof
            p_ = (char*) memchr(p_, '\n', num_bytes_in_buffer_); // try again
        }
        if (!*(++p_)) { // reached end of buffer, need Refill
            Refill();
        }
    }

    char NextChar() {
        if (*p_) return *(p_++);
        Refill();
        return *(p_++);
    }

private:
    // Refill buffer
    void Refill() {
        num_bytes_in_buffer_ = gzread(input_, buffer_, kReaderBufferSize);
        buffer_[num_bytes_in_buffer_] = 0;
        p_ = buffer_;
    }


private:
    static const int kReaderBufferSize = 4096;
    gzFile input_;
    char buffer_[kReaderBufferSize + 1];
    int num_bytes_in_buffer_;
    char* p_; // current pointer
};

class FastxReader {
public:
    enum FastxFormat {
        kFasta,
        kFastq,
        kNull,
    };

    FastxReader(): format_(kNull) {}
    FastxReader(gzFile input): buffer_reader_(input) {
        char first_char = buffer_reader_.NextChar();
        if (first_char == '>') {
            format_ = kFasta;
            read_terminator_ = '>';
        } else if (first_char == '@') {
            format_ = kFastq;
            read_terminator_ = '+';
        } else {
            assert(false);
        }
    }

    void init(gzFile input) {
        buffer_reader_.init(input);
        char first_char = buffer_reader_.NextChar();
        if (first_char == '>') {
            format_ = kFasta;
            read_terminator_ = '>';
        } else if (first_char == '@') {
            format_ = kFastq;
            read_terminator_ = '+';
        } else {
            assert(false);
        }
    }

    bool eof() {
        return buffer_reader_.eof();
    }

    size_t NextSeq(std::string &seq) {
        buffer_reader_.SkipLine();
        seq.clear();

        char c;
        while ((c = buffer_reader_.NextChar()) != read_terminator_ && !buffer_reader_.eof()) {
            if (c >= 'A') {
                // a quick check whether it is a letter
                seq.push_back(c);
            }
        }

        if (format_ == kFastq) {
            buffer_reader_.SkipLine();
            buffer_reader_.SkipLine();
        }

        return seq.size();
    }

private:
    BufferReader buffer_reader_;
    FastxFormat format_;
    char read_terminator_;
};

#endif