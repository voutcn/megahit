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

#ifndef SDBG_BUILDER_WRITERS_H_
#define SDBG_BUILDER_WRITERS_H_

#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdint.h>

#include "mem_file_checker-inl.h"

struct DBG_BinaryWriter {
    static const int kBufferSize = 4096;
    static const int kCharsPerWWord = 16;
    static const int kBitsPerWWord = 64;
    uint64_t _w_buffer[kBufferSize];
    uint64_t _last_buffer[kBufferSize];
    uint64_t _is_dollar_buffer[kBufferSize];
    uint64_t _is_low_cover_buffer[kBufferSize];

    int _w_index, _last_index, _is_dollar_index;
    int _w_index_in_word, _last_index_in_word, _is_dollar_index_in_word;
    FILE *_w_file, *_last_file, *_is_dollar_file;

    DBG_BinaryWriter(): _w_file(NULL), _last_file(NULL), _is_dollar_file(NULL) {}

    ~DBG_BinaryWriter() {
        if ((_w_index_in_word > 0) + _w_index > 0) {
            fwrite(_w_buffer, sizeof(uint64_t), (_w_index_in_word > 0) + _w_index, _w_file);
        }
        if ((_last_index_in_word > 0) + _last_index > 0) {
            fwrite(_last_buffer, sizeof(uint64_t), (_last_index_in_word > 0) + _last_index, _last_file);
        }
        if ((_is_dollar_index_in_word > 0) + _is_dollar_index > 0) {
            fwrite(_is_dollar_buffer, sizeof(uint64_t), (_is_dollar_index_in_word > 0) + _is_dollar_index, _is_dollar_file);
        }
        if (_w_file != NULL)
            fclose(_w_file);
        if (_last_file != NULL)
            fclose(_last_file);
        if (_is_dollar_file != NULL)
            fclose(_is_dollar_file);
    }

    void init(const char *w_file_name, const char *last_file_name, const char *is_dollar_file_name) {
        _w_file = OpenFileAndCheck(w_file_name, "wb");
        if (_w_file == NULL) {
            fprintf(stderr, "Open %s failed\n", w_file_name);
            exit(1);
        }

        _last_file = OpenFileAndCheck(last_file_name, "wb");
        if (_last_file == NULL) {
            fprintf(stderr, "Open %s failed\n", last_file_name);
            exit(1);
        }

        _is_dollar_file = OpenFileAndCheck(is_dollar_file_name, "wb");
        if (_is_dollar_file == NULL) {
            fprintf(stderr, "Open %s failed\n", is_dollar_file_name);
            exit(1);
        }

        memset(_w_buffer, 0, sizeof(_w_buffer));
        memset(_last_buffer, 0, sizeof(_last_buffer));
        memset(_is_dollar_buffer, 0, sizeof(_is_dollar_buffer));

        _w_index = 0;
        _w_index_in_word = 0;
        _last_index = 0;
        _last_index_in_word = 0;
        _is_dollar_index = 0;
        _is_dollar_index_in_word = 0;
    }

    void outputW(int c) {
        _w_buffer[_w_index] |= ((uint64_t)c << (4 * _w_index_in_word));
        _w_index_in_word++;
        if (_w_index_in_word >= kCharsPerWWord) {
            _w_index_in_word = 0;
            ++_w_index;
            if (_w_index >= kBufferSize) {
                fwrite(_w_buffer, sizeof(uint64_t), kBufferSize, _w_file);
                memset(_w_buffer, 0, sizeof(_w_buffer));
                _w_index = 0;
            }
        }
    }

    void _outputOneBit(uint64_t *buffer, FILE *file, int &index, int &index_in_word, int c) {
        buffer[index] |= ((uint64_t)c << index_in_word);
        index_in_word++;
        if (index_in_word >= kBitsPerWWord) {
            index_in_word = 0;
            ++index;
            if (index >= kBufferSize) {
                fwrite(buffer, sizeof(buffer[0]), kBufferSize, file);
                memset(buffer, 0, sizeof(buffer[0]) * kBufferSize);
                index = 0;
            }
        }
    }

    void outputLast(int c) {
        _outputOneBit(_last_buffer, _last_file, _last_index, _last_index_in_word, c);
    }

    void outputIsDollar(int c) {
        _outputOneBit(_is_dollar_buffer, _is_dollar_file, _is_dollar_index, _is_dollar_index_in_word, c);
    }

    void outputC(int c, long long number) {
        while (number > 0 && (_w_index_in_word != 0 || _w_index != 0)) {
            outputW(c);
            --number;
        }

        if (number >= kBufferSize * kCharsPerWWord) {
            memset(_w_buffer, (c << 4) | c, sizeof(_w_buffer));
            while (number >= kBufferSize * kCharsPerWWord) {
                fwrite(_w_buffer, sizeof(uint64_t), kBufferSize, _w_file);
                number -= kBufferSize * kCharsPerWWord;
            }
            memset(_w_buffer, 0, sizeof(_w_buffer));
        }

        while (number > 0) {
            outputW(c);
            --number;
        }
    }

    void _outputOnes(uint64_t *buffer, FILE *file, int &index, int &index_in_word, long long number) {
        while (number > 0 && (index_in_word != 0 || index != 0)) {
            _outputOneBit(buffer, file, index, index_in_word, 1);
            --number;
        }

        if (number >= kBufferSize * kBitsPerWWord) {
            memset(buffer, 0xFF, sizeof(buffer[0]) * kBufferSize);
            while (number >= kBufferSize * kBitsPerWWord) {
                fwrite(buffer, sizeof(uint64_t), kBufferSize, file);
                number -= kBufferSize * kBitsPerWWord;
            }
            memset(buffer, 0, sizeof(buffer[0]) * kBufferSize);
        }

        while (number > 0) {
            _outputOneBit(buffer, file, index, index_in_word, 1);
            --number;
        }
    }

    void outputOnesToLast(int64_t number) {
        _outputOnes(_last_buffer, _last_file, _last_index, _last_index_in_word, number);
    }

    void outputOnesToIsDollar(int64_t number) {
        _outputOnes(_is_dollar_buffer, _is_dollar_file, _is_dollar_index, _is_dollar_index_in_word, number);
    }
};

struct WordWriter {
    static const int kBufferSize = 4096;
    int buffer_index;
    uint32_t output_buffer[kBufferSize];
    FILE *file;

    WordWriter() {
        file = NULL;
    }
    ~WordWriter() {
        if (file != NULL)
            fprintf(stderr, "%p\n", file);
        destroy();
    }

    void init(const char *file_name) {
        file = OpenFileAndCheck(file_name, "wb");
        buffer_index = 0;
        assert(file != NULL);
    }

    void destroy() {
        if (file != NULL) {
            if (buffer_index > 0) {
                fwrite(output_buffer, sizeof(uint32_t), buffer_index, file);
            }
            fclose(file);
            file = NULL;
        }
    }

    void output(uint32_t w) {
        output_buffer[buffer_index++] = w;
        if (buffer_index == kBufferSize) {
            fwrite(output_buffer, sizeof(uint32_t), kBufferSize, file);
            buffer_index = 0;
        }
    }
};

#endif // SDBG_BUILDER_WRITERS_H_