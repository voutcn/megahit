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


/**
 * Functions for packed reads or edges (ACGT -> 0123, packed by uint32_t)
 */

#ifndef PACKED_READS_H__
#define PACKED_READS_H__

#include <algorithm>
#include "definitions.h"
#include "kseq.h"
#include "mem_file_checker-inl.h"

/*
 * Packs an ASCII read into 2-bit per base form. The last base goes to the MSB of the first word.
 * -Params-
 * read: the read, in ASCII ACGT
 * p: a pointer to the first edge_word_t of the packed sequence to be written to
 * read_length: number of bases in the read
 * last_word_shift: the number of empty bits in the last word. we need to shift the bits up by this amount. solely determined by read_length.
 */

static int dna_map[256];

#ifndef KSEQ_INITED
#define KSEQ_INITED
KSEQ_INIT(gzFile, gzread)
#endif

inline void PackReadFromAscii(char* read, edge_word_t* p, int read_length, int words_per_read) {
    // for de Bruijn graph construction, packing the reverse is more convenient
    edge_word_t w = 0;
    int i, j;
    for (i = 0, j = 0; j < read_length; ++j) {
        if (j % kCharsPerEdgeWord == 0 && j) { // TODO bitwise?
            *(p++) = w;
            w = 0;
            ++i;
        }
        while (read[j] == 'N') {
            break;
        }
        w = (w << kBitsPerEdgeChar) | dna_map[(int)read[j]];
    }

    int last_word_shift = j % kCharsPerEdgeWord;
    last_word_shift = last_word_shift ? (kCharsPerEdgeWord - last_word_shift) * kBitsPerEdgeChar : 0;
    *p = w << last_word_shift;

    while (++i < words_per_read) {
        *(++p) = 0;
    }

    *p |= read_length;
}

inline int64_t ReadFastxAndPack(edge_word_t *&packed_reads, const char *input_file, int words_per_read,
                                int min_read_len, int max_read_length, int64_t max_num_reads) {
    for (int i = 0; i < 10; ++i) {
        dna_map[int("ACGTNacgtn"[i])] = "0123101231"[i] - '0';
    }

    int64_t capacity = std::min(max_num_reads, int64_t(1048576)); // initial capacity 1M
    gzFile fp = strcmp(input_file, "-") ? gzopen(input_file, "r") : gzdopen(fileno(stdin), "r");
    kseq_t *seq = kseq_init(fp); // kseq to read files
    edge_word_t *packed_reads_p; // current pointer
    packed_reads_p = packed_reads = (edge_word_t*) MallocAndCheck(capacity * words_per_read * sizeof(edge_word_t), __FILE__, __LINE__);
    int64_t num_reads = 0;
    int read_length = 0;

    // --- main reading loop ---
    bool stop_reading = false;
    while ((read_length = kseq_read(seq)) >= 0 && !stop_reading) {
        std::reverse(seq->seq.s, seq->seq.s + read_length);
        char *next_p = seq->seq.s;
        while (read_length >= min_read_len) {
            int scan_len = 0;
            while (scan_len < read_length && next_p[scan_len] != 'N') {
                ++scan_len;
            }

            if (scan_len >= min_read_len && scan_len <= max_read_length) {
                if (num_reads >= capacity) {
                    if (capacity == max_num_reads) {
                        err("[%s WRANING] No enough memory to hold all the reads. Only the first %llu reads are kept.\n", __func__, num_reads);
                        stop_reading = true;
                        break;
                    }
                    capacity = std::min(capacity * 2, max_num_reads);
                    edge_word_t *new_ptr = (edge_word_t*) realloc(packed_reads, capacity * words_per_read * sizeof(edge_word_t));
                    if (new_ptr != NULL) {
                        packed_reads = new_ptr;
                        packed_reads_p = packed_reads + words_per_read * num_reads;
                    } else {
                        err("[%s WRANING] No enough memory to hold all the reads. Only the first %llu reads are kept.\n", __func__, num_reads);
                        stop_reading = true;
                        break;
                    }
                }
                // read length is ok! compress and store the packed read
                PackReadFromAscii(next_p, packed_reads_p, scan_len, words_per_read);
                packed_reads_p += words_per_read;
                ++num_reads;
            } else if (scan_len > max_read_length) { // this read length is wrong
                err("[%s WARNING] Found a read of length %d > max read length = %d\n, it will be discarded.", __func__, scan_len, max_read_length);
            }

            while (scan_len < read_length && next_p[scan_len] == 'N') {
                ++scan_len;
            }
            read_length -= scan_len;
            next_p += scan_len;
        }
    }

    kseq_destroy(seq);
    gzclose(fp);

    return num_reads;
}

/**
 * @param packed_reads the ptr of read0
 * @param i read id
 * @param words_per_read
 *
 * @return pointer pointing to the starting pos a read
 */
inline edge_word_t* GetReadPtr(edge_word_t *packed_reads, int64_t i, int words_per_read) {
    return packed_reads + i * words_per_read;
}

/**
 * @brief obtain the read length of a read
 *
 * @param read_p
 * @param words_per_read
 * @param mask 0000...11111 = (1 << bits_of_len) - 1
 * @return read length
 */
inline int GetReadLength(edge_word_t* read_p, int words_per_read, int mask) {
    return *(read_p + words_per_read - 1) & mask;
}

/**
 * @brief extract the nth char in a packed read/edge
 */
inline int ExtractNthChar(edge_word_t *read_ptr, int n) {
    int which_word = n / kCharsPerEdgeWord;
    int index_in_word = n % kCharsPerEdgeWord;
    return (read_ptr[which_word] >> (kBitsPerEdgeChar * (kCharsPerEdgeWord - 1 - index_in_word))) & kEdgeCharMask;
}

// 'spacing' is the strip length for read-word "coalescing"
/**
 * @brief copy src_read[offset...(offset+num_chars_to_copy-1)] to dest
 *
 * @param dest
 * @param src_read
 * @param offset
 * @param num_chars_to_copy
 */
inline void CopySubstring(edge_word_t* dest, edge_word_t* src_read, int offset, int num_chars_to_copy,
                          int64_t spacing, int words_per_read, int words_per_substring) {
    // copy words of the suffix to the suffix pool
    int which_word = offset / kCharsPerEdgeWord;
    int word_offset = offset % kCharsPerEdgeWord;
    edge_word_t *src_p = src_read + which_word;
    edge_word_t *dest_p = dest;
    int num_words_copied = 0;
    if (!word_offset) { // special case (word aligned), easy
        while (which_word < words_per_read && num_words_copied < words_per_substring) {
            *dest_p = *src_p; // write out
            dest_p += spacing;
            src_p++;
            which_word++;
            num_words_copied++;
        }
    } else { // not word-aligned
        int bit_shift = offset * kBitsPerEdgeChar;
        edge_word_t s = *src_p;
        edge_word_t d = s << bit_shift;
        which_word++;
        while (which_word < words_per_read) {
            s = *(++src_p);
            d |= s >> (kBitsPerEdgeWord - bit_shift);
            *dest_p = d; // write out
            if (++num_words_copied >= words_per_substring) goto here;
            dest_p += spacing;
            d = s << bit_shift;
            which_word++;
        }
        *dest_p = d; // write last word
here:
        ;
    }

    {
        // now mask the extra bits (TODO can be optimized)
        int num_bits_to_copy = num_chars_to_copy * 2;
        int which_word = num_bits_to_copy / kBitsPerEdgeWord;
        edge_word_t *p = dest + which_word * spacing;
        int bits_to_clear = kBitsPerEdgeWord - num_bits_to_copy % kBitsPerEdgeWord;
        if (bits_to_clear < kBitsPerEdgeWord) {
            *p >>= bits_to_clear;
            *p <<= bits_to_clear;
        } else if (which_word < words_per_substring) {
            *p = 0;
        }
        which_word++;
        while (which_word < words_per_substring) { // fill zero
            *(p+=spacing) = 0;
            which_word++;
        }
    }
}

/**
 * @brief copy the reverse complement of src_read[offset...(offset+num_chars_to_copy-1)] to dest
 *
 * @param dest [description]
 * @param src_read [description]
 * @param offset [description]
 * @param num_chars_to_copy [description]
 */
inline void CopySubstringRC(edge_word_t* dest, edge_word_t* src_read, int offset, int num_chars_to_copy,
                            int64_t spacing, int words_per_read, int words_per_substring) {
    int which_word = (offset + num_chars_to_copy - 1) / kCharsPerEdgeWord;
    int word_offset = (offset + num_chars_to_copy - 1) % kCharsPerEdgeWord;
    edge_word_t *dest_p = dest;

    if (word_offset == kCharsPerEdgeWord - 1) { // edge_word_t aligned
        for (int i = 0; i < words_per_substring && i <= which_word; ++i) {
            *dest_p = ~ mirror(src_read[which_word - i]);
            dest_p += spacing;
        }
    } else {
        int bit_offset = (kCharsPerEdgeWord - 1 - word_offset) * kBitsPerEdgeChar;
        int i;
        edge_word_t w;
        for (i = 0; i < words_per_substring - 1 && i < which_word; ++i) {
            w = (src_read[which_word - i] >> bit_offset) |
                (src_read[which_word - i - 1] << (kBitsPerEdgeWord - bit_offset));
            *dest_p = ~ mirror(w);
            dest_p += spacing;
        }
        // last word
        w = src_read[which_word - i] >> bit_offset;
        if (which_word >= i + 1) {
            w |= (src_read[which_word - i - 1] << (kBitsPerEdgeWord - bit_offset));
        }
        *dest_p = ~ mirror(w);
    }

    {
        // now mask the extra bits (TODO can be optimized)
        int num_bits_to_copy = num_chars_to_copy * 2;
        int which_word = num_bits_to_copy / kBitsPerEdgeWord;
        edge_word_t *p = dest + which_word * spacing;
        int bits_to_clear = kBitsPerEdgeWord - num_bits_to_copy % kBitsPerEdgeWord;
        if (bits_to_clear < kBitsPerEdgeWord) {
            *p >>= bits_to_clear;
            *p <<= bits_to_clear;
        } else if (which_word < words_per_substring) {
            *p = 0;
        }
        which_word++;
        while (which_word < words_per_substring) { // fill zero
            *(p+=spacing) = 0;
            which_word++;
        }
    }
}

#endif // PACKED_READS_H__