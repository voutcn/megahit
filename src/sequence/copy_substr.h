/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong & L3 Bioinformatics
 * Limited
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

#ifndef PACKED_READS_H
#define PACKED_READS_H

#include <algorithm>
#include <iterator>
#include "definitions.h"
#include "kmlib/kmbit.h"

template <typename Ptr>
struct is_uint32_ptr {
  using value_type = typename std::iterator_traits<Ptr>::value_type;
  static const bool value = std::is_same<value_type, uint32_t>::value;
};

// 'spacing' is the strip length for read-word "coalescing"
/**
 * @brief copy src_read[offset...(offset+num_chars_to_copy-1)] to dest
 *
 * @param dst
 * @param src
 * @param offset
 * @param num_chars_to_copy
 */

template <typename DstIt, typename SrcIt,
          typename std::enable_if<is_uint32_ptr<DstIt>::value, int>::type = 0,
          typename std::enable_if<is_uint32_ptr<SrcIt>::value, int>::type = 0>
inline void CopySubstring(DstIt dst, SrcIt src, unsigned offset,
                          unsigned num_chars_to_copy, uint64_t spacing,
                          unsigned words_per_read,
                          unsigned words_per_substring) {
  unsigned which_word = offset / kCharsPerEdgeWord;
  unsigned word_offset = offset % kCharsPerEdgeWord;
  auto src_ptr = src + which_word;

  if (!word_offset) {  // special case (word aligned), easy
    unsigned limit = std::min(words_per_read - which_word, words_per_substring);
    std::copy(src_ptr, src_ptr + limit, dst);
  } else {  // not word-aligned
    unsigned bit_shift = word_offset * kBitsPerEdgeChar;
    unsigned limit =
        std::min(words_per_read - which_word - 1, words_per_substring);
    assert(words_per_read - which_word - 1 >= 0);

    for (unsigned i = 0; i < limit; ++i) {
      dst[i] = (src_ptr[i] << bit_shift) |
               (src_ptr[i + 1] >> (kBitsPerEdgeWord - bit_shift));
    }
    if (limit != words_per_substring) {
      dst[limit] = src_ptr[limit] << bit_shift;
    }
  }

  {
    // now mask the extra bits (TODO can be optimized)
    unsigned num_bits_to_copy = num_chars_to_copy * kBitsPerEdgeChar;
    unsigned bits_to_clear =
        kBitsPerEdgeWord - num_bits_to_copy % kBitsPerEdgeWord;
    which_word = num_bits_to_copy / kBitsPerEdgeWord;
    auto dst_ptr = dst + which_word * spacing;

    if (bits_to_clear < kBitsPerEdgeWord) {
      *dst_ptr >>= bits_to_clear;
      *dst_ptr <<= bits_to_clear;
    } else if (which_word < words_per_substring) {
      *dst_ptr = 0;
    }

    ++which_word;

    while (which_word < words_per_substring) {  // fill zero
      *(dst_ptr += spacing) = 0;
      which_word++;
    }
  }
}

/**
 * @brief copy the reverse complement of
 * src_read[offset...(offset+num_chars_to_copy-1)] to dest
 *
 * @param dst [description]
 * @param src [description]
 * @param offset [description]
 * @param num_chars_to_copy [description]
 */
template <typename DstIt, typename SrcIt,
          typename std::enable_if<is_uint32_ptr<DstIt>::value, int>::type = 0,
          typename std::enable_if<is_uint32_ptr<SrcIt>::value, int>::type = 0>
inline void CopySubstringRC(DstIt dst, SrcIt src, unsigned offset,
                            unsigned num_chars_to_copy, uint64_t spacing,
                            unsigned words_per_read,
                            unsigned words_per_substring) {
  unsigned which_word = (offset + num_chars_to_copy - 1) / kCharsPerEdgeWord;
  unsigned word_offset = (offset + num_chars_to_copy - 1) % kCharsPerEdgeWord;
  auto dst_ptr = dst;

  if (word_offset == kCharsPerEdgeWord - 1) {  // uint32_t aligned
    unsigned limit = std::min(words_per_substring, which_word + 1);
    for (unsigned i = 0; i < limit; ++i) {
      dst[i] = src[which_word - i];
    }
    for (unsigned i = 0; i < words_per_substring && i <= which_word; ++i) {
      dst[i] = kmlib::bit::ReverseComplement<2>(dst[i]);
    }
  } else {
    unsigned bit_offset =
        (kCharsPerEdgeWord - 1 - word_offset) * kBitsPerEdgeChar;
    unsigned i;
    uint32_t w;

    for (i = 0; i < words_per_substring - 1 && i < which_word; ++i) {
      w = (src[which_word - i] >> bit_offset) |
          (src[which_word - i - 1] << (kBitsPerEdgeWord - bit_offset));
      w = kmlib::bit::ReverseComplement<2>(w);
      *dst_ptr = w;
      dst_ptr += spacing;
    }

    // last word
    w = src[which_word - i] >> bit_offset;

    if (which_word >= i + 1) {
      w |= (src[which_word - i - 1] << (kBitsPerEdgeWord - bit_offset));
    }
    *dst_ptr = kmlib::bit::ReverseComplement<2>(w);
  }

  {
    // now mask the extra bits (TODO can be optimized)
    unsigned num_bits_to_copy = num_chars_to_copy * kBitsPerEdgeChar;
    unsigned bits_to_clear =
        kBitsPerEdgeWord - num_bits_to_copy % kBitsPerEdgeWord;
    which_word = num_bits_to_copy / kBitsPerEdgeWord;
    dst_ptr = dst + which_word * spacing;

    if (bits_to_clear < kBitsPerEdgeWord) {
      *dst_ptr >>= bits_to_clear;
      *dst_ptr <<= bits_to_clear;
    } else if (which_word < words_per_substring) {
      *dst_ptr = 0;
    }

    which_word++;

    while (which_word < words_per_substring) {  // fill zero
      *(dst_ptr += spacing) = 0;
      which_word++;
    }
  }
}

#endif  // PACKED_READS_H
