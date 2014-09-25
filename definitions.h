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

#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

#include <stdint.h>

typedef uint16_t multi_t;
const int kBitsPerMulti_t = 16;
const int kMaxMulti_t = 65535;

typedef uint32_t edge_word_t;
static const int kBitsPerEdgeWord = 32;
static const int kBitsPerEdgeChar = 2;
static const int kCharsPerEdgeWord = 16;
static const int kEdgeCharMask = 0x3;
static const int kAlphabetSize = 4;
static const char dna_chars[] = "ACGT";

static const int kMaxNumCPUThreads = 255;

#endif // DEFINITIONS_H_