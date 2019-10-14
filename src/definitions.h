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

#ifndef MEGAHIT_DEFINITIONS_H
#define MEGAHIT_DEFINITIONS_H

#include <stdint.h>

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "v1.2.9"
#endif

#include "sdbg/sdbg_def.h"

static const unsigned kBitsPerEdgeWord = 32;
static const unsigned kBitsPerEdgeChar = 2;
static const unsigned kCharsPerEdgeWord = 16;
static const unsigned kEdgeCharMask = 0x3;

namespace contig_flag {

static const unsigned kStandalone = 0x1;
static const unsigned kLoop = 0x2;

}  // namespace contig_flag

static const int kUint32PerKmerMaxK = (kMaxK + 1 + 15) / 16;
static const int kUint64PerIdbaKmerMaxK = (kMaxK * 2 + 16 + 63) / 64;

#include "sequence/kmer.h"
typedef Kmer<kUint32PerKmerMaxK, uint32_t> GenericKmer;

#endif  // MEGAHIT_DEFINITIONS_H
