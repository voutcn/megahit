/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong & L3 Bioinformatics Limited
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

#ifndef ASSEMBLY_ALGORITHMS_H_
#define ASSEMBLY_ALGORITHMS_H_

#include <vector>
#include <string>
#include "succinct_dbg.h"

using std::vector;
using std::string;

namespace assembly_algorithms {

// traverse simple path
int64_t NextSimplePathNode(SuccinctDBG &dbg, int64_t curr_node); // return -1 if cannot extend, else the node index (whose last must be 1)
int64_t PrevSimplePathNode(SuccinctDBG &dbg, int64_t curr_node); // return -1 if cannot extend, else the incoming edge index (whose last may not be 1)

double SetMinDepth(SuccinctDBG &dbg);

// tips removal
int64_t Trim(SuccinctDBG &dbg, int len, int min_final_standalone);
int64_t RemoveTips(SuccinctDBG &dbg, int max_tip_len, int min_final_standalone);

// bubble merging
int64_t PopBubbles(SuccinctDBG &dbg, int max_bubble_len, double low_depth_ratio = 1);

// for experiments only
void MarkSubGraph(SuccinctDBG &dbg, const char* seq, int seq_len);

}

#endif // define ASSEMBLY_ALGORITHMS_H_