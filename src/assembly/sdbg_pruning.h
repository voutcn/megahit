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

#ifndef ASSEMBLY_ALGORITHMS_H_
#define ASSEMBLY_ALGORITHMS_H_

#include <string>
#include <vector>
#include "sdbg/sdbg.h"

using std::string;
using std::vector;

namespace sdbg_pruning {

double InferMinDepth(SDBG &dbg);

// tips removal
uint64_t RemoveTips(SDBG &dbg, int max_tip_len);
}  // namespace sdbg_pruning

#endif  // define ASSEMBLY_ALGORITHMS_H_