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

#ifndef SAFE_ALLOC_OPEN_INL_H
#define SAFE_ALLOC_OPEN_INL_H

#include <stdio.h>
#include <stdlib.h>

inline FILE *xfopen(const char *filename, const char *mode) {
  FILE *fp;
  if ((fp = fopen(filename, mode)) == nullptr) {
    pfprintf(stderr, "[ERROR] Cannot open {s}. Now exit to system...\n", filename);
    exit(-1);
  }
  return fp;
}

#endif  // SAFE_ALLOC_OPEN_INL_H