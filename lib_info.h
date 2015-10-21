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

#ifndef LIB_INFO_H__
#define LIB_INFO_H__

#include <stdint.h>
#include <string>
#include "sequence_package.h"

struct lib_info_t {
    SequencePackage *p;
    int64_t from;
    int64_t to;
    int max_read_len;
    bool is_pe;
    std::string metadata; // raw file names

    lib_info_t(SequencePackage *p = NULL, int64_t from = 0, int64_t to = 0,
               int max_read_len = 0, bool is_pe = false, const std::string &metadata = ""):
        p(p), from(from), to(to), max_read_len(max_read_len), is_pe(is_pe), metadata(metadata) { }
    ~lib_info_t() { }
};

#endif