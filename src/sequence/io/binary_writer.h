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

#ifndef READ_LIB_FUNCTIONS_INL_H
#define READ_LIB_FUNCTIONS_INL_H

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

#include "sequence/sequence_package.h"

inline void WriteBinarySequences(const SeqPackage &pkg, FILE *file, int64_t from = 0, int64_t to = -1) {
  if (to == -1) {
    to = pkg.seq_count() - 1;
  }

  uint32_t len;
  std::vector<uint32_t> s;

  for (int64_t i = from; i <= to; ++i) {
    len = pkg.GetSeqView(i).length();
    pkg.FetchSequence(i, &s);
    fwrite(&len, sizeof(uint32_t), 1, file);
    fwrite(&s[0], sizeof(uint32_t), s.size(), file);
  }
}

inline void GetBinaryLibSize(const std::string &file_prefix, int64_t &total_bases, int64_t &num_reads) {
  std::ifstream lib_info_file(file_prefix + ".lib_info");
  lib_info_file >> total_bases >> num_reads;
}

#endif