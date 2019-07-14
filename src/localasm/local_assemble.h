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

#ifndef LOCAL_ASSEMBLER_H
#define LOCAL_ASSEMBLER_H

#include <string>

struct LocalAsmOption {
  std::string contig_file;
  std::string lib_file_prefix;

  uint32_t kmin{11};
  uint32_t kmax{41};
  uint32_t step{6};
  uint32_t seed_kmer{31};

  uint32_t min_contig_len{200};
  uint32_t sparsity{8};
  double similarity{0.8};
  uint32_t min_mapping_len{75};

  uint32_t num_threads{0};
  std::string output_file;
};

void RunLocalAssembly(const LocalAsmOption &opt);

#endif