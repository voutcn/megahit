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

#ifndef READ_LIB_FUNCTIONS_INL_H__
#define READ_LIB_FUNCTIONS_INL_H__

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <memory>
#include <sequence/readers/binary_reader.h>

#include "utils/utils.h"
#include "lib_info.h"
#include "sequence/sequence_package.h"
#include "sequence/readers/fastx_reader.h"
#include "sequence/readers/pair_end_fastx_reader.h"
#include "utils/safe_alloc_open-inl.h"


inline void WriteBinarySequences(const SeqPackage & pkg, FILE *file, int64_t from = 0, int64_t to = -1) {
  if (to == -1) {
    to = pkg.Size() - 1;
  }

  uint32_t len;
  std::vector<uint32_t> s;

  for (int64_t i = from; i <= to; ++i) {
    len = pkg.SequenceLength(i);
    pkg.FetchSequence(i, &s);
    fwrite(&len, sizeof(uint32_t), 1, file);
    fwrite(&s[0], sizeof(uint32_t), s.size(), file);
  }
}

inline void ReadAndWriteMultipleLibs(const std::string &lib_file,
                                     const std::string &out_prefix) {
    std::ifstream lib_config(lib_file);

    if (!lib_config.is_open()) {
        xfatal("File to open read_lib file: %s\n", lib_file.c_str());
    }

    FILE *bin_file = xfopen(FormatString("%s.bin", out_prefix.c_str()), "wb");

    SeqPackage seq_batch;
    std::vector<lib_info_t> lib_info;

    std::string metadata;
    std::string type;
    std::string file_name1;
    std::string file_name2;

    int64_t total_reads = 0;
    int64_t total_bases = 0;
  seq_batch.Clear();
    lib_info.clear();

    while (std::getline(lib_config, metadata)) {
        lib_config >> type;
        std::unique_ptr<BaseSequenceReader> reader;
        if (type == "pe") {
            lib_config >> file_name1 >> file_name2;
            reader.reset(new PairEndFastxReader(file_name1, file_name2));
        } else if (type == "se") {
            lib_config >> file_name1;
            reader.reset(new FastxReader(file_name1));
        } else if (type == "interleaved") {
            lib_config >> file_name1;
            reader.reset(new FastxReader(file_name1));
        } else {
            xerr("Cannot identify read library type %s\n", type.c_str());
            xfatal("Valid types: pe, se, interleaved\n");
        }

        int64_t start = total_reads;
        int64_t num_read = 0;
        int reads_per_batch = 1u << 22;
        int bases_per_batch = 1u << 28;
        int max_read_len = 0;

        while (true) {
            num_read = reader->Read(&seq_batch, reads_per_batch, bases_per_batch, false, true);

            if (num_read == 0) {
                break;
            }

            total_reads += num_read;
            total_bases += seq_batch.BaseCount();
            WriteBinarySequences(seq_batch, bin_file);
            max_read_len = std::max(max_read_len, (int) seq_batch.MaxSequenceLength());
          seq_batch.Clear();
        }

        if (type == "pe" && (total_reads - start) % 2 != 0) {
            xerr("PE library number of reads is odd: %lld!\n", total_reads - start);
            xfatal("File(s): %s\n", metadata.c_str());
        }

        if (type == "interleaved" && (total_reads - start) % 2 != 0) {
            xerr("PE library number of reads is odd: %lld!\n", total_reads - start);
            xfatal("File(s): %s\n", metadata.c_str());
        }

        xinfo("Lib %d (%s): %s, %lld reads, %d max length\n",
            lib_info.size(), metadata.c_str(), type.c_str(), total_reads - start, max_read_len);

        lib_info.emplace_back(&seq_batch, start, total_reads - 1, max_read_len, type != "se", metadata);
        std::getline(lib_config, metadata); // eliminate the "\n"
    }

    FILE *lib_info_file = xfopen(FormatString("%s.lib_info", out_prefix.c_str()), "w");
    fprintf(lib_info_file, "%zu %zu\n", total_bases, total_reads);

    for (auto &i : lib_info) {
        fprintf(lib_info_file, "%s\n", i.metadata.c_str());
        fprintf(lib_info_file, "%" PRId64 " %" PRId64 " %d %s\n", i.from, i.to,
                i.max_read_len, i.is_pe ? "pe" : "se");
    }

    fclose(lib_info_file);
}

inline void GetBinaryLibSize(const std::string &file_prefix, int64_t &total_bases, int64_t &num_reads) {
    std::ifstream lib_info_file(file_prefix + ".lib_info");
    lib_info_file >> total_bases >> num_reads;
}

inline void ReadBinaryLibs(const std::string &file_prefix, SeqPackage &package, std::vector<lib_info_t> &lib_info,
                           bool is_reverse = false) {
    std::ifstream lib_info_file(file_prefix + ".lib_info");
    int64_t start, end;
    int max_read_len;
    int64_t total_bases, num_reads;
    std::string metadata, pe_or_se;

    lib_info_file >> total_bases >> num_reads;
    std::getline(lib_info_file, metadata); // eliminate the "\n"

    while (std::getline(lib_info_file, metadata)) {
        lib_info_file >> start >> end >> max_read_len >> pe_or_se;
        lib_info.emplace_back(&package, start, end, max_read_len, pe_or_se == "pe", metadata);
        std::getline(lib_info_file, metadata); // eliminate the "\n"
    }

    xinfo("Before reserve for %lu reads, %lu bases, sizeof seq_package: %lu\n", num_reads, total_bases,
          package.SizeInByte());

    package.ReserveSequences(num_reads);
    package.ReserveBases(total_bases);
    BinaryReader reader(file_prefix + ".bin");
    xinfo("Before reading, sizeof seq_package: %lld\n", package.SizeInByte());
  package.Clear();
    reader.ReadAll(&package, is_reverse);
    xinfo("After reading, sizeof seq_package: %lld\n", package.SizeInByte());
}

#endif