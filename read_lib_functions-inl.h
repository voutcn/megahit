/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong
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
#include <assert.h>
#include <stdint.h>
#include <inttypes.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "utils.h"
#include "lib_info.h"
#include "sequence_manager.h"
#include "sequence_package.h"
#include "mem_file_checker-inl.h"

inline void ReadMultipleLibs(const std::string &lib_file, SequencePackage &package,
                             std::vector<lib_info_t> &lib_info, bool is_reverse) {
    std::ifstream lib_config(lib_file);
    if (!lib_config.is_open()) {
        xerr_and_exit("File to open read_lib file: %s\n", lib_file.c_str());
    }

    std::string metadata;
    std::string type;
    std::string file_name1;
    std::string file_name2;

    bool append_to_package = true;
    bool trimN = true;
    SequenceManager seq_manager(&package);
    package.clear();
    lib_info.clear();

    while (std::getline(lib_config, metadata)) {
        assert(lib_config >> type);
        if (type == "pe") {
            assert(lib_config >> file_name1 >> file_name2);

            seq_manager.set_readlib_type(SequenceManager::kPaired);
            seq_manager.set_file_type(SequenceManager::kFastxReads);
            seq_manager.set_pe_files(file_name1, file_name2);

        } else if (type == "se") {
            assert(lib_config >> file_name1);

            seq_manager.set_readlib_type(SequenceManager::kSingle);
            seq_manager.set_file_type(SequenceManager::kFastxReads);
            seq_manager.set_file(file_name1);

        } else if (type == "interleaved") {
            assert(lib_config >> file_name1);

            seq_manager.set_readlib_type(SequenceManager::kInterleaved);
            seq_manager.set_file_type(SequenceManager::kFastxReads);
            seq_manager.set_file(file_name1);

        } else {
            xerr("Cannot identify read library type %s\n", type.c_str());
            xerr_and_exit("Valid types: pe, se, interleaved\n");
        }

        int64_t start = package.size();
        seq_manager.ReadShortReads(1LL << 60, 1LL << 60, append_to_package, is_reverse, trimN);
        seq_manager.clear();
        int64_t end = package.size() - 1;

        int max_read_len = 0;
        for (int64_t i = start; i <= end; ++i) {
            if (max_read_len < (int)package.length(i)) {
                max_read_len = package.length(i);
            }
        }

        if (type == "pe" && (end - start + 1) % 2 != 0) {
            xerr("PE library number of reads is odd: %lld!\n", end - start + 1);
            xerr_and_exit("File(s): %s\n", metadata.c_str())
        }

        if (type == "interleaved" && (end - start + 1) % 2 != 0) {
            xerr("PE library number of reads is odd: %lld!\n", end - start + 1);
            xerr_and_exit("File(s): %s\n", metadata.c_str())
        }

        lib_info.push_back(lib_info_t(&package, start, end, max_read_len, type != "se", metadata));
        std::getline(lib_config, metadata); // eliminate the "\n"
    }
}

inline void ReadBinaryLibs(const std::string &file_prefix, SequencePackage &package, std::vector<lib_info_t> &lib_info) {
    package.clear();
    lib_info.clear();

    std::ifstream lib_info_file(file_prefix + ".lib_info");
    int64_t start, end;
    int max_read_len;
    int64_t total_bases, num_reads;
    std::string metadata, pe_or_se;

    assert(lib_info_file >> total_bases >> num_reads);
    std::getline(lib_info_file, metadata); // eliminate the "\n"

    while (std::getline(lib_info_file, metadata)) {
        assert(lib_info_file >> start >> end >> max_read_len >> pe_or_se);
        lib_info.push_back(lib_info_t(&package, start, end, max_read_len, pe_or_se == "pe", metadata));
        std::getline(lib_info_file, metadata); // eliminate the "\n"
    }

    package.reserve_num_seq(num_reads);
    package.reserve_bases(total_bases);
    SequenceManager seq_manager(&package);
    seq_manager.set_file_type(SequenceManager::kBinaryReads);
    seq_manager.set_file(FormatString("%s.bin", file_prefix.c_str()));

    bool append_to_package = false;
    bool is_reverse = false;
    seq_manager.ReadShortReads(1LL << 60, 1LL << 60, append_to_package, is_reverse);

    // for (int i = 0; i < package.length(0); ++i) {
    // 	xerr_and_exit("%c", "ACGT"[package.get_base(0, i)]);
    // }
    // xerr_and_exit("\n");

    // for (int i = 0; i < package.length(1); ++i) {
    // 	xerr_and_exit("%c", "ACGT"[package.get_base(1, i)]);
    // }
    // xerr_and_exit("\n");
}

inline void WriteMultipleLibs(SequencePackage &package, std::vector<lib_info_t> &lib_info, const std::string &file_prefix, bool is_reverse) {
    SequenceManager seq_manager(&package);
    FILE *bin_file = OpenFileAndCheck(FormatString("%s.bin", file_prefix.c_str()), "wb");
    seq_manager.WriteBinarySequences(bin_file, is_reverse);
    fclose(bin_file);

    FILE *lib_info_file = OpenFileAndCheck(FormatString("%s.lib_info", file_prefix.c_str()), "w");
    fprintf(lib_info_file, "%zu %zu\n", package.base_size(), package.size());
    for (unsigned i = 0; i < lib_info.size(); ++i) {
        fprintf(lib_info_file, "%s\n", lib_info[i].metadata.c_str());
        fprintf(lib_info_file, "%" PRId64 " %" PRId64 " %d %s\n", lib_info[i].from, lib_info[i].to,
                lib_info[i].max_read_len, lib_info[i].is_pe ? "pe" : "se");
    }
    fclose(lib_info_file);
}

#endif