#ifndef READ_LIB_FUNCTIONS_INL_H__
#define READ_LIB_FUNCTIONS_INL_H__

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "lib_info.h"
#include "sequence_manager.h"
#include "sequence_package.h"
#include "mem_file_checker-inl.h"

inline void ReadMultipleLibs(const std::string &lib_file, SequencePackage &package,
	                         std::vector<lib_info_t> &lib_info, bool is_reverse) {
	std::ifstream lib_config(lib_file);
	if (!lib_config.is_open()) {
		fprintf(stderr, "File to open read_lib file: %s\n", lib_file.c_str());
		exit(1);
	}

	std::string buf;
	std::string type;
	std::string file_name1;
	std::string file_name2;

	bool append_to_package = false;
	SequenceManager seq_manager(&package);
	package.clear();
	lib_info.clear();

	while (lib_config >> type) {
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
			fprintf(stderr, "Cannot identify read library type %s\n", type.c_str());
			fprintf(stderr, "Valid types: pe, se, interleaved\n");
			exit(1);
		}

		int64_t start = package.size();
		seq_manager.ReadShortReads(1LL << 60, 1LL << 60, append_to_package, is_reverse);
		seq_manager.clear();
		int64_t end = package.size() - 1;

		int max_read_len = 0;
		for (int64_t i = start; i <= end; ++i) {
			if (max_read_len < (int)package.length(i)) {
				max_read_len = package.length(i);
			}
		}
		lib_info.push_back(lib_info_t(&package, start, end, max_read_len, type != "se"));
	}
}

inline void ReadBinaryLibs(const std::string &file_prefix, SequencePackage &package, std::vector<lib_info_t> &lib_info) {
	package.clear();
	lib_info.clear();

	FILE *lib_info_file = OpenFileAndCheck(FormatString("%s.lib_info", file_prefix.c_str()), "r");
	int64_t start, end;
	int max_read_len;
	while (fscanf(lib_info_file, "%lu", &start) == 1) {
		char is_pe[32];
		assert(fscanf(lib_info_file, "%lu", &end) == 1);
		assert(fscanf(lib_info_file, "%d", &max_read_len) == 1);
		assert(fscanf(lib_info_file, " %s", is_pe) == 1);
		lib_info.push_back(lib_info_t(&package, start, end, max_read_len, strcmp(is_pe, "pe") == 0));
	}
	fclose(lib_info_file);

	SequenceManager seq_manager(&package);
	seq_manager.set_file_type(SequenceManager::kBinaryReads);
	seq_manager.set_file(FormatString("%s.bin", file_prefix.c_str()));

	bool append_to_package = false;
	bool is_reverse = false;
	seq_manager.ReadShortReads(1LL << 60, 1LL << 60, append_to_package, is_reverse);
}

inline void WriteMultipleLibs(SequencePackage &package, std::vector<lib_info_t> &lib_info, const std::string &file_prefix, bool is_reverse) {
	SequenceManager seq_manager(&package);
	FILE *bin_file = OpenFileAndCheck(FormatString("%s.bin", file_prefix.c_str()), "wb");
	seq_manager.WriteBinarySequences(bin_file, is_reverse);
	fclose(bin_file);

	FILE *se_file = OpenFileAndCheck(FormatString("%s.lib_info", file_prefix.c_str()), "w");
	for (unsigned i = 0; i < lib_info.size(); ++i) {
		fprintf(se_file, "%ld %ld %d %s\n", lib_info[i].from, lib_info[i].to, 
			                                lib_info[i].max_read_len, lib_info[i].is_pe ? "pe" : "se");
	}
	fclose(se_file);
}

#endif