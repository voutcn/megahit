#ifndef READ_LIB_FUNCTIONS_INL_H__
#define READ_LIB_FUNCTIONS_INL_H__

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "SequencePackage.h"
#include "SequenceManager.h"
#include "mem_file_checker-inl.h"
#include "utils.h"

inline void ReadMultipleLibs(const std::string &lib_file, std::vector<SequencePackage*> &out_packages,
							 const std::string &prefix, bool is_reverse, bool is_binary) {
	std::ifstream lib_config(lib_file);
	if (!lib_config.is_open()) {
		fprintf(stderr, "File to open read_lib file: %s\n", lib_file.c_str());
		exit(1);
	}

	int lib_id = 0;
	std::string buf;
	std::string type;
	std::string file_name1;
	std::string file_name2;

	bool append_to_package = false;
	SequenceManager seq_manager;

	while (lib_config >> type) {
		if (type == "pe") {
			assert(lib_config >> file_name1 >> file_name2);
			SequencePackage *p = new SequencePackage;
			assert(p != NULL);

			seq_manager.set_read_type(SequenceManager::kPaired);

			if (is_binary) {
				seq_manager.set_file_type(SequenceManager::kBinaryReads);
				seq_manager.set_file(FormatString("%s.%d.bin", prefix.c_str(), lib_id));
			} else {
				seq_manager.set_file_type(SequenceManager::kFastxReads);
				seq_manager.set_pe_files(file_name1, file_name2);
			}

			seq_manager.ReadShortReads(1LL << 60, 1LL << 60, append_to_package, is_reverse);
			seq_manager.clear();
			out_packages.push_back(p);

		} else if (type == "se") {
			assert(lib_config >> file_name1);
			SequencePackage *p = new SequencePackage;
			assert(p != NULL);

			seq_manager.set_read_type(SequenceManager::kSingle);

			if (is_binary) {
				seq_manager.set_file_type(SequenceManager::kBinaryReads);
				seq_manager.set_file(FormatString("%s.%d.bin", prefix.c_str(), lib_id));
			} else {
				seq_manager.set_file_type(SequenceManager::kFastxReads);
				seq_manager.set_file(file_name1);
			}

			seq_manager.ReadShortReads(1LL << 60, 1LL << 60, append_to_package, is_reverse);
			seq_manager.clear();
			out_packages.push_back(p);

		} else if (type == "interleaved") {
			assert(lib_config >> file_name1);
			SequencePackage *p = new SequencePackage;
			assert(p != NULL);

			seq_manager.set_read_type(SequenceManager::kInterleaved);

			if (is_binary) {
				seq_manager.set_file_type(SequenceManager::kBinaryReads);
				seq_manager.set_file(FormatString("%s.%d.bin", prefix.c_str(), lib_id));
			} else {
				seq_manager.set_file_type(SequenceManager::kFastxReads);
				seq_manager.set_file(file_name1);
			}

			seq_manager.ReadShortReads(1LL << 60, 1LL << 60, append_to_package, is_reverse);
			seq_manager.clear();
			out_packages.push_back(p);

		} else {
			fprintf(stderr, "Cannot identify read library type %s\n", type.c_str());
			fprintf(stderr, "Valid types: pe, se, interleaved\n", type.c_str());
			exit(1);
		}

		++lib_id;
	}
}

inline void WriteMultipleLibs(std::vector<SequencePackage*> &packages, const std::string &prefix, bool is_reverse) {
	SequenceManager seq_manager;
	for (int lib_id = 0; lib_id = packages.size(); ++lib_id) {
		seq_manager.set_package(packages[i]);
		FILE *file = fopen(FormatString("%s.%d.bin", prefix.c_str(), lib_id), "wb");
		assert(file != NULL);
		seq_manager.WriteBinarySequences("wb")
		fclose(file);
	}
}

#endif