//
// Created by vout on 6/29/19.
//

#include "sequence_lib.h"

void SequenceLibCollection::Build(const std::string &lib_file, const std::string &out_prefix) {
  std::ifstream lib_config(lib_file);

  if (!lib_config.is_open()) {
    xfatal("File to open read_lib file: {}\n", lib_file.c_str());
  }

  FILE *bin_file = xfopen((out_prefix + ".bin").c_str(), "wb");

  SeqPackage seq_batch;
  std::string metadata;
  std::string type;
  std::string file_name1;
  std::string file_name2;

  int64_t total_reads = 0;
  int64_t total_bases = 0;
  seq_batch.Clear();

  std::vector<SequenceLib> libs;

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
      xerr("Cannot identify read library type {}\n", type.c_str());
      xfatal("Valid types: pe, se, interleaved\n");
    }

    int64_t start = total_reads;
    int64_t num_read = 0;
    const int reads_per_batch = 1u << 22;
    const int bases_per_batch = 1u << 28;
    int max_read_len = 0;

    while (true) {
      num_read = reader->Read(&seq_batch, reads_per_batch, bases_per_batch);

      if (num_read == 0) {
        break;
      }

      total_reads += num_read;
      total_bases += seq_batch.base_count();
      WriteBinarySequences(seq_batch, bin_file);
      max_read_len = std::max(max_read_len, (int) seq_batch.max_length());
      seq_batch.Clear();
    }

    if (type == "pe" && (total_reads - start) % 2 != 0) {
      xerr("PE library number of reads is odd: {}!\n", total_reads - start);
      xfatal("File(s): {}\n", metadata.c_str());
    }

    if (type == "interleaved" && (total_reads - start) % 2 != 0) {
      xerr("PE library number of reads is odd: {}!\n", total_reads - start);
      xfatal("File(s): {}\n", metadata.c_str());
    }

    xinfo("Lib {} ({s}): {s}, {} reads, {} max length\n", libs.size(), metadata.c_str(), type.c_str(),
          total_reads - start, max_read_len);

    libs.emplace_back(nullptr, start, total_reads, max_read_len, type != "se", metadata);
    std::getline(lib_config, metadata);  // eliminate the "\n"
  }

  std::ofstream lib_info_file(out_prefix + ".lib_info");
  lib_info_file << total_bases << ' ' << total_reads << '\n';

  for (auto &lib : libs) {
    lib.DumpMetadata(lib_info_file);
  }
  lib_info_file.close();
}

void SequenceLibCollection::Read(const std::string &file_prefix, SeqPackage *pkg, bool reverse_seq) {
  std::ifstream lib_info_file(file_prefix + ".lib_info");
  int64_t total_bases, num_reads;
  bool is_paired;
  std::string metadata;

  lib_info_file >> total_bases >> num_reads;
  std::getline(lib_info_file, metadata);  // eliminate the "\n"

  while (std::getline(lib_info_file, metadata)) {
    int64_t start, end;
    int max_read_len;
    lib_info_file >> start >> end >> max_read_len >> is_paired;
    libs_.emplace_back(pkg, start, end, max_read_len, is_paired, metadata);
    std::getline(lib_info_file, metadata);  // eliminate the "\n"
  }

  pkg->Clear();
  pkg->ReserveSequences(num_reads);
  pkg->ReserveBases(total_bases);
  BinaryReader reader(file_prefix + ".bin");

  xinfo("Before reading, sizeof seq_package: {}\n", pkg->size_in_byte());
  reader.ReadAll(pkg, reverse_seq);
  xinfo("After reading, sizeof seq_package: {}\n", pkg->size_in_byte());
}
