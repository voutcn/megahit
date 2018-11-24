//
// Created by vout on 11/24/18.
//

#ifndef MEGAHIT_ASYNC_SEQUENCE_READER_H
#define MEGAHIT_ASYNC_SEQUENCE_READER_H

#include <future>
#include "sequence_package.h"
#include "sequence_manager.h"

template <class PackageType>
class AsyncSequenceReader {
 public:
  using package_type = PackageType;
  AsyncSequenceReader() = default;
  virtual ~AsyncSequenceReader() = default;
  package_type &Next() {
    auto &ret = next_pkg_.get();
    AsyncReadNextBatch();
    return ret;
  }
 protected:
  void AsyncReadNextBatch() {
    next_pkg_ = std::async(
        std::launch::async,
        [this]() -> package_type & {
          ReadOneBatch(&packages_[batch_index_]);
          return packages_[batch_index_];
        });
    batch_index_ ^= 1;
  }
  void StopReading() {
    next_pkg_.wait();
  }
  virtual void ReadOneBatch(package_type *pkg) = 0;
 private:
  unsigned batch_index_{0};
  package_type packages_[2];
  std::future<package_type &> next_pkg_;
};

class AsyncContigReader : public AsyncSequenceReader<std::pair<SequencePackage, std::vector<float>>> {
 public:
  AsyncContigReader(const std::string &file_name) {
    seq_manager_.set_file_type(SequenceManager::kMegahitContigs);
    seq_manager_.set_file(file_name);
    AsyncReadNextBatch();
  }
  virtual ~AsyncContigReader() {
    StopReading();
  }
 protected:
  void ReadOneBatch(package_type *pkg) {
    seq_manager_.set_package(&pkg->first);
    seq_manager_.set_float_multiplicity_vector(&pkg->second);
    int64_t kMaxNumContigs = 1 << 22;
    int64_t kMaxNumBases = 1 << 28;
    bool append = false;
    bool reverse = false;
    int discard_flag = contig_flag::kLoop | contig_flag::kStandalone;
    bool extend_loop = false;
    bool calc_depth = false;
    seq_manager_.ReadMegahitContigs(
        kMaxNumContigs, kMaxNumBases, append, reverse, discard_flag, extend_loop, calc_depth);
  }
 private:
  SequenceManager seq_manager_;
};

class AsyncReadReader : public AsyncSequenceReader<SequencePackage> {
 public:
  AsyncReadReader(const std::string &file_name) {
    seq_manager_.set_file_type(SequenceManager::kBinaryReads);
    seq_manager_.set_readlib_type(SequenceManager::kSingle); // PE info not used
    seq_manager_.set_file(file_name);
    AsyncReadNextBatch();
  }
  virtual ~AsyncReadReader() {
    StopReading();
  }
 protected:
  void ReadOneBatch(SequencePackage *seq_pkg) {
    seq_manager_.set_package(seq_pkg);
    int64_t kMaxNumReads = 1 << 22;
    int64_t kMaxNumBases = 1 << 28;
    bool append = false;
    bool reverse = false;
    seq_manager_.ReadShortReads(kMaxNumReads, kMaxNumBases, append, reverse);
  }
 private:
  SequenceManager seq_manager_;
};

#endif //MEGAHIT_ASYNC_SEQUENCE_READER_H
