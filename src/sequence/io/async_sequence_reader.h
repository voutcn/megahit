//
// Created by vout on 11/24/18.
//

#ifndef MEGAHIT_ASYNC_SEQUENCE_READER_H
#define MEGAHIT_ASYNC_SEQUENCE_READER_H

#include <future>
#include <string>
#include "sequence/io/binary_reader.h"
#include "sequence/io/contig/contig_reader.h"
#include "sequence/sequence_package.h"

template <class PackageType>
class BaseAsyncSequenceReader {
 public:
  using package_type = PackageType;
  BaseAsyncSequenceReader() = default;
  virtual ~BaseAsyncSequenceReader() = default;
  package_type &Next() {
    auto ret = next_pkg_.get();
    AsyncReadNextBatch();
    return *ret;
  }

 protected:
  void AsyncReadNextBatch() {
    auto pkg = &packages_[batch_index_];
    next_pkg_ = std::async(std::launch::async, [this, pkg]() -> package_type * {
      ReadOneBatch(pkg);
      return pkg;
    });
    batch_index_ ^= 1;
  }
  void StopReading() { next_pkg_.wait(); }
  virtual void ReadOneBatch(package_type *pkg) = 0;

 private:
  unsigned batch_index_{0};
  package_type packages_[2];
  std::future<package_type *> next_pkg_;
};

class AsyncSequenceReader : public BaseAsyncSequenceReader<SeqPackage> {
 public:
  static const int64_t kDefaultNumSeqPerBatch = 1u << 22u;
  static const int64_t kDefaultNumBasesPerBatch = 1u << 28u;

  explicit AsyncSequenceReader(
      BaseSequenceReader *reader, bool reverse = false,
      int64_t sequences_per_batch = kDefaultNumSeqPerBatch,
      int64_t bases_per_batch = kDefaultNumBasesPerBatch)
      : reader_(reader),
        reverse_(reverse),
        max_sequence_per_batch_(sequences_per_batch),
        max_bases_per_batch_(bases_per_batch) {
    AsyncReadNextBatch();
  }
  ~AsyncSequenceReader() override { StopReading(); }

 protected:
  void ReadOneBatch(SeqPackage *seq_pkg) override {
    seq_pkg->Clear();
    reader_->Read(seq_pkg, max_sequence_per_batch_, max_bases_per_batch_,
                  reverse_);
  }

 private:
  BaseSequenceReader *reader_;
  bool reverse_;
  int64_t max_sequence_per_batch_{1u << 22u};
  int64_t max_bases_per_batch_{1u << 28u};
};

class AsyncContigReader : public BaseAsyncSequenceReader<
                              std::pair<SeqPackage, std::vector<float>>> {
 public:
  explicit AsyncContigReader(const std::string &file_name)
      : reader_(file_name) {
    reader_.SetDiscardFlag(contig_flag::kLoop | contig_flag::kStandalone);
    AsyncReadNextBatch();
  }
  ~AsyncContigReader() override { StopReading(); }

 protected:
  void ReadOneBatch(package_type *pkg) override {
    pkg->first.Clear();
    pkg->second.clear();
    const int64_t kMaxNumContigs = 1u << 22u;
    const int64_t kMaxNumBases = 1u << 28u;
    const bool reverse = false;
    reader_.ReadWithMultiplicity(&pkg->first, &pkg->second, kMaxNumContigs,
                                 kMaxNumBases, reverse);
  }

 private:
  ContigReader reader_;
};

#endif  // MEGAHIT_ASYNC_SEQUENCE_READER_H
