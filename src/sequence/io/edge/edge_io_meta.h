//
// Created by vout on 6/24/19.
//

#ifndef MEGAHIT_EDGE_IO_META_H
#define MEGAHIT_EDGE_IO_META_H

/**
 * A bucket is a length-L prefix of an edge
 */
struct EdgeIoBucketInfo {
  int file_id{-1};
  int64_t file_offset{0};
  int64_t total_number{0};
};

struct EdgeIoMetadata {
  uint32_t kmer_size{};
  uint32_t words_per_edge{};
  uint32_t num_files{};
  int64_t num_edges{};
  bool is_sorted{true};
  std::vector<EdgeIoBucketInfo> buckets;

  void Serialize(std::ofstream &os) {
    if (is_sorted) {
      num_edges = 0;
      for (const auto &b : buckets) {
        num_edges += b.total_number;
      }
    }

    os << "kmer_size " << kmer_size << '\n'
       << "words_per_edge " << words_per_edge << '\n'
       << "num_files " << num_files << '\n'
       << "num_buckets " << buckets.size() << '\n'
       << "num_edges " << num_edges << '\n'
       << "is_sorted " << is_sorted << '\n';

    for (unsigned i = 0; i < buckets.size(); ++i) {
      os << i << ' ' << buckets[i].file_id << ' ' << buckets[i].file_offset
         << ' ' << buckets[i].total_number << '\n';
    }
  }

  void Deserialize(std::ifstream &is) {
    unsigned num_buckets;
    ScanField(is, "kmer_size", kmer_size);
    ScanField(is, "words_per_edge", words_per_edge);
    ScanField(is, "num_files", num_files);
    ScanField(is, "num_buckets", num_buckets);
    ScanField(is, "num_edges", num_edges);
    ScanField(is, "is_sorted", is_sorted);

    buckets.resize(num_buckets);

    for (unsigned i = 0; i < num_buckets; ++i) {
      unsigned b_id;
      is >> b_id >> buckets[i].file_id >> buckets[i].file_offset >>
          buckets[i].total_number;
      if (b_id != i) {
        xfatal("Invalid format: bucket id not matched!\n");
      }
      if (buckets[i].file_id >= static_cast<int>(num_files)) {
        xfatal("Record ID {} is greater than number of files {}\n",
               buckets[i].file_id, num_files);
      }
    }
  }
};

#endif  // MEGAHIT_EDGE_IO_META_H
