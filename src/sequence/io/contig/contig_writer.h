//
// Created by vout on 6/24/19.
//

#ifndef MEGAHIT_CONTIG_WRITER_H
#define MEGAHIT_CONTIG_WRITER_H

#include <definitions.h>
#include <atomic>
#include <string>
#include "utils/utils.h"

class ContigWriter {
 public:
  explicit ContigWriter(const std::string file_name) : file_name_(file_name) {
    file_ = xfopen(file_name.c_str(), "w");
  }

  ~ContigWriter() {
    std::FILE *info_file = xfopen((file_name_ + ".info").c_str(), "w");
    pfprintf(info_file, "{} {}\n", n_contigs_.load(), n_bases_.load());
    fclose(info_file);
    fclose(file_);
  }

  void WriteContig(const std::string &ascii_contig, unsigned k_size,
                   long long id, int flag, double multi) {
    pfprintf(file_, ">k{}_{} flag={} multi={.4} len={}\n{s}\n", k_size, id,
             flag, multi, ascii_contig.length(), ascii_contig.c_str());
    n_contigs_.fetch_add(1, std::memory_order_relaxed);
    n_bases_.fetch_add(
        ascii_contig.length() + (flag & contig_flag::kLoop) ? 28 : 0,
        std::memory_order_relaxed);
  }

  void WriteLocalContig(const std::string &ascii_contig,
                        int64_t origin_contig_id, int strand,
                        int64_t contig_id) {
    pfprintf(file_, ">lc_{}_strand_{}_id_{} flag=0 multi=1\n{s}\n",
             origin_contig_id, strand, contig_id, ascii_contig.c_str());
    n_contigs_.fetch_add(1, std::memory_order_relaxed);
    n_bases_.fetch_add(ascii_contig.length(), std::memory_order_relaxed);
  }

 private:
  std::string file_name_;
  std::FILE *file_;
  std::atomic<int64_t> n_contigs_{0};
  std::atomic<int64_t> n_bases_{0};
};

#endif  // MEGAHIT_CONTIG_WRITER_H
