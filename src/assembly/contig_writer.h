//
// Created by vout on 11/21/18.
//

#ifndef MEGAHIT_CONTIG_WRITER_H
#define MEGAHIT_CONTIG_WRITER_H

#include <string>
#include <cstdint>
#include "histgram.h"

class UnitigGraph;

// TODO refactor
inline void WriteContig(const std::string &label, unsigned k_size,
                        long long id, int flag, double multi, FILE *file) {
  fprintf(file, ">k%d_%lld flag=%d multi=%.4lf len=%lu\n%s\n",
          k_size, id, flag, multi, label.length(), label.c_str());
}

void OutputContigs(UnitigGraph &graph, FILE *contig_file, FILE *final_file,
                   bool change_only, uint32_t min_standalone);

#endif //MEGAHIT_CONTIG_WRITER_H
