//
// Created by vout on 11/21/18.
//

#ifndef MEGAHIT_CONTIG_WRITER_H
#define MEGAHIT_CONTIG_WRITER_H

#include <string>
#include <cstdint>
#include <cassert>
#include "histgram.h"

class UnitigGraph;

static inline char Complement(char c) {
  if (c >= 0 && c < 4) {
    return 3 - c;
  }
  switch (c) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    default: assert(false);
  }
  return 0;
}

inline void ReverseComplement(std::string &s) {
  int i, j;
  for (i = 0, j = s.length() - 1; i < j; ++i, --j) {
    std::swap(s[i], s[j]);
    s[i] = Complement(s[i]);
    s[j] = Complement(s[j]);
  }
  if (i == j) {
    s[i] = Complement(s[i]);
  }
}

// TODO refactor
inline void WriteContig(const std::string &ascii_contig, unsigned k_size,
                        long long id, int flag, double multi, FILE *file) {

  std::string rev_comp = ascii_contig;
  ReverseComplement(rev_comp);

  fprintf(file, ">k%d_%lld flag=%d multi=%.4lf len=%lu\n%s\n",
          k_size, id, flag, multi, ascii_contig.length(),
          ascii_contig > rev_comp ? ascii_contig.c_str() : rev_comp.c_str());
}

void OutputContigs(UnitigGraph &graph, FILE *contig_file, FILE *final_file,
                   bool change_only, uint32_t min_standalone);

#endif //MEGAHIT_CONTIG_WRITER_H
