/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong & L3 Bioinformatics
 * Limited
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

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "definitions.h"
#include "utils/cpu_dispatch.h"
#include "utils/utils.h"

int main_assemble(int argc, char **argv);
int main_local(int argc, char **argv);
int main_iterate(int argc, char **argv);
int main_build_lib(int argc, char **argv);

int main_kmer_count(int argc, char **argv);
int main_read2sdbg(int argc, char **argv);
int main_seq2sdbg(int argc, char **argv);

int main_contig2fastg(int argc, char **argv);
int main_read_stat(int argc, char **argv);
int main_filter_by_len(int argc, char **argv);

void show_help(const char *program_name) {
  pfprintf(
      stderr,
      "Usage: {s} <sub_program> [sub options]\n"
      "    sub-programs:\n"
      "       assemble       assemble from SdBG\n"
      "       local          local asssembly\n"
      "       iterate        extract iterative edges\n"
      "       buildlib       build read library\n"
      "       count          kmer counting\n"
      "       read2sdbg      build sdbg from reads\n"
      "       seq2sdbg       build sdbg from megahit contigs + edges\n"
      "       contig2fastg   convert MEGAHIT's k*.contigs.fa to fastg format\n"
      "       readstat       calculate read stats (# of reads, bases, longest, "
      "shortest, average)\n"
      "       filterbylen    filter contigs by length\n"
      "       checkcpu       check whether the run-time CPU supports POPCNT "
      "and BMI2\n"
      "       checkpopcnt    check whether the run-time CPU supports POPCNT\n"
      "       checkbmi2      check whether the run-time CPU supports BMI2\n"
      "       dumpversion    dump version\n"
      "       kmax           the largest k value supported\n",
      program_name);
}

int main(int argc, char **argv) {
  if (argc < 2) {
    show_help(argv[0]);
    exit(1);
  }

  if (strcmp(argv[1], "assemble") == 0) {
    return main_assemble(argc - 1, argv + 1);
  } else if (strcmp(argv[1], "local") == 0) {
    return main_local(argc - 1, argv + 1);
  } else if (strcmp(argv[1], "iterate") == 0) {
    return main_iterate(argc - 1, argv + 1);
  } else if (strcmp(argv[1], "buildlib") == 0) {
    return main_build_lib(argc - 1, argv + 1);
  } else if (strcmp(argv[1], "count") == 0) {
    return main_kmer_count(argc - 1, argv + 1);
  } else if (strcmp(argv[1], "read2sdbg") == 0) {
    return main_read2sdbg(argc - 1, argv + 1);
  } else if (strcmp(argv[1], "seq2sdbg") == 0) {
    return main_seq2sdbg(argc - 1, argv + 1);
  } else if (strcmp(argv[1], "contig2fastg") == 0) {
    return main_contig2fastg(argc - 1, argv + 1);
  } else if (strcmp(argv[1], "readstat") == 0) {
    return main_read_stat(argc - 1, argv + 1);
  } else if (strcmp(argv[1], "filterbylen") == 0) {
    return main_filter_by_len(argc - 1, argv + 1);
  } else if (strcmp(argv[1], "checkcpu") == 0) {
    pprintf("{}\n", HasPopcnt() && HasBmi2());
  } else if (strcmp(argv[1], "checkpopcnt") == 0) {
    pprintf("{}\n", HasPopcnt());
  } else if (strcmp(argv[1], "checkbmi2") == 0) {
    pprintf("{}\n", HasBmi2());
  } else if (strcmp(argv[1], "dumpversion") == 0) {
    pprintf("{s}\n", PACKAGE_VERSION);
  } else if (strcmp(argv[1], "kmax") == 0) {
    pprintf("{}\n", kMaxK);
  } else {
    show_help(argv[0]);
    return 1;
  }

  return 0;
}
