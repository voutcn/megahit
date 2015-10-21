/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong & L3 Bioinformatics Limited
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "definitions.h"

int main_contig2fastg(int argc, char **argv);
int main_read_stat(int argc, char **argv);
int main_trim_lowq_tail(int argc, char **argv);
int main_filter_by_len(int argc, char **argv);
int main_extract_pe(int argc, char **argv);

void show_help(const char *program_name) {
    fprintf(stderr, "Usage: %s <sub_program> [sub options]\n"
            "    sub-programs:\n"
            "       contig2fastg          convert MEGAHIT's k*.contigs.fa to fastg format that can be viewed by Bandage\n"
            "       readstat              calculate read stats (# of reads, bases, longest, shortest, average)\n"
            "       trim                  trim low quality tail of fastq reads\n"
            "       filterbylen           filter contigs by length\n"
            "       extractpe             extract pe reads and se reads from fasta/fastq files\n"
            "       dumpversion           dump version\n",
            program_name);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        show_help(argv[0]);
        exit(1);
    }

    if (strcmp(argv[1], "contig2fastg") == 0) {
        return main_contig2fastg(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "readstat") == 0) {
        return main_read_stat(argc - 1 , argv + 1);
    }
    else if (strcmp(argv[1], "trim") == 0) {
        return main_trim_lowq_tail(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "filterbylen") == 0) {
        return main_filter_by_len(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "dumpversion") == 0) {
        printf("%s\n", PACKAGE_VERSION);
        return 0;
    }
    else if (strcmp(argv[1], "extractpe") == 0) {
        return main_extract_pe(argc - 1, argv + 1);
    }
    else {
        show_help(argv[0]);
        return 1;
    }
}