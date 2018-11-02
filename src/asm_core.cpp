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

#include "utils.h"
#include "definitions.h"

int main_assemble(int argc, char **argv);
int main_local(int argc, char **argv);
int main_iterate(int argc, char **argv);
int main_build_lib(int argc, char **argv);

void show_help(const char *program_name) {
    fprintf(stderr, "Usage: %s <sub_program> [sub options]\n"
            "    sub-programs:\n"
            "       assemble              assemble from SdBG\n"
            "       local                 local asssembly\n"
            "       iterate               extract iterative edges\n"
            "       buildlib              build read library\n"
            "       dumpversion           dump version\n"
            "       kmax                  the largest k value supported\n",
            program_name);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        show_help(argv[0]);
        exit(1);
    }

    if (strcmp(argv[1], "assemble") == 0) {
        return main_assemble(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "local") == 0) {
        return main_local(argc - 1 , argv + 1);
    }
    else if (strcmp(argv[1], "iterate") == 0) {
        return main_iterate(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "buildlib") == 0) {
        return main_build_lib(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "dumpversion") == 0) {
        printf("%s\n", PACKAGE_VERSION);
        return 0;
    }
    else if (strcmp(argv[1], "kmax") == 0) {
        printf("%d\n", kMaxK);
    }
    else {
        show_help(argv[0]);
        exit(1);
    }

    return 0;
}