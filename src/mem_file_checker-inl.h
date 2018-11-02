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

#ifndef MEM_FILE_CHECKER_INL_H__
#define MEM_FILE_CHECKER_INL_H__

#include <stdlib.h>
#include <stdio.h>

inline FILE *OpenFileAndCheck(const char *filename, const char *mode) {
    FILE *fp;

    if ((fp = fopen(filename, mode)) == NULL) {
        fprintf(stderr, "[ERROR] Cannot open %s. Now exit to system...\n", filename);
        exit ( -1 );
    }

    return fp;
}

inline void *MallocAndCheck(size_t size_in_byte,
                            const char *malloc_from_which_file = __FILE__,
                            int malloc_from_which_line = __LINE__) {
    void *ptr = malloc(size_in_byte);

    if (ptr == NULL && size_in_byte != 0) {
        fprintf(stderr, "[ERROR] Ran out of memory while applying %llubytes\n", (unsigned long long)size_in_byte);
        fprintf(stderr, "In file: %s, line %d\n", malloc_from_which_file, malloc_from_which_line);
        fprintf(stderr, "There may be errors as follows:\n");
        fprintf(stderr, "1) Not enough memory.\n");
        fprintf(stderr, "2) The ARRAY may be overrode.\n");
        fprintf(stderr, "3) The wild pointers.\n");
        exit(-1);
    }

    return ptr;
}

inline void *ReAllocAndCheck(void *ptr,
                             size_t size_in_byte,
                             const char *realloc_from_which_file = __FILE__,
                             int realloc_from_which_line = __LINE__) {
    void *new_ptr = realloc(ptr, size_in_byte);

    if (size_in_byte != 0 && new_ptr == NULL) {
        fprintf(stderr, "[ERROR] Ran out of memory while re-applying %llubytes\n", (unsigned long long)size_in_byte);
        fprintf(stderr, "In file: %s, line %d\n", realloc_from_which_file, realloc_from_which_line);
        fprintf(stderr, "There may be errors as follows:\n");
        fprintf(stderr, "1) Not enough memory.\n");
        fprintf(stderr, "2) The ARRAY may be overrode.\n");
        fprintf(stderr, "3) The wild pointers.\n");
        free(ptr);
        exit(-1);
    }

    return new_ptr;
}

#endif // MEM_FILE_CHECKER_INL_H__