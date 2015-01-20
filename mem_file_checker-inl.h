/*
 *  mem_file_checker-inl.h
 *  
 *  This file is a part of MEGAHIT
 *  
 *  Copyright (C) 2014 The University of Hong Kong
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

#ifndef MEM_FILE_CHECKER_INL_H__
#define MEM_FILE_CHECKER_INL_H__

#include <stdlib.h>
#include <stdio.h>
#include "helper_functions-inl.h"

inline FILE *OpenFileAndCheck(const char *filename, const char * mode) {
    FILE *fp;
    if ((fp = fopen(filename, mode)) == NULL){
        err("[ERROR] Cannot open %s. Now exit to system...\n", filename);
        exit ( -1 );
    }

    return fp;
}

inline void* MallocAndCheck(size_t size_in_byte,
                            const char *malloc_from_which_file = __FILE__,
                            int malloc_from_which_line = __LINE__) {
    void *ptr = malloc(size_in_byte);
    if (ptr == NULL && size_in_byte != 0) {
        err("[ERROR] Ran out of memory while applying %llubytes\n", (unsigned long long)size_in_byte);
        err("In file: %s, line %d\n", malloc_from_which_file, malloc_from_which_line);
        err("There may be errors as follows:\n");
        err("1) Not enough memory.\n");
        err("2) The ARRAY may be overrode.\n");
        err("3) The wild pointers.\n");
        exit(-1);
    }

    return ptr;
}

inline void* ReAllocAndCheck(void *ptr,
                            size_t size_in_byte,
                            const char *realloc_from_which_file = __FILE__,
                            int realloc_from_which_line = __LINE__) {
    void *new_ptr = realloc(ptr, size_in_byte);
    if (size_in_byte == 0 || ptr != NULL) {
        return new_ptr;
    } else {
        err("[ERROR] Ran out of memory while re-applying %llubytes\n", (unsigned long long)size_in_byte);
        err("In file: %s, line %d\n", realloc_from_which_file, realloc_from_which_line);
        err("There may be errors as follows:\n");
        err("1) Not enough memory.\n");
        err("2) The ARRAY may be overrode.\n");
        err("3) The wild pointers.\n");
        free(ptr);
        exit(-1);
    }
}

#endif // MEM_FILE_CHECKER_INL_H__