#-------------------------------------------------------------------------------
# MEGAHIT
# Copyright (C) 2014 The University of Hong Kong
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#
# Makefile usage
#
# make <target> [use_gpu=<0|1>] [sm=<XXX,...>] [abi=<0|1>] [open64=<0|1>] [verbose=<0|1>] [keep=<0|1>]
#
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# use GPU or not
#-------------------------------------------------------------------------------

ifeq ($(use_gpu), 1)
	NVCC = "$(shell which nvcc)"
	NVCC_VERSION = $(strip $(shell nvcc --version | grep release | sed 's/.*release //' |  sed 's/,.*//'))
endif 

# detect OS
OSUPPER = $(shell uname -s 2>/dev/null | tr [:lower:] [:upper:])

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------

CUDA_INC = "$(shell dirname $(NVCC))/../include"
INC = -I$(CUDA_INC) -I.

#-------------------------------------------------------------------------------
# SM Arch
#-------------------------------------------------------------------------------
COMMA = ,
ifdef sm
	SM_ARCH = $(subst $(COMMA),-,$(sm))
else 
    SM_ARCH = 300
endif

# Only one arch per tuning binary
ifeq (350, $(findstring 350, $(SM_ARCH)))
    SM_TARGETS = -arch=sm_35
    SM_ARCH = 350
endif
ifeq (300, $(findstring 300, $(SM_ARCH)))
    SM_TARGETS = -arch=sm_30
    SM_ARCH = 300
endif
ifeq (200, $(findstring 200, $(SM_ARCH)))
    SM_TARGETS = -arch=sm_20
    SM_ARCH = 200
endif
ifeq (130, $(findstring 130, $(SM_ARCH)))
    SM_TARGETS = -arch=sm_13
    SM_ARCH = 130
endif
ifeq (110, $(findstring 110, $(SM_ARCH)))
    SM_TARGETS = -arch=sm_11 
    SM_ARCH = 110
endif
ifeq (100, $(findstring 100, $(SM_ARCH)))
    SM_TARGETS = -arch=sm_10 
    SM_ARCH = 100
endif

#-------------------------------------------------------------------------------
# Compiler Flags
#-------------------------------------------------------------------------------

NVCCFLAGS = -Xptxas -v -Xcudafe -\# -cuda --ptxas-options=-v

# force 64bits
CPU_ARCH = -m64
CPU_ARCH_SUFFIX = x86_64

ifeq (,$(findstring 3.0, $(NVCC_VERSION)))
ifneq ($(abi), 1)
# Disable the ABI by default for 3.1+
	NVCCFLAGS += -Xptxas -abi=no
endif
endif

# CUDA ABI enable/disable (enabled by default) 
ifneq ($(abi), 0)
    ABI_SUFFIX = abi
else 
    NVCCFLAGS += -Xptxas -abi=no
    ABI_SUFFIX = noabi
endif

# NVVM/Open64 middle-end compiler (nvvm by default)
ifeq ($(open64), 1)
    NVCCFLAGS += -open64
    PTX_SUFFIX = open64
else 
    PTX_SUFFIX = nvvm
endif

# Verbose toolchain output from nvcc
ifeq ($(verbose), 1)
    NVCCFLAGS += -v
endif

# Keep intermediate compilation artifacts
ifeq ($(keep), 1)
    NVCCFLAGS += -keep
endif

# Data type size to compile a schmoo binary for
ifdef tunesize
    TUNE_SIZE = $(tunesize)
else 
    TUNE_SIZE = 4
endif

SUFFIX = $(TUNE_SIZE)B_sm$(SM_ARCH)_$(PTX_SUFFIX)_$(NVCC_VERSION)_$(ABI_SUFFIX)_$(CPU_ARCH_SUFFIX)

#-------------------------------------------------------------------------------
# Dependency Lists
#-------------------------------------------------------------------------------

rwildcard=$(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2) $(filter $(subst *,%,$2),$d))

DEPS =   ./Makefile \
        $(call rwildcard,cub/,*.cuh)

#-------------------------------------------------------------------------------
# g++ and its options
#-------------------------------------------------------------------------------
# CC = /nas1/dhli/gcc/4.8.3/rtf/bin/g++
CC = g++
CUDALIBFLAG = -L/usr/local/cuda/lib64/ -lcuda -lcudart
CFLAGS = -O3 -Wall -funroll-loops -march=native -fomit-frame-pointer -maccumulate-outgoing-args -fprefetch-loop-arrays -lm -static-libgcc -mpopcnt -fopenmp -g -std=c++0x
DEPS = Makefile
BIN_DIR = ./bin/

#-------------------------------------------------------------------------------
# CPU & GPU version
#-------------------------------------------------------------------------------

ifeq ($(use_gpu), 1)
all: make_bin_dir $(BIN_DIR)assembler iterate_edges_all $(BIN_DIR)sdbg_builder_cuda_$(SUFFIX) $(BIN_DIR)sdbg_builder_cpu
else
all: make_bin_dir $(BIN_DIR)assembler iterate_edges_all $(BIN_DIR)sdbg_builder_cpu  
endif

#-------------------------------------------------------------------------------
# CPU objectives
#-------------------------------------------------------------------------------
%.o: %.cpp %.h $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

.cx1_functions_cpu.o: cx1_functions.cpp $(DEPS)
	$(CC) $(CFLAGS) -c cx1_functions.cpp -o .cx1_functions_cpu.o -D DISABLE_GPU

#-------------------------------------------------------------------------------
# CPU Applications
#-------------------------------------------------------------------------------
$(BIN_DIR)sdbg_builder_cpu: sdbg_builder.cpp .cx1_functions_cpu.o lv2_cpu_sort.h $(DEPS)
	$(CC) $(CFLAGS) sdbg_builder.cpp .cx1_functions_cpu.o options_description.o -o $(BIN_DIR)sdbg_builder_cpu -D DISABLE_GPU -lz

$(BIN_DIR)assembler: assembler.cpp succinct_dbg.o rank_and_select.o assembly_algorithms.o branch_group.o options_description.o unitig_graph.o $(DEPS)
	$(CC) $(CFLAGS) assembler.cpp rank_and_select.o succinct_dbg.o assembly_algorithms.o branch_group.o options_description.o unitig_graph.o -o $(BIN_DIR)assembler

iterate_edges_all: $(BIN_DIR)iterate_edges_k61 $(BIN_DIR)iterate_edges_k92 $(BIN_DIR)iterate_edges_k124

$(BIN_DIR)iterate_edges_k61: iterate_edges.cpp iterate_edges.h options_description.o $(DEPS)
	$(CC) $(CFLAGS) -lz iterate_edges.cpp options_description.o -o $(BIN_DIR)iterate_edges_k61 -D KMER_NUM_UINT64=2

$(BIN_DIR)iterate_edges_k92: iterate_edges.cpp iterate_edges.h options_description.o $(DEPS)
	$(CC) $(CFLAGS) -lz iterate_edges.cpp options_description.o -o $(BIN_DIR)iterate_edges_k92 -D KMER_NUM_UINT64=3

$(BIN_DIR)iterate_edges_k124: iterate_edges.cpp iterate_edges.h options_description.o $(DEPS)
	$(CC) $(CFLAGS) -lz iterate_edges.cpp options_description.o -o $(BIN_DIR)iterate_edges_k124 -D KMER_NUM_UINT64=4

$(BIN_DIR)query_sdbg: query_sdbg.cpp succinct_dbg.o rank_and_select.o assembly_algorithms.o branch_group.o unitig_graph.o $(DEPS)
	$(CC) $(CFLAGS) query_sdbg.cpp rank_and_select.o succinct_dbg.o assembly_algorithms.o branch_group.o unitig_graph.o -o $(BIN_DIR)query_sdbg

$(BIN_DIR)rank_and_select_sample: rank_and_select_sample.cpp rank_and_select.o $(DEPS)
	$(CC) $(CFLAGS) rank_and_select.o rank_and_select_sample.cpp -o $(BIN_DIR)rank_and_select_sample

ifeq ($(use_gpu), 1)
#-------------------------------------------------------------------------------
# GPU objectives
#-------------------------------------------------------------------------------
# cu -> cpp
.lv2_gpu_functions_$(SUFFIX).cpp: lv2_gpu_functions.cu lv2_gpu_functions.h $(DEPS)
	$(NVCC) $(DEFINES) $(SM_TARGETS) -o .lv2_gpu_functions_$(SUFFIX).cpp lv2_gpu_functions.cu $(NVCCFLAGS) $(CPU_ARCH) $(INC) $(LIBS) -O3 -DTUNE_ARCH=$(SM_ARCH) -DTUNE_SIZE=$(TUNE_SIZE)

# cpp -> o
.lv2_gpu_functions_$(SUFFIX).o: .lv2_gpu_functions_$(SUFFIX).cpp $(DEPS)
	$(CC) $(CFLAGS) -c .lv2_gpu_functions_$(SUFFIX).cpp -o .lv2_gpu_functions_$(SUFFIX).o

.cx1_functions.o: cx1_functions.cpp lv2_cpu_sort.h $(DEPS)
	$(CC) $(CFLAGS) -c cx1_functions.cpp -o .cx1_functions.o

#-------------------------------------------------------------------------------
# GPU Applications
#-------------------------------------------------------------------------------
$(BIN_DIR)sdbg_builder_cuda_$(SUFFIX): sdbg_builder.cpp .cx1_functions.o .lv2_gpu_functions_$(SUFFIX).o $(DEPS)
	$(CC) $(CFLAGS) $(CUDALIBFLAG) sdbg_builder.cpp .lv2_gpu_functions_$(SUFFIX).o .cx1_functions.o options_description.o -o $(BIN_DIR)sdbg_builder_cuda_$(SUFFIX) -lz
	[[ ! -e $(BIN_DIR)sdbg_builder_gpu ]] || rm $(BIN_DIR)sdbg_builder_gpu
	mv $(BIN_DIR)sdbg_builder_cuda_$(SUFFIX) $(BIN_DIR)sdbg_builder_gpu 
endif

#-------------------------------------------------------------------------------
# Build binary directory
#-------------------------------------------------------------------------------
make_bin_dir:
	[[ -d $(BIN_DIR) ]] || mkdir $(BIN_DIR)

clean:
	rm -f *.i* *.cubin *.cu.c *.cudafe* *.fatbin.c *.ptx *.hash *.cu.cpp *.o .*.cpp
