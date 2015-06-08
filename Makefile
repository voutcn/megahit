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
# make <target>[use_gpu=<0|1>] [disablempopcnt=<0|1>] [sm=<XXX,...>] [abi=<0|1>] [open64=<0|1>] [verbose=<0|1>] [keep=<0|1>]
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
    SM_ARCH = 350
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
CUDALIBFLAG = -L/usr/local/cuda/lib64/ -lcuda -lcudart
CFLAGS = -I. -O2 -Wall -funroll-loops -fprefetch-loop-arrays -fopenmp -std=c++0x -static-libgcc
LIB = -lpthread -lm -lz
ifneq ($(disablempopcnt), 1)
	CFLAGS += -mpopcnt
endif
DEPS = Makefile

#-------------------------------------------------------------------------------
# CPU & GPU version
#-------------------------------------------------------------------------------
ifeq ($(use_gpu), 1)
all:  megahit_assemble megahit_iter sdbg_builder_gpu sdbg_builder_cpu megahit_local_asm
	chmod +x ./megahit
else
all:  megahit_assemble megahit_iter sdbg_builder_cpu megahit_local_asm
	chmod +x ./megahit
endif

#-------------------------------------------------------------------------------
# IDBA library
#-------------------------------------------------------------------------------
LIB_IDBA_DIR = lib_idba
LIB_IDBA = $(LIB_IDBA_DIR)/contig_graph.o
LIB_IDBA += $(LIB_IDBA_DIR)/contig_graph_branch_group.o
LIB_IDBA += $(LIB_IDBA_DIR)/contig_info.o
LIB_IDBA += $(LIB_IDBA_DIR)/hash_graph.o
LIB_IDBA += $(LIB_IDBA_DIR)/sequence.o

#-------------------------------------------------------------------------------
# CPU objectives
#-------------------------------------------------------------------------------
%.o: %.cpp %.h $(DEPS)
	$(CXX) $(CFLAGS) -c $< -o $@
%.o: %.cpp $(DEPS)
	$(CXX) $(CFLAGS) -c $< -o $@

#-------------------------------------------------------------------------------
# CPU Applications
#-------------------------------------------------------------------------------
sdbg_builder_cpu: sdbg_builder.cpp cx1.h lv2_cpu_sort.h cx1_kmer_count.o cx1_edge2sdbg.o cx1_read2sdbg_s1.o cx1_read2sdbg_s2.o cx1_sequences2sdbg.o options_description.o sequence_manager.o $(DEPS)
	$(CXX) $(CFLAGS) sdbg_builder.cpp cx1_kmer_count.o options_description.o cx1_edge2sdbg.o cx1_read2sdbg_s1.o cx1_read2sdbg_s2.o cx1_sequences2sdbg.o sequence_manager.o $(LIB) -o sdbg_builder_cpu

# sdbg_builder_cpu_1pass: sdbg_builder_1pass.cpp cx1_functions_1pass.o lv2_cpu_sort.h options_description.o $(DEPS)
# 	$(CXX) $(CFLAGS) sdbg_builder_1pass.cpp cx1_functions_1pass.o options_description.o $(LIB) -o sdbg_builder_cpu_1pass

megahit_assemble: assembler.cpp succinct_dbg.o rank_and_select.h assembly_algorithms.o branch_group.o options_description.o unitig_graph.o $(DEPS)
	$(CXX) $(CFLAGS) assembler.cpp succinct_dbg.o assembly_algorithms.o branch_group.o options_description.o unitig_graph.o $(LIB) -o megahit_assemble

megahit_iter: iterate_edges.cpp iterate_edges.h options_description.o city.o $(DEPS)
	$(CXX) $(CFLAGS) iterate_edges.cpp options_description.o city.o $(LIB) -o megahit_iter

megahit_local_asm: local_assembler.o local_assemble.cpp city.o options_description.o $(LIB_IDBA) $(DEPS)
	$(CXX) $(CFLAGS) local_assemble.cpp local_assembler.o options_description.o city.o $(LIB_IDBA) $(LIB) -o megahit_local_asm

#-------------------------------------------------------------------------------
# Applications for debug usage
#-------------------------------------------------------------------------------
query_sdbg: query_sdbg.cpp succinct_dbg.o rank_and_select.h assembly_algorithms.o branch_group.o unitig_graph.o $(DEPS)
	$(CXX) $(CFLAGS) query_sdbg.cpp succinct_dbg.o assembly_algorithms.o branch_group.o unitig_graph.o -o query_sdbg

#-------------------------------------------------------------------------------
# GPU objectives
#-------------------------------------------------------------------------------
# cu -> cpp
.lv2_gpu_functions_$(SUFFIX).cpp: lv2_gpu_functions.cu lv2_gpu_functions.h $(DEPS)
	$(NVCC) $(DEFINES) $(SM_TARGETS) lv2_gpu_functions.cu $(NVCCFLAGS) $(CPU_ARCH) $(INC) $(LIBS) -O3 -DTUNE_ARCH=$(SM_ARCH) -DTUNE_SIZE=$(TUNE_SIZE) -o .lv2_gpu_functions_$(SUFFIX).cpp 

# cpp -> o
lv2_gpu_functions_$(SUFFIX).o: .lv2_gpu_functions_$(SUFFIX).cpp $(DEPS)
	$(CXX) $(CFLAGS) -c .lv2_gpu_functions_$(SUFFIX).cpp -o lv2_gpu_functions_$(SUFFIX).o

cx1_kmer_count_gpu.o: cx1_kmer_count.cpp $(DEPS)
	$(CXX) $(CFLAGS) -D USE_GPU -c cx1_kmer_count.cpp -o cx1_kmer_count_gpu.o

cx1_read2sdbg_s1_gpu.o: cx1_read2sdbg_s1.cpp $(DEPS)
	$(CXX) $(CFLAGS) -D USE_GPU -c cx1_read2sdbg_s1.cpp -o cx1_read2sdbg_s1_gpu.o

cx1_read2sdbg_s2_gpu.o: cx1_read2sdbg_s2.cpp $(DEPS)
	$(CXX) $(CFLAGS) -D USE_GPU -c cx1_read2sdbg_s2.cpp -o cx1_read2sdbg_s2_gpu.o

cx1_edge2sdbg_gpu.o: cx1_edge2sdbg.cpp $(DEPS)
	$(CXX) $(CFLAGS) -D USE_GPU -c cx1_edge2sdbg.cpp -o cx1_edge2sdbg_gpu.o

cx1_sequences2sdbg_gpu.o: cx1_edge2sdbg.cpp $(DEPS)
	$(CXX) $(CFLAGS) -D USE_GPU -c cx1_sequences2sdbg.cpp -o cx1_sequences2sdbg_gpu.o

#-------------------------------------------------------------------------------
# GPU Applications
#-------------------------------------------------------------------------------
sdbg_builder_gpu: sdbg_builder.cpp cx1_kmer_count_gpu.o cx1_edge2sdbg_gpu.o cx1_read2sdbg_s1_gpu.o cx1_read2sdbg_s2_gpu.o cx1_sequences2sdbg_gpu.o lv2_gpu_functions_$(SUFFIX).o options_description.o sequence_manager.o $(DEPS)
	$(CXX) $(CFLAGS) $(CUDALIBFLAG) -D USE_GPU sdbg_builder.cpp lv2_gpu_functions_$(SUFFIX).o cx1_kmer_count_gpu.o cx1_edge2sdbg_gpu.o cx1_read2sdbg_s1_gpu.o cx1_read2sdbg_s2_gpu.o cx1_sequences2sdbg_gpu.o options_description.o sequence_manager.o $(LIB) -o sdbg_builder_gpu

#-------------------------------------------------------------------------------
# Build binary directory
#-------------------------------------------------------------------------------

.PHONY:
test: megahit_assemble megahit_iter sdbg_builder_cpu
	-rm -fr example/megahit_out
	./megahit -m 0.9 -l 100 -r example/readsInterleaved.fa -o example/megahit_out

test_gpu: megahit_assemble megahit_iter sdbg_builder_gpu
	-rm -fr example/megahit_gpu_out
	./megahit -m 0.9 -l 100 -r example/readsInterleaved.fa --use-gpu -o example/megahit_gpu_out

.PHONY:
clean:
	-rm -fr *.i* *.cubin *.cu.c *.cudafe* *.fatbin.c *.ptx *.hash *.cu.cpp *.o .*.o .*.cpp \
		$(LIB_IDBA) \
		example/megahit_*out \
		megahit_assemble megahit_iter sdbg_builder_cpu sdbg_builder_gpu sdbg_builder_cpu_1pass sdbg_builder_gpu_1pass
