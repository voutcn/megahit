#-------------------------------------------------------------------------------
# MEGAHIT
# Copyright (C) 2014-2015 The University of Hong Kong & L3 Bioinformatics Limited
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
# make <target>[use_gpu=<0|1>] [version=xxx] [disablempopcnt=<0|1>] [sm=<XXX,...>] [abi=<0|1>] [open64=<0|1>] [verbose=<0|1>] [keep=<0|1>]
#
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# use GPU or not
#-------------------------------------------------------------------------------

ifeq ($(use_gpu), 1)
	NVCC = "$(shell which nvcc)"
	NVCC_VERSION = $(strip $(shell nvcc --version | grep release | sed 's/.*release //' |  sed 's/,.*//'))
endif

version = $(shell git describe --tag 2>/dev/null || echo "git_not_found" 2>/dev/null)

# detect OS
OSUPPER = $(shell uname -s 2>/dev/null | tr [:lower:] [:upper:])

# force 64bits
CPU_ARCH = -m64
CPU_ARCH_SUFFIX = x86_64

IS_PPC64 := $(shell echo `$(CXX) -v 2>&1 | grep powerpc64 | wc -l`)
ifneq (0, $(IS_PPC64))
	CPU_ARCH_SUFFIX = ppc64
	CPU_ARCH = -mpowerpc64
endif

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------

CUDA_INC = "$(shell dirname $(NVCC))/../include"
CUDA_LIB = "$(shell dirname $(NVCC))/../lib64"
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
CUDALIBFLAG = -L$(CUDA_LIB) -lcuda -lcudart
GCC_VER := $(shell echo `$(CXX) -dumpversion | cut -f1-2 -d.`)

CXXFLAGS = -g -O2 -Wall -Wno-unused-function -Wno-array-bounds -D__STDC_FORMAT_MACROS -funroll-loops -fprefetch-loop-arrays -fopenmp -I. -std=c++0x -static-libgcc $(CPU_ARCH)
LIB = -lm -lz -lpthread

ifeq "4.5" "$(word 1, $(sort 4.5 $(GCC_VER)))"
	CXXFLAGS += -static-libstdc++
endif

ifneq ($(disablempopcnt), 1)
	ifeq (0, $(IS_PPC64))
		CXXFLAGS += -mpopcnt
	else
		CXXFLAGS += -mpopcntd
	endif
endif

ifneq ($(version), git_not_found)
	CXXFLAGS += -DPACKAGE_VERSION="\"$(version)\""
endif

#-------------------------------------------------------------------------------
# standalone headers
#-------------------------------------------------------------------------------
STANDALONE_H = rank_and_select.h kmer_plus.h kmer.h lib_info.h \
			   bit_operation.h atomic_bit_vector.h functional.h \
			   khash.h kseq.h pool.h packed_reads.h sequence_package.h \
			   utils.h mem_file_checker-inl.h read_lib_functions-inl.h \
			   edge_io.h histgram.h definitions.h lv2_cpu_sort.h sdbg_multi_io.h \
			   cx1.h

DEPS = Makefile $(STANDALONE_H)

#-------------------------------------------------------------------------------
# CPU & GPU version
#-------------------------------------------------------------------------------
ifeq ($(use_gpu), 1)
all:  megahit_asm_core megahit_sdbg_build_gpu megahit_sdbg_build megahit_toolkit
	chmod +x ./megahit
else
all:  megahit_asm_core megahit_sdbg_build megahit_toolkit
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
# Tookits
#-------------------------------------------------------------------------------
TOOLS_DIR = tools
TOOLKIT = $(TOOLS_DIR)/toolkit.cpp
TOOLKIT += $(TOOLS_DIR)/contigs_to_fastg.cpp
TOOLKIT += $(TOOLS_DIR)/read_stat.cpp
TOOLKIT += $(TOOLS_DIR)/trim_low_qual_tail.cpp
TOOLKIT += $(TOOLS_DIR)/filter_by_len.cpp
TOOLKIT += $(TOOLS_DIR)/extract_pe_reads.cpp

#-------------------------------------------------------------------------------
# CPU objectives
#-------------------------------------------------------------------------------
%.o: %.cpp %.h $(DEPS)
	$(CXX) $(CXXFLAGS) -c $< -o $@
%.o: %.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

#-------------------------------------------------------------------------------
# asm_core objectives
#-------------------------------------------------------------------------------
LIB_ASM = succinct_dbg.o assembly_algorithms.o options_description.o \
				  unitig_graph.o sequence_manager.o local_assembler.o city.o

#-------------------------------------------------------------------------------
# CPU Applications
#-------------------------------------------------------------------------------
megahit_sdbg_build: sdbg_builder.cpp cx1_kmer_count.o cx1_read2sdbg_s1.o cx1_read2sdbg_s2.o cx1_seq2sdbg.o options_description.o sequence_manager.o $(DEPS)
	$(CXX) $(CXXFLAGS) kthread.cpp sdbg_builder.cpp cx1_kmer_count.o options_description.o cx1_read2sdbg_s1.o cx1_read2sdbg_s2.o cx1_seq2sdbg.o sequence_manager.o $(LIB) -o megahit_sdbg_build

megahit_asm_core: $(LIB_ASM) $(LIB_IDBA) asm_core.o assembler.o local_assemble.o iterate_edges.o build_read_lib.o $(DEPS)
	$(CXX) $(CXXFLAGS) asm_core.o assembler.o local_assemble.o iterate_edges.o build_read_lib.o $(LIB_IDBA) $(LIB_ASM) $(LIB) -o megahit_asm_core

megahit_toolkit: $(TOOLKIT) $(DEPS)
	$(CXX) $(CXXFLAGS) $(TOOLKIT) $(LIB) -o megahit_toolkit

#-------------------------------------------------------------------------------
# GPU objectives
#-------------------------------------------------------------------------------
# cu -> cpp
.lv2_gpu_functions_$(SUFFIX).cpp: lv2_gpu_functions.cu lv2_gpu_functions.h $(DEPS)
	$(NVCC) $(DEFINES) $(SM_TARGETS) lv2_gpu_functions.cu $(NVCCFLAGS) $(CPU_ARCH) $(INC) $(LIBS) -O3 -DTUNE_ARCH=$(SM_ARCH) -DTUNE_SIZE=$(TUNE_SIZE) -o .lv2_gpu_functions_$(SUFFIX).cpp 

# cpp -> o
lv2_gpu_functions_$(SUFFIX).o: .lv2_gpu_functions_$(SUFFIX).cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c .lv2_gpu_functions_$(SUFFIX).cpp -o lv2_gpu_functions_$(SUFFIX).o

cx1_kmer_count_gpu.o: cx1_kmer_count.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -D USE_GPU -c cx1_kmer_count.cpp -o cx1_kmer_count_gpu.o

cx1_read2sdbg_s1_gpu.o: cx1_read2sdbg_s1.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -D USE_GPU -c cx1_read2sdbg_s1.cpp -o cx1_read2sdbg_s1_gpu.o

cx1_read2sdbg_s2_gpu.o: cx1_read2sdbg_s2.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -D USE_GPU -c cx1_read2sdbg_s2.cpp -o cx1_read2sdbg_s2_gpu.o

cx1_seq2sdbg_gpu.o: cx1_seq2sdbg.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -D USE_GPU -c cx1_seq2sdbg.cpp -o cx1_seq2sdbg_gpu.o

#-------------------------------------------------------------------------------
# GPU Applications
#-------------------------------------------------------------------------------
megahit_sdbg_build_gpu: sdbg_builder.cpp cx1_kmer_count_gpu.o cx1_read2sdbg_s1_gpu.o cx1_read2sdbg_s2_gpu.o cx1_seq2sdbg_gpu.o lv2_gpu_functions_$(SUFFIX).o options_description.o sequence_manager.o $(DEPS)
	$(CXX) $(CXXFLAGS) $(CUDALIBFLAG) -D USE_GPU kthread.cpp sdbg_builder.cpp lv2_gpu_functions_$(SUFFIX).o cx1_kmer_count_gpu.o cx1_read2sdbg_s1_gpu.o cx1_read2sdbg_s2_gpu.o cx1_seq2sdbg_gpu.o options_description.o sequence_manager.o $(LIB) -o megahit_sdbg_build_gpu

#-------------------------------------------------------------------------------
# Build binary directory
#-------------------------------------------------------------------------------

.PHONY:
test: megahit_asm_core megahit_sdbg_build megahit_toolkit
	-rm -fr example/megahit_out
	./megahit --12 example/readsInterleaved1.fa.gz,example/readsInterleaved2.fa.bz2,example/readsInterleaved3.fa -o example/megahit_out -t 1
	-rm -fr example/megahit_out
	./megahit --12 example/readsInterleaved1.fa.gz,example/readsInterleaved2.fa.bz2,example/readsInterleaved3.fa -o example/megahit_out -t 1 --kmin-1pass

test_gpu: megahit_asm_core megahit_sdbg_build_gpu megahit_toolkit
	-rm -fr example/megahit_gpu_out
	./megahit --12 example/readsInterleaved1.fa.gz,example/readsInterleaved2.fa.bz2,example/readsInterleaved3.fa --use-gpu -o example/megahit_gpu_out -t 1
	-rm -fr example/megahit_gpu_out
	./megahit --12 example/readsInterleaved1.fa.gz,example/readsInterleaved2.fa.bz2,example/readsInterleaved3.fa --use-gpu -o example/megahit_gpu_out -t 1 --kmin-1pass

release_dir = megahit_$(version)_$(OSUPPER)_CUDA$(NVCC_VERSION)_sm$(SM_ARCH)_$(CPU_ARCH_SUFFIX)-bin

.PHONY:
release: megahit_asm_core megahit_sdbg_build megahit_sdbg_build_gpu megahit_toolkit megahit README.md ChangeLog.md
	rm -rf $(release_dir) $(release_dir).tar.gz
	mkdir -p $(release_dir)
	cp megahit_asm_core megahit_sdbg_build megahit_sdbg_build_gpu megahit_toolkit megahit README.md ChangeLog.md \
	   $(release_dir)
	tar zvcf $(release_dir).tar.gz \
	         $(release_dir)

cpu_release_dir = megahit_$(version)_$(OSUPPER)_CPUONLY_$(CPU_ARCH_SUFFIX)-bin

.PHONY:
release_cpu: megahit_asm_core megahit_sdbg_build megahit_toolkit megahit README.md ChangeLog.md
	rm -rf $(cpu_release_dir) $(cpu_release_dir).tar.gz
	mkdir -p $(cpu_release_dir)
	cp megahit_asm_core megahit_sdbg_build megahit_toolkit megahit README.md ChangeLog.md\
	   $(cpu_release_dir)
	tar zvcf $(cpu_release_dir).tar.gz $(cpu_release_dir)

.PHONY:
clean:
	-rm -fr *.i* *.cubin *.cu.c *.cudafe* *.fatbin.c *.ptx *.hash *.cu.cpp *.o .*.o .*.cpp \
		$(LIB_IDBA) \
		example/megahit_*out \
		megahit_asm_core megahit_sdbg_build megahit_sdbg_build_gpu megahit_toolkit
