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
# make <target> [version=xxx] [disablempopcnt=<0|1>]
#
#-------------------------------------------------------------------------------

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

DEPS =   ./Makefile

#-------------------------------------------------------------------------------
# g++ and its options
#-------------------------------------------------------------------------------
GCC_VER := $(shell echo `$(CXX) -dumpversion | cut -f1-2 -d.`)

CXXFLAGS = -O2 -g -Wall -Wno-unused-function -Wno-array-bounds -D__STDC_FORMAT_MACROS -DXXH_INLINE_ALL -funroll-loops -fprefetch-loop-arrays -fopenmp -I. -std=c++0x -static-libgcc $(CPU_ARCH)
LIB = -lm -lz -lpthread

ifeq "4.5" "$(word 1, $(sort 4.5 $(GCC_VER)))"
	CXXFLAGS += -static-libstdc++
endif

ifeq ($(static), 1)
	CXXFLAGS += -static
endif

ifneq ($(disablempopcnt), 1)
	ifeq (0, $(IS_PPC64))
		CXXFLAGS += -mpopcnt -mbmi2
	else
		CXXFLAGS += -mpopcntd -mbmi2
	endif
endif

ifneq ($(version), git_not_found)
	CXXFLAGS += -DPACKAGE_VERSION="\"$(version)\""
endif

#-------------------------------------------------------------------------------
# standalone headers
#-------------------------------------------------------------------------------
STANDALONE_H = kmlib/kmrns.h kmer_plus.h kmer.h lib_info.h \
			   bit_operation.h \
			   kseq.h packed_reads.h sequence_package.h \
			   utils.h mem_file_checker-inl.h read_lib_functions-inl.h \
			   edge_io.h histgram.h definitions.h lv2_cpu_sort.h \
			   cx1.h

DEPS = Makefile $(STANDALONE_H)

all:  megahit_asm_core megahit_sdbg_build megahit_toolkit
	chmod +x ./megahit
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
LIB_ASM_ALGO = assembly/unitig_graph.o assembly/sdbg_pruning.o \
	           assembly/bubble_remover.o assembly/tip_remover.o \
	           assembly/weak_link_remover.o assembly/low_depth_remover.o \
	           assembly/contig_writer.o
LIB_ASM = sdbg/sdbg_writer.o sdbg/sdbg_raw_content.o \
          sdbg/sdbg_meta.o options_description.o \
	      sequence_manager.o local_assembler.o \
	      asm_core.o assembler.o local_assemble.o iterate_edges.o build_read_lib.o \
	      $(LIB_ASM_ALGO)

#-------------------------------------------------------------------------------
# CPU Applications
#-------------------------------------------------------------------------------
megahit_sdbg_build: sdbg_builder.cpp cx1_kmer_count.o cx1_read2sdbg_s1.o cx1_read2sdbg_s2.o cx1_seq2sdbg.o options_description.o sequence_manager.o  sdbg/sdbg_writer.o sdbg/sdbg_meta.o $(DEPS)
	$(CXX) $(CXXFLAGS) kthread.cpp sdbg_builder.cpp cx1_kmer_count.o options_description.o cx1_read2sdbg_s1.o cx1_read2sdbg_s2.o cx1_seq2sdbg.o sequence_manager.o  sdbg/sdbg_writer.o sdbg/sdbg_meta.o $(LIB) -o megahit_sdbg_build

megahit_asm_core: $(LIB_ASM) $(LIB_IDBA) $(DEPS)
	$(CXX) $(CXXFLAGS) $(LIB_IDBA) $(LIB_ASM) $(LIB) -o megahit_asm_core

megahit_toolkit: $(TOOLKIT) $(DEPS)
	$(CXX) $(CXXFLAGS) $(TOOLKIT) $(LIB) -o megahit_toolkit

#-------------------------------------------------------------------------------
# Build binary directory
#-------------------------------------------------------------------------------

.PHONY:
test: megahit_asm_core megahit_sdbg_build megahit_toolkit
	-rm -fr example/megahit_out
	./megahit --12 example/readsInterleaved1.fa.gz,example/readsInterleaved2.fa.bz2,example/readsInterleaved3.fa -o example/megahit_out -t 2
	-rm -fr example/megahit_out
	./megahit --12 example/readsInterleaved1.fa.gz,example/readsInterleaved2.fa.bz2,example/readsInterleaved3.fa -o example/megahit_out -t 2 --kmin-1pass

.PHONY:
release: megahit_asm_core megahit_sdbg_build megahit_toolkit megahit README.md ChangeLog.md
	rm -rf $(cpu_release_dir) $(cpu_release_dir).tar.gz
	mkdir -p $(cpu_release_dir)
	cp megahit_asm_core megahit_sdbg_build megahit_toolkit megahit README.md ChangeLog.md\
	   $(cpu_release_dir)
	tar zvcf $(cpu_release_dir).tar.gz $(cpu_release_dir)

.PHONY:
clean:
	-rm -fr *.o .*.o .*.cpp \
		$(LIB_IDBA) \
		$(LIB_ASM) \
		example/megahit_*out \
		megahit_asm_core megahit_sdbg_build megahit_toolkit

INSTALL	?= install
PREFIX	?= /usr/local
DESTDIR	?= .
STRIP	?= strip

.PHONY:
install:
	$(INSTALL) -d $(DESTDIR)$(PREFIX)/bin
	$(INSTALL) -c megahit megahit_* $(DESTDIR)$(PREFIX)/bin

.PHONY:
install-strip: install
	$(STRIP) $(DESTDIR)$(PREFIX)/bin/megahit_*
