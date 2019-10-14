### 1.2.9 / 2019-10-13
-   Fix segfault triggered by length-zero sequences
-   Fix memory detection problem for some outdated MacOS versions
-   Fix an incorrect assertion in unitig graph refreshing
-   Added `--verbose` to output full log to the screen

### 1.2.8 / 2019-08-10
-   Add intermediate `megahit_core_popcnt` for CPUs that have ABM but not BMI2
-   Allow new assembly task with `--continue`

### 1.2.7 / 2019-07-28
-   Symbol link `megahit_core_no_hw_accel` to `megahit_toolkit` for backward compatibility
-   Better logging and for memory adjustment during SDBG building
-   Attempt to continue SDBG building even user-specified memory size is not sufficient


### 1.2.6 / 2019-07-13
-   Refactored and fixed a bug in local assembler
-   Refactored `megahit` script
-   Obtain total memory size from `os.sysconf`
-   Fixed segmentation fault in Mac OS with clang 4.0
-   Added `--cleaning-rounds` and `--disconnect-ratio` options for more flexible graph cleaning control

### 1.2.5-beta / 2019-06-28 PST
-   Fixed a bug that causes higher memory usage in seq2sdbg
-   Refactor on sequence sorters, edge I/O and contig I/O

### 1.2.4-beta / 2019-05-25 PST
-   Fixed a few memory leaks
-   Use std::vector to replace malloc in SDBG builders
-   Try to fix potential problem caused by benign data race in unitig graph refreshing
-   Faster by using phmap and xxh3

### 1.2.3-beta / 2019-05-12 PST
-   Refactored sequence readers
-   Fixed a bug in SDBG building of large k-mer sizes

### 1.2.2-beta / 2019-04-16 PST
-   Automatically detect POPCNT/BMI2 and select the correct megahit_core binary

### 1.2.1-beta / 2019-03-30 PST
-   Added `--no-hw-accel` option for users whose CPUs do not support BMI2/POPCNT
-   Added `--test` option for testing
-   Compilable with CMake 2.8 and g++4.8

### 1.2.0-beta / 2019-03-24 PST

Heavily refactored the whole project:

-   Remove GPU support
-   Use cmake
-   Use sparsepp to replace IDBA's hash map for better performance in both speed and memory efficiency
-   Use pdep instruction to speed up rand and select
-   Rewrite unitig graph
-   Rewrite the iterate-edge component
-   Rewrite the SDBG library, except for the builder
-   Fixed a bug which may involve too many reads into local assembly

The changes result in a faster and more memory-efficient tool, but have little effect on assembly quality.

### 1.1.4 / 2018-11-01 PST

-   Fixed a bug in mercy edge stage in 1-pass mode

### 1.1.3 / 2018-03-02 PST

-   Fix a bug in atomic bit vector that may cause a race condition

### 1.1.2 / 2017-08-01 HKT

-   Hotfix of an integer overflow bug

### 1.1.1 / 2016-12-08 HKT

-   Added `-f` option to force overwrite output directory
-   Added `--bubble-level` option to control bubble merging; though level 3 (i.e. super bubble) is not mature (default level is 2)
-   Optimized the speed of tips removal

### 1.1-beta / 2016-11-30 HKT

-   Added components to better handle high depth errors
-   Added components to merge super bubbles
-   Fine tuning k-mer sizes to support longer reads (150bp)
-   In general, it produces longer contigs compared to previous versions

### 1.0.6 / 2016-06-06 HKT

-   Fixed a bug that ignores edge multiplicity at all. This bug existed since v1.0.4-beta

### 1.0.5 / 2016-05-17 HKT

-   Removed the requirement for CPU\_thread &gt;= 2.
-   Added the `--tmp-dir` option
-   More user-friendly error message when seeing bad pair-end files
-   Fixed a bug that may stop an assembly earlier

### 1.0.4-beta / 2016-02-16 HKT

-   Faster index of succinct de Bruijn graph via prefix look up
-   Add `--prune-level` 3 and `--prune-depth` options for more aggressive pruning
-   Tune `bulk` parameters
-   Support reads with length &gt;= 65536 bp

### 1.0.3 / 2015-10-11 HKT

-   Hotfix of number of tip nodes in SdBG 32 bit integer overflow

### 1.0.2 / 2015-08-15 HKT

-   Fixed a bug when number of large multiplicities &gt; INT\_MAX
-   Fixed dead loop in local assembler when \# of reads is 0
-   Use `mmap` for edge/sdbg IO
-   Correct the rounding of edge multiplicity from float to integer

### 1.0.1 / 2015-07-31 HKT

-   Fixed number of SdBG edges 32-bit integer overflow.

### 1.0.0-beta / 2015-07-23 HKT

-   `--presets` option: preset parameters for different types of assembly
-   New CPU sorting design: faster kmer counting & graph construction
-   New unitig graph and edge multiplicity design: more accurate assembly
-   Merge bubble carefully at small *k*: reduce the occurrences of bubble collapsing

### 0.3.3-a / 2015-07-04 HKT

-   Hotfix of incorrect max read length of multiple library
-   Hotfix of `--input-cmd` idling

### 0.3.3 / 2015-07-04 HKT

-   Fixed segmentation fault when a read is all N
-   Fixed continue mode: check continue mode before writing binary reads
-   Slightly improve SdBG traversal functions

### 0.3.2-beta / 2015-06-28 HKT

-   fine tune local assembly multi-thread schedules
-   fine tune SdBG builder memory usage

### 0.3.0-beta3 / 2015-06-22 HKT

-   `--verbose` option
-   fixed a bug when reading paired-end reads of different length
-   print assembly stats to the end of the screen message
-   print read libraries info to the screen message, and number of reads in each libraries to the log file

### 0.3.0-beta2 / 2015-06-20 HKT

-   added the missing file `citycrc.h`

### 0.3.0-beta / 2015-06-18 HKT

New features:

-   `--max-read-len` parameter no longer required
-   `--memory` option set to 0.9 by default
-   make use of PE informations (with local assembly)
-   `--prune-level` and `--merge-level` for setting pruning and merging intensity
-   `--kmin-1pass` option for assembling ultra low-coverage datasets in less memory
-   supporting bzip2 input files
-   useful tools in `megahit_toolkit`, including contig2fastg for conversion of contig files into SPAdes-like fastg

### 0.2.1 / 2015-03-18

Bug Fixes:

-   Fixed incorrect mercy kmers searching when read length &gt;= 255
-   Minor bugs (huge log in some edge cases, etc.)

New features:

-   Semi-auto memory setting: when 0 &lt; "-m" &lt; 1, use fraction of the machine's memory
-   `--out-prefix` option
-   `--cpu-only` turn on by default, and `--use-gpu` option to enable GPU
-   python3 compatibility

### 0.2.0 / 2015-01-30

Bug Fixes:

-   Fixed "option --num-cpu-threads not recognized"

Enhancements:

-   `--mem-flag` option for memory control
-   `--continue` option to resume an interrupted run
-   support mixed fasta/fastq input via `kseq.h`

### 0.1.4 / 2015-01-20

Bug Fixes:

-   Fixed crashes related to OpenMP

### 0.1.3 / 2014-12-01

Enhancements:

-   MAC OS X support
-   Minor improvement to reduce memory usage of the SdBG builder

Bug Fixes:

-   Fixed crashes in some edge cases

### 0.1.2 / 2014-10-10

Enhancements:

-   Update Makefile and python wrapper to improve compatibility
-   Output exit code when subprogram exit abnormally
-   Improve memory usage of the subprogram `assembler`
-   Fix the issue of minor stat differences caused by loops

Bug Fixes:

-   Use `get_ref_with_lock()` to ensure `hash_map` being thread-safe when updating values
-   Correct the computation of edge multiplicity when iterate from small *k* to large *k*
-   Fix a bug that cause segfault in the subprogram `sdbg_builder`

### 0.1.1 beta / 2014-10-02

Enhancements:

-   Add change log
-   More detailed README for input format
-   Use `CompactSequence` in `UnitigGraph`
-   Remove unused parallel sorting codes

Bug Fixes:

-   Fixed wrong computation of `word_per_read` in `cx1_functions.cpp`
-   Fixed crash in `FastxReader` if the file is empty
-   Fixed floating point error in `assembly_algorithms.cpp`

