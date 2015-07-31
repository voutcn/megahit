### 1.0.0-beta / 2015-07.23 HKT

* `--presets` option: preset parameters for different types of assembly
* New CPU sorting design: faster kmer counting & graph construction
* New unitig graph and edge multiplicity design: more accurate assembly
* Merge bubble carefully at small *k*: reduce the occurrences of bubble collapsing

### 0.3.3 / 2015-07-04 HKT

* Fixed segmentation fault when a read is all N
* Fixed continue mode: check continue mode before writing binary reads
* Slightly improve SdBG traversal functions

### 0.3.2-beta / 2015-06-28 HKT

* fine tune local assembly multi-thread schedules
* fine tune SdBG builder memory usage

### 0.3.0-beta3 / 2015-06-22 HKT

* `--verbose` option
* fixed a bug when reading paired-end reads of different length
* print assembly stats to the end of the screen message
* print read libraries info to the screen message, and number of reads in each libraries to the log file

### 0.3.0-beta2 / 2015-06-20 HKT

* added the missing file `citycrc.h`

### 0.3.0-beta / 2015-06-18 HKT

New features:

* `--max-read-len` parameter no longer required
* `--memory` option set to 0.9 by default
* make use of PE informations (with local assembly)
* `--prune-level` and `--merge-level` for setting pruning and merging intensity
* `--kmin-1pass` option for assembling ultra low-coverage datasets in less memory
* supporting bzip2 input files
* useful tools in `megahit_toolkit`, including contig2fastg for conversion of contig files into SPAdes-like fastg

### 0.2.1 / 2015-03-18
Bug Fixes:

* Fixed incorrect mercy kmers searching when read length >= 255
* Minor bugs (huge log in some edge cases, etc.)

New features:

* Semi-auto memory setting: when 0 < "-m" < 1, use fraction of the machine's memory
* `--out-prefix` option
* `--cpu-only` turn on by default, and `--use-gpu` option to enable GPU
* python3 compatibility

### 0.2.0 / 2015-01-30
Bug Fixes:

* Fixed "option --num-cpu-threads not recognized"

Enhancements:

* `--mem-flag` option for memory control
* `--continue` option to resume an interrupted run
* support mixed fasta/fastq input via `kseq.h`

### 0.1.4 / 2015-01-20
Bug Fixes:

* Fixed crashes related to OpenMP

### 0.1.3 / 2014-12-01

Enhancements:

* MAC OS X support
* Minor improvement to reduce memory usage of the SdBG builder

Bug Fixes:

* Fixed crashes in some edge cases

### 0.1.2 / 2014-10-10

Enhancements:

* Update Makefile and python wrapper to improve compatibility
* Output exit code when subprogram exit abnormally
* Improve memory usage of the subprogram `assembler`
* Fix the issue of minor stat differences caused by loops

Bug Fixes:

* Use `get_ref_with_lock()` to ensure `hash_map` being thread-safe when updating values
* Correct the computation of edge multiplicity when iterate from small *k* to large *k*
* Fix a bug that cause segfault in the subprogram `sdbg_builder`


### 0.1.1 beta / 2014-10-02

Enhancements:

* Add change log
* More detailed README for input format
* Use `CompactSequence` in `UnitigGraph`
* Remove unused parallel sorting codes

Bug Fixes:

* Fixed wrong computation of `word_per_read` in `cx1_functions.cpp`
* Fixed crash in `FastxReader` if the file is empty
* Fixed floating point error in `assembly_algorithms.cpp`