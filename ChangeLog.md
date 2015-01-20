### 0.1.4 / 2015-01-20
Bug Fixes:

* Fixed crashes related to OpenMP

### 0.1.3 / 2014-12-01

Enhancements:

* MAC OS X support (Thank [Laurent Gautier](https://github.com/lgautier))
* Minor improvement to reduce memory usage of the SdBG builder

Bug Fixes:

* Fixed crashes in some edge cases

### 0.1.2 / 2014-10-10

Enhancements:

* Update Makefile and python wrapper to improve compatibility (Thank [C. Titus Brown](https://github.com/ctb) and [Ben J Woodcroft](https://github.com/wwood))
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