MEGAHIT
=========

MEGAHIT is a single node assembler for large and complex metagenomics NGS reads, such as soil. It makes use of succinct *de Bruijn* graph to achieve low memory usage, whereas its goal is not to make memory usage as low as possible. It automatically adjusts itself to use all/moderate available memory (specified by `-m` and `--mem-flag`) to build succinct *de Bruijn* graphs. Both CPU-only and GPU-accelerated version of MEGAHIT are provided. The GPU-accelerated version of MEGAHIT has been tested on NVIDIA GTX680 (4G memory) and Tesla K40c (12G memory).

[![Build Status](https://travis-ci.org/voutcn/megahit.svg)](https://travis-ci.org/voutcn/megahit)
[![Build Status](https://drone.io/github.com/voutcn/megahit/status.png)](https://drone.io/github.com/voutcn/megahit/latest)

Getting Started
----------------

### Dependency
MEGAHIT is suitable for 64-bit Linux and MAC OS X. It requires [zlib](http://www.zlib.net/), python 2.6 or greater and GCC 4.4 or greater (with `-std=c++0x` and [OpenMP](http://openmp.org) support). The GPU version of MEGAHIT further requires [CUDA](https://developer.nvidia.com/cuda-toolkit) 5.5 or greater.

### Compiling from Source Codes
```
git clone https://github.com/voutcn/megahit.git
cd megahit
make
```

Notably, for MAC OS X, the `g++` in the path is probably the sym-link of `clang`, which do not support OpenMP. Users should have the "real" GCC installed and use `make CXX=/PATH/TO/G++` to specify the compiler.

### Running MEGAHIT
If MEGAHIT is successfully compiled, it can be run by the following command:

```
./megahit [options] --cpu-only -m <memory_to_use> -l <max_read_len> {-r <reads.fa> | --input_cmd <command>}
```

User can also run `./megahit -h` for the usage message.

### Using GPU Version
To use the GPU version, run `make use_gpu=1` to compile MEGAHIT, and run MEGAHIT without `--cpu-only` option. At this stage, the GPU version only supports Linux.


Memory Control
----------------
We recommend to set `-m` as large as possible. For example, for a server with 64GB free memory, you may try `-m 60e9`, which is about 56GB.
Since v0.1.4, the `-m` is used to set the maximum memory. It is not necessary for the SdBG builder to use up all the memories. The option `--mem-flag` specifies the ways to utilize memory. `--mem-flag 0` to use minimum memory, `--mem-flag 1` to use moderate memory and `--mem-flag 2` to use all memory.

Input Files
--------------

MEGAHIT accepts one fasta or fastq file as input. The input file can be gzip'ed. Alternatively, you can use the option `--input-cmd` to input reads from multiple files. Following the `--input-cmd` should be a command that output all reads to `STDOUT` in fasta or fastq format. A mix of fasta and fastq is also supported. Pair-end information is not used by MEGAHIT of current version. Therefore pair-end files can be input to MEGAHIT as multiple single-end files. Some examples are shown below.

###Correct Examples
* Input from one gzip'ed fastq file named *reads.fastq.gz*:
```
-r reads.fastq.gz
```
* Input from two fasta files *sample_1.fa* and *sample_2.fa*:
```
--input-cmd "cat sample_[12].fa"
```
* Input from all gzip'ed fastq files in current directory:
```
--input-cmd "zcat *.fastq.gz"
```
* Assume that fastq-dump is installed, input from a sra file *xxx.sra*:
```
--input-cmd "fastq-dump -Z --fasta xxx.sra"
```
* Mixed fastq and fasta:
```
--input-cmd "cat 1.fa 2.fq"
```

Options
------------------------
###Choosing *k*
MEGAHIT uses multiple *k*-mer strategy. Minimum *k*, maximum *k* and the step for iteration can be set by options `--k-min`, `--k-max` and `--k-step` respectively. *k* must be odd numbers while the step must be an even number.

For ultra complex metagenomics data such as soil, a larger *k<sub>min</sub>*, say 27, is recommended to reduce the complexity of the *de Bruijn* graph.

###Filtering (*k<sub>min</sub>*+1)-mer
(*k<sub>min</sub>*+1)-mer with multiplicity lower than *d* (default 2, specified by `--min-count` option) will be discarded. You should be cautious to set *d* less than 2, which will lead to a much larger and noisy graph. We recommend using the default value 2.

###Mercy *k*-mer
This is specially designed for metagenomics assembly to recover low coverage sequence. You can disable this strategy by adding `--no-mercy` option.

License
-----------------------
```
  MEGAHIT
  
  Copyright (C) 2014 The University of Hong Kong

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
```