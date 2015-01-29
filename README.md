MEGAHIT
=========

MEGAHIT is a single node assembler for large and complex metagenomics NGS reads, such as soil. It makes use of succinct *de Bruijn* graph (SdBG) to achieve low memory assembly. The graph construction algorithm can self-adjust to use all available or moderate memory, and can be accelerated if a CUDA-enable GPU is provided. The GPU-accelerated version of MEGAHIT has been tested on NVIDIA GTX680 (4G memory) and Tesla K40c (12G memory).

[![Build Status](https://travis-ci.org/voutcn/megahit.svg)](https://travis-ci.org/voutcn/megahit)
[![Build Status](https://drone.io/github.com/voutcn/megahit/status.png)](https://drone.io/github.com/voutcn/megahit/latest)

Getting Started
----------------

### Dependency
MEGAHIT is suitable for 64-bit Linux and MAC OS X. It requires [zlib](http://www.zlib.net/), python 2.6 or greater and G++ 4.4 or greater (with `-std=c++0x` and [OpenMP](http://openmp.org) support). The GPU version further requires [CUDA](https://developer.nvidia.com/cuda-toolkit) 5.5 or greater.

### Compiling from Source Codes
```
git clone https://github.com/voutcn/megahit.git
cd megahit
make
```

Notably, for MAC OS X, the `g++` in the path is probably the sym-link of `clang`, which do not support OpenMP. Users should have the "real" G++ installed and use `make CXX=/PATH/TO/G++` to specify the compiler.

### Running MEGAHIT
If MEGAHIT is successfully compiled, it can be run by the following command:

```
./megahit [options] --cpu-only -m <max_memory_to_use> -l <max_read_len> {-r <reads.fa> | --input_cmd <command>}
```

User can also run `./megahit -h` for the usage message.

### Using GPU Version
To use the GPU version, run `make use_gpu=1` to compile MEGAHIT, and run MEGAHIT without `--cpu-only` option. GPU version has only been tested in Linux.


Memory Setting
----------------
Users are requried to set a memory parameter `-m` for MEGAHIT. This parameter specifies the maximum memory that can be used by the SdBG constrution conponent of MEGAHIT. `--mem-flag` is another option for memory control.

### Quick recommendation
Set the `-m` parameter to be 90-95% of the available memory and leave the `--mem-flag` default. For example if your machine have 64G available memory, a proper setting is `-m 60e9`.

### Detail explanation 
Please refer to [this wiki page](https://github.com/voutcn/megahit/wiki/MEGAHIT-Memory-setting).


Input Files
--------------

MEGAHIT accepts one fasta or fastq file as input. The input file can be gzip'ed. Alternatively, you can use the option `--input-cmd` to input reads from multiple files. Following the `--input-cmd` should be a command that outputs all reads to `STDOUT` in fasta or fastq format. A mix of fasta and fastq is also supported from version 0.2.0. Currently pair-end information is not used by MEGAHIT. Therefore pair-end files can be input to MEGAHIT as multiple single-end files. Some examples are shown on [this wiki page](https://github.com/voutcn/megahit/wiki/Input-examples).

Options
------------------------
###Choosing *k*
MEGAHIT uses multiple *k*-mer strategy. Minimum *k*, maximum *k* and the step for iteration can be set by options `--k-min`, `--k-max` and `--k-step` respectively. *k* must be odd numbers while the step must be an even number.

For ultra complex metagenomics data such as soil, a larger *k<sub>min</sub>*, say 27, is recommended to reduce the complexity of the *de Bruijn* graph. Quality trimming is also recommended.

###Filtering (*k<sub>min</sub>*+1)-mer
(*k<sub>min</sub>*+1)-mer with multiplicity lower than *d* (default 2, specified by `--min-count` option) will be discarded. You should be cautious to set *d* less than 2, which will lead to a much larger and noisy graph. We recommend using the default value 2 for metagenomics assembly. If you want to use MEGAHIT to do generic assembly, please change this value according to the sequencing depth.

###Mercy *k*-mer
This is specially designed for metagenomics assembly to recover low coverage sequence. You can disable it with `--no-mercy` option.

Issue Report
-----------------------
If you have any questions or suggestions, please [report an issus](https://github.com/voutcn/megahit/issues) on Github.

Cite MEGAHIT
-----------------------
* Li, D., *et al*. (2015) MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics*, doi: 10.1093/bioinformatics/btv033 [PMID: [25609793](http://www.ncbi.nlm.nih.gov/pubmed/25609793)].

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