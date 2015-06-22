MEGAHIT
=========

MEGAHIT is a single node assembler for large and complex metagenomics NGS reads, such as soil. It makes use of succinct *de Bruijn* graph (SdBG) to achieve low memory assembly. The graph construction algorithm can self-adjust to use all available or moderate memory, and can be accelerated if a CUDA-enable GPU is provided. The GPU-accelerated version of MEGAHIT has been tested on NVIDIA GTX680 (4G memory) and Tesla K40c (12G memory) with CUDA 5.5, 6.0 and 6.5.

[![Build Status](https://travis-ci.org/voutcn/megahit.svg)](https://travis-ci.org/voutcn/megahit)
[![Build Status](https://drone.io/github.com/voutcn/megahit/status.png)](https://drone.io/github.com/voutcn/megahit/latest)

##Getting Started

#### Dependency
MEGAHIT is suitable for 64-bit Linux and MAC OS X. It requires [zlib](http://www.zlib.net/), python 2.6 or greater and G++ 4.4 or greater (with `-std=c++0x` and [OpenMP](http://openmp.org) support). The GPU version further requires [CUDA](https://developer.nvidia.com/cuda-toolkit) 5.5 or greater.

#### Compiling from Source Codes
```
git clone https://github.com/voutcn/megahit.git
cd megahit
make
```

Notably, for MAC OS X, the `g++` in the path is probably the sym-link of `clang`, which do not support OpenMP. Users should have the "real" G++ installed and use `make CXX=/PATH/TO/G++` to specify the compiler.

If you need to install MEGAHIT to another directory after compilation, please copy `megahit`, `megahit_asm_core`, `megahit_toolkit` and `megahit_sdbg_build` (and `megahit_sdbg_build_gpu` for GPU counterpart) to the destination.

#### Stable binary releases
Please refer to the [release page](https://github.com/voutcn/megahit/releases).

#### Running MEGAHIT
If MEGAHIT is successfully compiled, it can be run by the following command:

```
./megahit [options] {-1 <pe_1.fq> -2 <pe_2.fq> | --12 <pe12.fq> | -r <se.fq>}
```

The `-1/-2/--12` and `-r` parameters accept files in fasta (`.fasta`, `.fa`, `.fna`) or fastq (`.fastq`, `.fq`) formats. They also supports gzip files (with `.gz` extensions) and bzip2 files (with `.bz2` extensions). Please run `./megahit -h` for detailed usage message.

#### Using GPU Version
To use the GPU version, run `make use_gpu=1` to compile MEGAHIT, and run MEGAHIT with `--use-gpu`. GPU version has only been tested in Linux.

##Assembly Tips

####Choosing *k*
MEGAHIT uses multiple *k*-mer strategy. Minimum *k*, maximum *k* and the step for iteration can be set by options `--k-min`, `--k-max` and `--k-step` respectively. *k* must be odd numbers while the step must be an even number.

* for ultra complex metagenomics data such as soil, a larger *k<sub>min</sub>*, say 27, is recommended to reduce the complexity of the *de Bruijn* graph. Quality trimming is also recommended
* for high-depth generic data, large `--k-min` (25 to 31) is recommended
* smaller `--k-step`, say 10, is more friendly to low-coverage datasets

####Filtering (*k<sub>min</sub>*+1)-mer
(*k<sub>min</sub>*+1)-mer with multiplicity lower than *d* (default 2, specified by `--min-count` option) will be discarded. You should be cautious to set *d* less than 2, which will lead to a much larger and noisy graph. We recommend using the default value 2 for metagenomics assembly. If you want to use MEGAHIT to do generic assembly, please change this value according to the sequencing depth. (recommend `--min-count 3` for >40x).

####Mercy *k*-mer
This is specially designed for metagenomics assembly to recover low coverage sequence. For generic dataset >= 30x,  MEGAHIT may generate better results with `--no-mercy` option.

####*k*-min 1pass mode
This mode can be activated by option `--kmin-1pass`. It is more memory efficient for ultra low-depth datasets, such as soil metagenomics data.


##FAQ & Reporting issues

For other questions, please first refer to [our wiki](https://github.com/voutcn/megahit/wiki). Please [report an issue](https://github.com/voutcn/megahit/issues) in github when necessary.

##Citing MEGAHIT

* Li, D., Liu, C-M., Luo, R., Sadakane, K., and Lam, T-W., (2015) MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics*, doi: 10.1093/bioinformatics/btv033 [PMID: [25609793](http://www.ncbi.nlm.nih.gov/pubmed/25609793)].

##License

```
    MEGAHIT
    Copyright (C) 2014-2015  The University of Hong Kong

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

Third-party codes used in MEGAHIT:

* [IDBA package](https://github.com/loneknightpy/idba) under GPLv2
* [CUB](https://github.com/NVlabs/cub) under "New BSD" license
* [klib](https://github.com/attractivechaos/klib) under MIT license
* [CityHash](https://code.google.com/p/cityhash/) under MIT license