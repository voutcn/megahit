[![Build Status](https://travis-ci.org/voutcn/megahit.svg?branch=master)](https://travis-ci.org/voutcn/megahit)
[![Build Status](https://drone.io/github.com/voutcn/megahit/status.png)](https://drone.io/github.com/voutcn/megahit/latest)

## Getting Started

```
git clone https://github.com/voutcn/megahit.git
cd megahit
make
./megahit -1 pe_1.fq.gz -2 pe_2.fq.gz -o megahit_out
```

## Introduction
MEGAHIT is a single node assembler for large and complex metagenomics NGS reads, such as soil. It makes use of succinct *de Bruijn* graph (SdBG) to achieve low memory assembly. MEGAHIT can **optionally** utilize a CUDA-enabled GPU to accelerate its SdBG contstruction. The GPU-accelerated version of MEGAHIT has been tested on NVIDIA GTX680 (4G memory) and Tesla K40c (12G memory) with CUDA 5.5, 6.0 and 6.5. MEGAHIT v1.0 or greater also supports IBM Power PC and has been tested on IBM POWER8.

## Dependency & Installation
MEGAHIT is suitable for 64-bit Linux and MAC OS X. It requires [zlib](http://www.zlib.net/), python 2.6 or greater and G++ 4.4 or greater (with `-std=c++0x` and [OpenMP](http://openmp.org) support).
Notably, for MAC OS X, the `g++` in the path is probably the sym-link of `clang`, which do not support OpenMP. Users should have the "real" G++ installed and use `make CXX=/PATH/TO/G++` to specify the compiler.

The GPU counterpart further requires [CUDA](https://developer.nvidia.com/cuda-toolkit) 5.5 or greater. Please use `make use_gpu=1` to compile it, and turn on `--use-gpu` to activate GPU acceleration when running MEGAHIT.

Binary release can be found at the [release page](https://github.com/voutcn/megahit/releases). 

To install MEGAHIT to another directory, please copy *megahit*, *megahit_asm_core*, *megahit_toolkit* and *megahit_sdbg_build* (and *megahit_sdbg_build_gpu* for GPU counterpart) to the destination.

## Running MEGAHIT
If MEGAHIT is successfully compiled, it can be run by the following command:

```
./megahit [options] {-1 <pe_1.fq> -2 <pe_2.fq> | --12 <pe12.fq> | -r <se.fq>}
```

`-1/-2`, `--12` and `-r` are parameters for inputting paired-end, interleaved-paired-end and single-end files. They accept files in fasta (*.fasta*, *.fa*, *.fna*) or fastq (*.fastq*, *.fq*) formats. They also supports gzip files (with *.gz* extensions) and bzip2 files (with *.bz2* extensions). Please run `./megahit -h` for detailed usage message.

##Assembly Tips
MEGAHIT `--presets` option provides five preset parameter combinations:

| Presets | Targeting applications |
|---|---|
| `meta` | General metagenome assembly, such as guts |
| `meta-sensitive` | More sensitive metagenome assembly, but slower |
| `meta-large` | Large and complex metagenome assembly, such as soil |
| `bulk` | Experimental; assembly of standard bulk sequencing with sufficient depth |
| `single-cell` | Experimental; single-cell sequence assembly |

To fine tune parameters for specific datasets, please find our suggestions on [this wiki page](https://github.com/voutcn/megahit/wiki/Assembly-Tips).

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