MEGAHIT
=========

MEGAHIT is a single node assembler for large and complex metagenomics NGS reads, such as soil. It makes use of succinct *de Bruijn* graph to achieve low memory usage, whereas its goal is not to make memory usage as low as possible. It leverages all available memory (assigned by `-m` option) to build succinct *de Bruijn* graphs. CPU-only and GPU-accelerated version of MEGAHIT are provided. The GPU-accelerated version of MEGAHIT has been tested on NVIDIA GTX680 (4G memory) and Tesla K40c (12G memory).

Quick Start
----------------
Use CPU-only version of MEGAHIT

```
% make
% python ./megahit [options] --cpu-only -m <memory_to_use> -l <max_read_len> {-r <reads.fa> | --input_cmd <command>}
```

Use GPU-accelerated version of MEGAHIT, with a CUDA-enabled GPU and NVCC version 5.5 or higher.

```
% make use_gpu=1
% python ./megahit [options] -m <memory_to_use> -l <max_read_len> {-r <reads.fa> | --input_cmd <command>}
```

To show the usage message, type the command

```
% python ./megahit -h # show the helping manual
```


Memory control
----------------
We recommend to set `-m` as large as possible. But remember to leave some space for your server. For example, for a server with 64GB free memory, you may try `-m 60000000000`, which is about 56GB. Typically, 56GB memory is quite enough for human guts samples containing 15-30G base-pairs.

Input files
--------------

MEGAHIT accepts **one** fasta or fastq file as input. The input file can be gzip'ed. Alternatively, you can use the option `--input-cmd` to input reads. Following the `--input-cmd` should be a command that output all reads to `STDOUT` in fasta or fastq format. A mixed of fasta and fastq is NOT supported. Some correct/wrong examples below.

###Correct examples
* Input from one fastq file named *reads.fastq*:
```
-r read.fastq
```
* Input from two fasta files with prefix *sample_1.fa* and  *sample_2.fa*:
```
--input-cmd "cat sample_[12].fa"
```
* Input from all gzip'ed fastq files in current directory:
```
--input-cmd "zcat *.fastq.gz"
```
* Assume fastq-dump is installed, input from a sra file *xxx.sra*:
```
--input-cmd "fastq-dump -Z --fasta xxx.sra"
```

###Wrong examples
* Mixed fastq and fasta files to the input:
```
--input-cmd "cat *.fa *.fq"
```

Options
------------------------
###Choose *k*
MEGAHIT uses multiple *k*-mer strategy. Minimum *k*, maximum *k* and the step for iteration can be set by options `--k-min`, `--k-max` and `--k-step` respectively. *k* must be odd numbers while the step must be an even number.

###Filter (*k_min*+1)-mer
(*k_min*+1)-mer with multiplicity lower than *d* (default 2, assigned by `--min-count` option) will be discarded. You should be cautious to set *d* less than 2, which will lead to a much larger and noisy graph. We recommend use the default value 2.

###Mercy *k*-mer
This is specially designed for metagenomics assembly. You can disable this stategy by adding `--no-mercy` option.

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