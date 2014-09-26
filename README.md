MEGAHIT
=========

MEGAHIT is a single node assembler for large and complex metagenomics assembly, such as soil NGS reads. It makes used of succinct *de Bruijn* graph to achieve low memory usage, that can fit the whole assembly graph within a single node. However, its goal is not to make memory usage as low as possible, because it leverages all available memory (assigned by `-m` option) to build succinct *de Bruijn* graphs.

Quick Start
----------------
Use CPU-only version of MEGAHIT

```
make
python megahit -h # show the helping manual
python megahit [options] --cpu-only -m <memory_to_use> -l <max_read_len> {-r <reads.fa> | --input_cmd <command>}
```

Use CUDA version of MEGAHIT, with NVCC version 5.5 or higher.

```
make use_gpu=1
python megahit -h # show the helping manual
python megahit [options] -m <memory_to_use> -l <max_read_len> {-r <reads.fa> | --input_cmd <command>}
```

We recommend to set `-m` as large as possible, but remember to leave some space for your server. For example, for a 64G server, use `-m 60000000000`, which is about 56GB. Typically, 56GB memory is quite enough for human guts samples containing 15-30G base-pairs.

Options
------------------------
###Choose *k*
MEGAHIT uses iterative *k*-mer strategy. Minimum *k*, maximum *k* and the step for iteration can be set by options `--k-min`, `--k-max` and `--k-step`. *k* must be odd numbers while the step must be an even number.

###Filter (*k_min*+1)-mer
(*k_min*+1)-mer with multiplicity lower than *d* (default 2, assigned by `--min-count` option) will be discarded. You should be cautious to set *d* less than 2. This will lead to a much larger and noisy graph. We recommend use the default value 2.

###Mercy *k*-mer
This is specially designed for metagenomics assembly. You can disuse it by adding `--no-mercy` option.

License
------------------------

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
