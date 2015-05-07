MEGAHIT
=========

MEGAHIT is a single node assembler for large and complex metagenomics NGS reads, such as soil. It makes use of succinct *de Bruijn* graph (SdBG) to achieve low memory assembly. The graph construction algorithm can self-adjust to use all available or moderate memory, and can be accelerated if a CUDA-enable GPU is provided. The GPU-accelerated version of MEGAHIT has been tested on NVIDIA GTX680 (4G memory) and Tesla K40c (12G memory).

[![Build Status](https://travis-ci.org/voutcn/megahit.svg)](https://travis-ci.org/voutcn/megahit)
[![Build Status](https://drone.io/github.com/voutcn/megahit/status.png)](https://drone.io/github.com/voutcn/megahit/latest)

Getting Started
----------------

The objective of this branch is using SdBG builder to generate sorted kmers (including reverse complement and dummies) to a file.

### Example (run with 12 threads, 60GB memory)
```
git clone https://github.com/voutcn/megahit.git
cd megahit
git checkout keep-sorted-kmer
make sdbg_builder_cpu edge_reader

# counting
gunzip -cd YOUR_READS.fq.gz | sdbg_builder_cpu count -k 61 -m 1 --host_mem 60000000000 --max_read_length 100 --num_cpu_threads 12 --num_output_threads 4 --output_prefix test 61 --input_file -
# sorting (RC & dummy included); --input_prefix should be the same as --output_prefix above; --num_edge_files should be the same as --num_output_threads above
sdbg_builder_cpu build --host_mem 60000000000 -t 12 --input_prefix test --num_edge_files 4 -o test
# take a look
edge_reader 61 test.sortedEdges | less
```

Citing MEGAHIT
-----------------------
* Li, D., Liu, C-M., Luo, R., Sadakane, K., and Lam, T-W., (2015) MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics*, doi: 10.1093/bioinformatics/btv033 [PMID: [25609793](http://www.ncbi.nlm.nih.gov/pubmed/25609793)].

License
-----------------------
MEGAHIT is released under GPLv3. Several third-party libs are used, including:

* [CUB](https://github.com/NVlabs/cub) under "New BSD"" license
* kseq.h in [klib](https://github.com/attractivechaos/klib) under MIT license
* hash_{table, set, map}.h in [IDBA package](http://i.cs.hku.hk/~alse/hkubrg/projects/idba/) under GPLv2