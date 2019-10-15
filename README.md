MEGAHIT
=======

[![BioConda Install](https://img.shields.io/conda/dn/bioconda/megahit.svg?style=flat-square&label=BioConda%20install)](https://anaconda.org/bioconda/megahit)
[![Downloads](https://img.shields.io/github/downloads/voutcn/megahit/total?style=flat-square)](https://github.com/voutcn/megahit/releases)
[![Build Status](https://img.shields.io/travis/voutcn/megahit?style=flat-square)](https://travis-ci.org/voutcn/megahit)
[![codecov](https://img.shields.io/codecov/c/github/voutcn/megahit?style=flat-square)](https://codecov.io/gh/voutcn/megahit)

MEGAHIT is an ultra-fast and memory-efficient NGS assembler. It is optimized for metagenomes, but also works well on generic single genome assembly (small or mammalian size) and single-cell assembly.

Installation
---------------

### Conda
```sh
conda install -c bioconda megahit
```

### Pre-built binaries for x86_64 Linux

```sh
wget https://github.com/voutcn/megahit/releases/download/v1.2.9/MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz
tar zvxf MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz
cd MEGAHIT-1.2.9-Linux-x86_64-static/bin/
./megahit --test  # run on a toy dataset
./megahit -1 MY_PE_READ_1.fq.gz -2 MY_PE_READ_2.fq.gz -o MY_OUTPUT_DIR
```

### Pre-built docker image
``` sh
# in the directory with the input reads
docker run -v $(pwd):/workspace -w /workspace --user $(id -u):$(id -g) vout/megahit \
  megahit -1 MY_PE_READ_1.fq.gz -2 MY_PE_READ_2.fq.gz -o MY_OUTPUT_DIR
```

### Building from source

#### Prerequisites

-   For building: zlib, cmake &gt;= 2.8, g++ &gt;= 4.8.4
-   For running: gzip and bzip2

```sh
git clone https://github.com/voutcn/megahit.git
cd megahit
git submodule update --init
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release  # add -DCMAKE_INSTALL_PREFIX=MY_PREFIX if needed
make -j4
make simple_test  # will test MEGAHIT with a toy dataset
# make install if needed
```

Usage
-----

### Basic usage
```sh
megahit -1 pe_1.fq -2 pe_2.fq -o out  # 1 paired-end library
megahit --12 interleaved.fq -o out # one paired & interleaved paired-end library
megahit -1 a1.fq,b1.fq,c1.fq -2 a2.fq,b2.fq,c2.fq -r se1.fq,se2.fq -o out # 3 paired-end libraries + 2 SE libraries
megahit_core contig2fastg 119 out/intermediate_contigs/k119.contig.fa > k119.fastg # get FASTG from the intermediate contigs of k=119
```
The contigs can be found `final.contigs.fa` in the output directory.

### Advanced usage
- `--kmin-1pass`: if sequencing depth is low and too much memory used when build the graph of k_min
- `--presets meta-large`: if the metagenome is complex (i.e., bio-diversity is high, for example soil metagenomes)
- `--cleaning-rounds 1 --disconnect-ratio 0`: get less pruned assembly (usually shorter contigs)
- `--continue -o out`: resume an interrupted job from `out`

To see the full manual, run `megahit` without parameters or with `-h`.

Also, our [wiki](https://github.com/voutcn/megahit/wiki) may be helpful.

Publications
------------

-   Li, D., Liu, C-M., Luo, R., Sadakane, K., and Lam, T-W., (2015) MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics*, doi: 10.1093/bioinformatics/btv033 \[PMID: [25609793](http://www.ncbi.nlm.nih.gov/pubmed/25609793)\].
-   Li, D., Luo, R., Liu, C.M., Leung, C.M., Ting, H.F., Sadakane, K., Yamashita, H. and Lam, T.W., 2016. MEGAHIT v1.0: A Fast and Scalable Metagenome Assembler driven by Advanced Methodologies and Community Practices. Methods.

License
-------

This project is licensed under the GPLv3 License - see the [LICENSE](LICENSE) file for details
