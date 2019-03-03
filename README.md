[![BioConda Install](https://img.shields.io/conda/dn/bioconda/megahit.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/megahit)
[![Build Status](https://travis-ci.org/voutcn/megahit.svg?branch=master)](https://travis-ci.org/voutcn/megahit)

## Getting Started

### Run with Docker
```bash
# in the directory with your input reads
docker run -v $(pwd):/workspace -w /workspace --user $(id -u):$(id -g) vout/megahit \
  megahit -1 pe_1.fq.gz -2 pe_2.fq.gz -o megahit_out
```


### Compile from source
#### Requirement
- zlib
- bzip2
- gzip
- cmake >= 3
- gcc >= 5

```
git clone https://github.com/voutcn/megahit.git
cd megahit
git submodule update --init
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=release ..
make simple_test
./megahit -1 pe_1.fq.gz -2 pe_2.fq.gz -o megahit_out
```

## Introduction
MEGAHIT is a single node metagenome assembler for NGS reads. It makes use of succinct *de Bruijn* graph (SdBG) to achieve low memory assembly.

## What's new in version 1.2
TBD

## Publications
If you use MEGAHIT v0.x or want to cite MEGAHIT for general purpose (e.g. review), please cite:
- Li, D., Liu, C-M., Luo, R., Sadakane, K., and Lam, T-W., (2015) MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics*, doi: 10.1093/bioinformatics/btv033 [PMID: [25609793](http://www.ncbi.nlm.nih.gov/pubmed/25609793)].

If you use MEGAHIT v1.0 or higher version, or assemblies in [MEGABOX](http://hku-bal.github.io/megabox/), please also cite:
- Li, D., Luo, R., Liu, C.M., Leung, C.M., Ting, H.F., Sadakane, K., Yamashita, H. and Lam, T.W., 2016. MEGAHIT v1.0: A Fast and Scalable Metagenome Assembler driven by Advanced Methodologies and Community Practices. Methods.
