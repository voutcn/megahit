# MEGAHIT
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/megahit.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/megahit)
[![Build Status](https://travis-ci.org/voutcn/megahit.svg?branch=master)](https://travis-ci.org/voutcn/megahit)

MEGAHIT is an ultra-fast and memory-efficient NGS assembler. It is optimized for metagenomes, but also works well generic single genome assembly (small or mammalian size) and single-cell sequencing assembly.

## News
MEGAHIT v1.2.0-beta is released. Main changes include
- faster and more memory-efficient than before, by using [BMI2 instructions](https://en.wikipedia.org/wiki/Bit_Manipulation_Instruction_Sets), [sparsepp](https://github.com/greg7mdp/sparsepp) and [xxhash](https://github.com/Cyan4973/xxHash).
- refactored with C++11 features
- removal of GPU support

It is highly recommended to use v1.2.0-beta. Past versions can be found at the [release](https://github.com/voutcn/megahit/releases) page.

## Getting Started

### Run with docker (recommended)
```bash
# in the directory with your input reads
docker run -v $(pwd):/workspace -w /workspace --user $(id -u):$(id -g) vout/megahit \
  megahit -1 YOUR_PE_READ_1.gz -2 YOUR_PE_READ_2.fq.gz -o YOUR_OUTPUT_DIR
```

### Build from source
#### Prerequisites
- For building: zlib, cmake, gcc/g++ >= 5
- For running: gzip and bzip2

#### Build and test

1. Obtain the source code
```bash
git clone https://github.com/voutcn/megahit.git
cd megahit
git submodule update --init
```

2. Create the build directory
```bash
mkdir build && cd build
```
3. Run cmake
```bash
cmake -DCMAKE_BUILD_TYPE=release ..
```
If your CPU does not support BMI2 instructions (uncommon), run the following command instead
```
cmake -DUSE_BMI2=OFF -DCMAKE_BUILD_TYPE=release ..
```
4. Compile & test
```bash
make -j4
make simple_test  # will test MEGAHIT with a toy dataset
```
If you need to install Megahit to your PATH, run `make install` in the build directory.

## Usage

To run MEGAHIT with default parameters:
```bash
megahit -1 YOUR_PE_READ_1.fq.gz -2 YOUR_PE_READ_2.fq.gz -r YOUR_SE_READ.fq.gz -o YOUR_OUTPUT_DIR
```

If you did not install Megahit to your PATH, just run Megahit with full-path, e.g.
```bash
/PATH/TO/MEGAHIT/build/megahit
```

To see the full manual of Megahit, run the program without parameters or with `-h`.

Also, our [wiki](https://github.com/voutcn/megahit/wiki) may be helpful.

## Publications
- Li, D., Liu, C-M., Luo, R., Sadakane, K., and Lam, T-W., (2015) MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics*, doi: 10.1093/bioinformatics/btv033 [PMID: [25609793](http://www.ncbi.nlm.nih.gov/pubmed/25609793)].
- Li, D., Luo, R., Liu, C.M., Leung, C.M., Ting, H.F., Sadakane, K., Yamashita, H. and Lam, T.W., 2016. MEGAHIT v1.0: A Fast and Scalable Metagenome Assembler driven by Advanced Methodologies and Community Practices. Methods.

## License
This project is licensed under the GPLv3 License - see the [LICENSE.md](LICENSE.md) file for details
