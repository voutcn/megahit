MEGAHIT
=======

[![BioConda Install](https://img.shields.io/conda/dn/bioconda/megahit.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/megahit) [![Build Status](https://travis-ci.org/voutcn/megahit.svg?branch=master)](https://travis-ci.org/voutcn/megahit) [![codecov](https://codecov.io/gh/voutcn/megahit/branch/master/graph/badge.svg)](https://codecov.io/gh/voutcn/megahit)

MEGAHIT is an ultra-fast and memory-efficient NGS assembler. It is optimized for metagenomes, but also works well on generic single genome assembly (small or mammalian size) and single-cell assembly.

*News: try v1.2.x!*
------

MEGAHIT v1.2.x (beta) is released. Compared to v1.1.x, its changes include

-   faster and more memory-efficient than before, by using [BMI2 instructions](https://en.wikipedia.org/wiki/Bit_Manipulation_Instruction_Sets), [sparsepp](https://github.com/greg7mdp/sparsepp) and [xxhash](https://github.com/Cyan4973/xxHash)
-   refactored with C++11 features
-   use CMake to build the project
-   removal of GPU support

Please follow the instructions in [Getting Started](#gst) to try this new version.
Past versions can be found at the [release](https://github.com/voutcn/megahit/releases) page.

<a name="gst"></a>Getting Started
---------------

### Running with Linux binaries or docker images (recommended)

``` sh
wget https://github.com/voutcn/megahit/releases/download/v1.2.3-beta/MEGAHIT-1.2.3-beta-Linux-static.tar.gz
tar zvxf MEGAHIT-1.2.3-beta-Linux-static.tar.gz
cd MEGAHIT-1.2.3-beta-Linux-static/bin/
./megahit --test  # run on a toy dataset
./megahit -1 YOUR_PE_READ_1.gz -2 YOUR_PE_READ_2.fq.gz -o YOUR_OUTPUT_DIR
```

You can also run MEGAHIT with its docker images.

``` sh
# in the directory with your input reads
docker run -v $(pwd):/workspace -w /workspace --user $(id -u):$(id -g) vout/megahit \
  megahit -1 YOUR_PE_READ_1.gz -2 YOUR_PE_READ_2.fq.gz -o YOUR_OUTPUT_DIR
```

### Building from source

#### Prerequisites

-   For building: zlib, cmake &gt;= 2.8, g++ &gt;= 4.8.4
-   For running: gzip and bzip2

``` sh
git clone https://github.com/voutcn/megahit.git
cd megahit
git submodule update --init
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release  # add -DCMAKE_INSTALL_PREFIX=YOUR_PREFIX if needed
make -j4
make simple_test  # will test MEGAHIT with a toy dataset
# make install if needed
```

Usage
-----

To run MEGAHIT with default parameters:

``` sh
megahit -1 YOUR_PE_READ_1.fq.gz -2 YOUR_PE_READ_2.fq.gz -r YOUR_SE_READ.fq.gz -o YOUR_OUTPUT_DIR
```

If you did not install Megahit to your PATH, just run Megahit with full-path, e.g.

``` sh
/PATH/TO/MEGAHIT/build/megahit
```

To see the full manual of Megahit, run the program without parameters or with `-h`.

Also, our [wiki](https://github.com/voutcn/megahit/wiki) may be helpful.

Publications
------------

-   Li, D., Liu, C-M., Luo, R., Sadakane, K., and Lam, T-W., (2015) MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics*, doi: 10.1093/bioinformatics/btv033 \[PMID: [25609793](http://www.ncbi.nlm.nih.gov/pubmed/25609793)\].
-   Li, D., Luo, R., Liu, C.M., Leung, C.M., Ting, H.F., Sadakane, K., Yamashita, H. and Lam, T.W., 2016. MEGAHIT v1.0: A Fast and Scalable Metagenome Assembler driven by Advanced Methodologies and Community Practices. Methods.

License
-------

This project is licensed under the GPLv3 License - see the [LICENSE](LICENSE) file for details
