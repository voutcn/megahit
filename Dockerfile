FROM ubuntu:18.04

COPY . /root/megahit
WORKDIR /root/megahit

RUN apt-get update && apt-get install -y --no-install-recommends \
    bzip2 \
    cmake \
    gzip \
    g++ \
    libgomp1 \
    make \
    python \
    zlib1g-dev && \
    rm -rf build && \
    mkdir -p build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release .. && \
    make -j4 install && \
    apt-get autoremove --purge -y \
    cmake \
    g++ \
    make \
    zlib1g-dev

RUN megahit --test && megahit --test --kmin-1pass
ENTRYPOINT ["megahit"]
