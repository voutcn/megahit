FROM ubuntu:18.04
MAINTAINER Dinghua Li (voutcn@gmail.com)
RUN apt-get update
RUN apt-get install -y g++ zlib1g-dev gzip bzip2 cmake
RUN apt-get install -y python
ADD . /root/megahit
RUN cd /root/megahit && mkdir -p build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && \
    make -j4 simple_test && make install
