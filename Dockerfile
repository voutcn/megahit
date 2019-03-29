FROM ubuntu:18.04
RUN apt-get update
RUN apt-get install -y g++ make zlib1g-dev gzip bzip2 cmake python --no-install-recommends
COPY . /root/megahit
WORKDIR /root/megahit
RUN rm -rf build
RUN mkdir -p build
WORKDIR build
RUN cmake -DCMAKE_BUILD_TYPE=Release ..
RUN make -j4
RUN make install
RUN megahit --test --no-hw-accel
RUN megahit --test --no-hw-accel --kmin-1pass
