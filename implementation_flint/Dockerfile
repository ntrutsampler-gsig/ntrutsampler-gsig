# syntax=docker/dockerfile:1

FROM ubuntu:jammy-20260509

# Copy code and libraries to docker image
RUN mkdir -p /home/ntru_tsampler
COPY code /home/ntru_tsampler

# Install compiler and required build tools
RUN apt-get update
RUN apt-get install -y gcc g++ m4 autoconf autotools-dev make libtool cmake xz-utils

# download and unpack GMP, MPFR and FLINT
WORKDIR /home/ntru_tsampler/libs
ADD https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz /home/ntru_tsampler/libs/gmp.tar.xz
ADD https://www.mpfr.org/mpfr-current/mpfr-4.2.2.tar.xz /home/ntru_tsampler/libs/mpfr.tar.xz
ADD https://flintlib.org/download/flint-3.5.0.tar.gz /home/ntru_tsampler/libs/flint.tar.gz
RUN tar -xf gmp.tar.xz
RUN tar -xf mpfr.tar.xz
RUN tar -xf flint.tar.gz

# Setup GMP library
WORKDIR /home/ntru_tsampler/libs/gmp-6.3.0
RUN ./configure
RUN make -j
RUN make check
RUN make install

# Setup MPFR library
WORKDIR /home/ntru_tsampler/libs/mpfr-4.2.2
RUN ./configure
RUN make -j
RUN make check
RUN make install

# Setup FLINT library
WORKDIR /home/ntru_tsampler/libs/flint-3.5.0
RUN ./bootstrap.sh
RUN ./configure --enable-avx2 --disable-pthread CFLAGS="-O3 -Wall -march=native" CC="gcc"
RUN make -j
RUN make check
RUN make install

# Build ntru_tsampler
RUN mkdir -p /home/ntru_tsampler/_build
WORKDIR /home/ntru_tsampler/_build
RUN cmake ..
RUN make -j