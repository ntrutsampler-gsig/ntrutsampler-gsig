# Practical Group Signatures from Tag-Based NTRU Sampler
---

This repository contains the code for the paper **Practical Group Signatures from Tag-Based NTRU Sampler**. Part of the code was taken from the implementation of [tag-friendly sampler](https://github.com/tagfriendlysampler/tag-friendly_sampler).  

The repository contains
- the parameter selection and security estimation scripts in python for the group signature of Section 5, and
- the two implementations in C of NTRU-TSampler

## Parameter Scripts

In the `/scripts` subdirectory, we provide a parameter estimation script which depends on the lattice-estimator.
The lattice-estimator is included as submodule, so please do a recursive pull.
```shell
git submodule update --init --recursive
```

You can then run the python script `parameters_group_signature.py` which gives the detailed parameter set for the group signature of Section 5. The scripts estimate all the security assumptions using the lattice estimator (without rough flag) which can take around 3 minutes on a standard laptop. 


## C Implementations

The repository contains two full C implementations of our sampler (NTRU-TSampler) with the parameters of our group signature. Each implementation can be built with Docker (see instructions below) but we also provide all the manual build instructions here. Below, the directory `implementationX` should be replaced by `implementation_flint` (first implementation referenced in the paper) or `implementation` (second implementation referenced in the paper).

Requirements:
- FLINT, which depends on GMP and MPFR
- AES-NI instructions
- cmake
- gcc >= 13

For library integrity purposes, we do not include the library dependencies in this repository and instead give the instructions to download and build them. The implementation depends on Flint, which itself depends on GMP and MPFR. To install these libraries, run the following commands (from the root directory of this repository). It only needs to be installed once.

```shell
cd implementationX/code
mkdir libs
cd libs
curl --output gmp.tar.xz "https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz"
curl --output mpfr.tar.xz "https://www.mpfr.org/mpfr-current/mpfr-4.2.2.tar.xz"
curl --output flint.tar.gz "https://flintlib.org/download/flint-3.5.0.tar.gz"
tar -xf gmp.tar.xz
tar -xf mpfr.tar.xz
tar -xf flint.tar.xz

cd gmp-6.3.0
./configure
make -j
make check
make install

cd ../mpfr-4.2.2
./configure
make -j
make check
make install

cd ../flint-3.5.0
./bootstrap.sh 
./configure --enable-avx2 --disable-pthread CFLAGS="-O3 -Wall -march=native" CC="gcc-13"
make -j
make check
make install
```
To build each implementation, run the following commands (from the root directory of this repository).
```shell
cd implementationX/code
mkdir _build
cd _build
cmake ..
make -j
```

Then, running `./test` or `./bench` runs the tests or the benchmarks, respectively. Each change of the source code requires to run the last command (`make` or `make -j`).

## Building with Docker

We provide a Dockerfile to build a docker image corresponding to our implementations, and give the instructions to run our implementations with Docker here. The same instructions apply to each implementation, by changing `implementationX` with either `implementation_flint` or `implementation`.
For information on how to install docker on your system visit the [docker docs](https://docs.docker.com/).
Depending on your setup, you may need to prefix the following commands with `sudo`.

### 1. Build the docker image
Navigate to the top-level directory which contains this README file. Then run (notice the trailing dot for the second command!)
```shell
cd implementationX
docker build -t ntru_tsampler-docker .
```
This takes some time as docker needs to pull the base-image for Ubuntu and then builds all dependencies
(i.e. GMP, MPFR and FLINT) as well as the test and benchmark executables from scratch.
In particular the `make check` instructions are *very* time consuming.

The libraries are built with standard options and fine-tuning is possible by adjusting the corresponding
`make` and `configure` commands in the Dockerfile. However, installation paths *must not* be changed.

Note that this step needs to be repeated for any changes in the code base or Dockerfile but not
for subsequent runs of the code.

### 2. Run the docker image
Issue the following command to start the docker container, where `YYY` is either `./test` or `./bench` to select
the tests or benchmarks respectively.
```shell
docker run -t ntru_tsampler-docker YYY
```
You should now see the output of either the test or the benchmarking executables.

### 3. Stopping the running docker image
To get a list of running docker containers use `docker ps`. Then you can use `docker stop <container_id>` to stop
the execution. For more information see the [docker docs](https://docs.docker.com/engine/reference/builder/).

### 4. Cleaning up
To remove all stopped docker containers use `docker container prune`.
To remove docker images use `docker images` to get a list of all images and then `docker rm <image_id>` to remove
the image in question.