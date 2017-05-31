# Installing GooFit on different systems

The Docker command for each system is listed at the beginning of each command, to show setup from scratch. Ignore that line if you are not using Docker.

The following commands show you how to get a *minimal* install of GooFit on a vanilla system; you will probably want to add ROOT for your system, and possibly CUDA if you have a graphics card. If you do not have ROOT, some functionality, such as the Minuit1 version of the fitter, will not be available, and most of the examples will not be included in the build.

## CentOS 7

For simplicity, this uses EPEL to get access to `python-pip`, and uses the pip version of CMake. Feel free to download CMake directly from Kitware instead.

```bash
docker run -it centos
yum install epel-release
yum install python-pip git gcc-c++ make -y
pip install cmake plumbum
git clone --recursive https://github.com/GooFit/GooFit.git
cd GooFit
mkdir build
cd build
cmake ..
make
make test
```

If you'd like to add ROOT, add the following lines before running CMake:

```bash
mkdir root-6 && curl https://root.cern.ch/download/root_v6.08.06.Linux-centos7-x86_64-gcc4.8.tar.gz | tar --strip-components=1 -xz -C root-6
source root-6/bin/thisroot.sh
```

## Alpine Linux 3.5

A truly minimal system, Alpine gives you a working Docker system under 3 MB.

```bash
docker run -it alpine
apk add --no-cache make cmake g++ git libexecinfo-dev
git clone --recursive https://github.com/GooFit/GooFit.git
cd GooFit
mkdir build
cd build
cmake ..
make
ctest
```

In the spirit of minimality, this is less instructive and contains more magic, but also would also work:

```bash
docker run -it alpine
apk add --no-cache make cmake g++ git
git clone https://github.com/GooFit/GooFit.git
cd GooFit
make
```

## Ubuntu 16.04

Ubiquitous Ubuntu works also. Ubuntu was used for the NVidia docker solution due to better support from NVidia. The following example uses ninja-build instead of make.

```bash
docker run -it ubuntu
apt-get update && apt-get install -y git cmake ninja-build g++
git clone --recursive https://github.com/GooFit/GooFit.git
cd GooFit
mkdir build
cd build
cmake ..
cmake --build .
ctest
```

If you'd like to add ROOT, add the following lines before running cmake:
```bash
mkdir root-6 && curl https://root.cern.ch/download/root_v6.08.06.Linux-ubuntu16-x86_64-gcc5.4.tar.gz | tar --strip-components=1 -xz -C root-6
source root-6/bin/thisroot.sh
```
