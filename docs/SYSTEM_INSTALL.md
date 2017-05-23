# Installing GooFit on different systems

The docker command for each system is listed at the beginning of each command, to show setup from scratch. Ignore that line if you are not using docker.

The following commands show you how to get a *minimal* install of GooFit on a vanilla system; you will probably want to add ROOT for your system, and possibly CUDA if you have a graphics card. If you do not have ROOT, some functionality, such as the Minuit1 version of the fitter, will not be available, and most of the examples will not be included in the build.

## CentOS 7

For simplicity, this uses the EPEL version of CMake; feel free to download CMake directly from Kitware instead. The EPEL version has the odd `cmake3` name to distinguish it from CentOS default.

```bash
docker run -it centos
yum install epel-release -y
yum install cmake3 git -y
yum group install -y "Development Tools"
git clone --recursive https://github.com/GooFit/GooFit.git
cd GooFit
mkdir build
cd build
cmake3 ..
make
make test
```

## Alpine Linux 3.5

A truly minimal system, Alpine gives you a working Docker system under 5 MB. To get a minimal install of GooFit:

```bash
docker run -it alpine sh
apk add --no-cache make cmake g++ git libexecinfo-dev
git clone --recursive https://github.com/GooFit/GooFit.git
cd GooFit
mkdir build
cd build
cmake ..
make
make test
```

## Ubuntu 17.04

Ubiquitous Ubuntu works also. Ubuntu was used for the NVidia docker solution due to better support from NVidia. 

```bash
docker run -it ubuntu
apt-get update && apt-get install -y git cmake make g++
git clone --recursive https://github.com/GooFit/GooFit.git
cd GooFit
mkdir build
cd build
cmake ..
make
make test
```
