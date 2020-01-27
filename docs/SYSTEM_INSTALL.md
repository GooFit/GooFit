# Installing GooFit on different systems

The following commands show you how to get a *minimal* install of GooFit on a vanilla system; you will probably want to add ROOT for your system, and possibly CUDA if you have a graphics card. If you do not have ROOT, some functionality, such as the Minuit1 version of the fitter, will not be available, and most of the examples will not be included in the build. Since ROOT is such a common addition, some of the system also have a note on installing ROOT.

<details><summary>CentOS 7: (click to expand)</summary><p>

For simplicity, this uses EPEL to get access to `python-pip`, and uses the pip version of CMake. Feel free to download CMake directly from Kitware instead. If you want to use Docker, you can start with `docker run --rm -it centos`.
You can also use this recipe with an [nvidia-docker CentOS image](https://hub.docker.com/r/nvidia/cuda/).

```bash
yum install epel-release -y
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
mkdir root-6 && curl https://root.cern.ch/download/root_v6.14.02.Linux-centos7-x86_64-gcc4.8.tar.gz | tar --strip-components=1 -xz -C root-6
source root-6/bin/thisroot.sh
```
</p></details>

<details><summary>Alpine Linux 3.8: (click to expand)</summary><p>

A truly minimal system, Alpine gives you a working Docker system under 3 MB. Since it is unlikely that you'll be running Alpine outside of Docker, the Docker command is included.

```bash
docker run --rm -it alpine
apk add make cmake g++ git libexecinfo-dev
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
docker run --rm -it alpine
apk add build-base cmake git
git clone https://github.com/GooFit/GooFit.git
cd GooFit
make auto
```

If you'd like to use the Python version:

```bash
docker run --rm -it alpine
apk add python3-dev build-base cmake ninja git libexecinfo-dev
pip3 install scikit-build
pip3 -v install git+https://github.com/GooFit/GooFit.git
```

If you want Python2 instead, either add the `py2-pip` package or the line `python2 -m ensurepip`. Python 3 comes with Pip in Alpine. You can also `pip install -v goofit` if you want the latest released PyPI version.

</p></details>

<details><summary>Clang-style checking on Alpine (click to expand)</summary><p>

If you'd like to use LLVM's clang-format to check the style, probably the easiest way to do that is with Docker and Alpine. The following lines will run the style check for you:

```bash
docker run --rm -it alpine
apk add clang git
git clone https://github.com/GooFit/GooFit.git
cd GooFit
./scripts/check_style.sh
```

</p></details>

<details><summary>Ubuntu 18.04 LTS (click to expand)</summary><p>

Ubiquitous Ubuntu works also. Ubuntu was used for the NVidia Docker solution due to better support from NVidia. The following example uses `ninja-build` instead of make, but make works if you prefer it. You should also be able to use this recipe with  `docker run -it ubuntu` or
an [nvidia-docker Ubuntu image](https://hub.docker.com/r/nvidia/cuda/) if you don't have Ubuntu installed.

```bash
apt update && apt install -y git cmake ninja-build g++
git clone --recursive https://github.com/GooFit/GooFit.git
cd GooFit
mkdir build
cd build
cmake ..
cmake --build .
ctest
```

If you want Python bindings, add `python-dev` or `python3-dev` to the list.

If you'd like to add ROOT, add the following lines before running CMake:
```bash
mkdir root-6 && curl https://root.cern.ch/download/root_v6.14.02.Linux-ubuntu18-x86_64-gcc7.3.tar.gz | tar --strip-components=1 -xz -C root-6
source root-6/bin/thisroot.sh
```

</p></details>

<details><summary>OpenSUSE (click to expand)</summary><p>

If you use `make`, adding `-jN` where `N` is the number of cores will make builds much faster on multicore systems! `ninja` does this automatically.
If you'd like to use Docker to provide an OpenSUSE environment, add
`docker run -it opensuse sh`
to the beginning of these commands.

```bash
zypper install -y git cmake gcc-c++
git clone --recursive https://github.com/GooFit/GooFit.git
cd GooFit
make
```

If you'd like to add ROOT:

As always, https://root.cern.ch/build-prerequisites#opensuse is a good resource, though `glu-devel` seems to now be required as well. I'm using `ninja` to build in parallel automatically (and it seems to be faster as well).

```bash
zypper install -y git bash cmake gcc-c++ gcc binutils xorg-x11-libX11-devel xorg-x11-libXpm-devel xorg-x11-devel xorg-x11-proto-devel xorg-x11-libXext-devel glu-devel
zypper install -y ninja tar
git clone --depth=1 --branch=v6-14-02 http://github.com/root-project/root.git root_git
mkdir root_build
cd root_build
cmake ../root_git -GNinja -DCMAKE_INSTALL_PREFIX=/opt/root-6-14-02
cmake --build . --target install
```

Then, you can activated that copy of root with:

```bash
source /opt/root-6-14-02/bin/thisroot.sh
```

You might want to add useful extra ROOT library: `-Droofit=ON -Dmathmore=ON -Dminuit2=ON`

</p></details>

<details><summary>Docker (click to expand)</summary><p>

If you are interested in actually running on Docker, you can use the official GooFit Docker images:

```bash
docker run --rm -it goofit/goofit-omp
nvidia-docker run --rm -it goofit/goofit-cuda
```

The CUDA version will need to build on your computer; the OMP version is prebuilt. You can also start with the ROOT Docker instance:

```
docker run --rm -it rootproject/root-ubuntu16 bash
cd
git clone --recursive https://github.com/GooFit/GooFit.git
cd GooFit
make
```

</p></details>

<details><summary>SLC 6 with CVMFS (LxPlus 6): (click to expand)</summary><p>

```bash
# If you have not run this already (automatic on LxPlus):
source /cvmfs/lhcb.cern.ch/group_login.sh

# Set up LCG releases
. /cvmfs/sft.cern.ch/lcg/views/LCG_94python3/x86_64-slc6-gcc62-opt/setup.sh

# Download and build as usual (ssh download doesn't seem to
# work properly due to bug on LxPlus)
git clone --recursive https://github.com/GooFit/GooFit.git
cd GooFit
mkdir build
cd build
cmake ..
make -j4
```

</p></details>

<details><summary>NVidia (click to expand)</summary><p>

```bash
docker run --rm -v $PWD:/goofit -it nvidia/cuda:9.2-devel-ubuntu18.04 bash
apt update && apt install -y cmake
mkdir build
cd build
cmake ../goofit -DGOOFIT_ARCH=3.5
make -j4
```

Or, to select a newer version of CMake:

```bash
docker run --rm -v $PWD:/goofit -it nvidia/cuda:9.2-devel-ubuntu18.04 bash
apt update && apt install -y wget
wget -qO- "https://cmake.org/files/v3.16/cmake-3.16.3-Linux-x86_64.tar.gz" | tar --strip-components=1 -xz -C /usr/local
cmake -S goofit -B build -DGOOFIT_ARCH=3.5
cmake --build build -j4
```

</p></details>

<details><summary>Note about installing CMake: (click to expand)</summary><p>

While other install methods for CMake, like `pip`, are easier, this way should always work. On Linux, you can manually get a version of CMake using:

```bash
mkdir cmake && wget -qO- "https://cmake.org/files/v3.16/cmake-3.16.3-Linux-x86_64.tar.gz" | tar --strip-components=1 -xz -C cmake
export PATH=`pwd`/cmake/bin:$PATH
```

The second line will need to be rerun whenever use a new shell. Feel free to make your updated CMake default; CMake is insanely backward compatible and will even "dumb itself down" when it sees a lower version in the `minimum_required` line in  `CMakeLists.txt`.

If you are a fan of using `~/.local` and already have `~/.local/bin` in your path, you can instead use:

```bash
wget -qO- "https://cmake.org/files/v3.16/cmake-3.16.3-Linux-x86_64.tar.gz" | tar --strip-components=1 -xz -C ~/.local
```

</p></details>
