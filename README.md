[![Build Status][travis-badge]][travis-link]
[![Code Coverage][codecov-badge]][codecov-link]
[![Join the chat at https://gitter.im/GooFit/Lobby][gitter-badge]][gitter-link]
[![License: LGPL v3][license-badge]](./LICENSE)

![GooFit logo][./docs/GooFitLogo.png]

GooFit is a massively-parallel framework, written in CUDA, for
doing maximum-likelihood fits with a comfortable syntax.
It is also possible to build
GooFit using OpenMP.

* [Changelog](./CHANGELOG.md)
* [Contributing](./CONTRIBUTING.md)
* [API documentation]

## Requirements

* A new version of CMake if using the CMake build. Like ROOT, the minimum is 3.4, but tested primarily with 3.7. CMake is incredibly easy to install.
  * With CMake, Thrust is downloaded automatically for OpenMP if not found
  * GoogleTest is downloaded automatically
* A ROOT 6 build highly recommended.
* If using CUDA:
  * CUDA 7.5 and 8.0 tested, older versions should work for now
  * An nVidia GPU supporting compute capability at least 2.0 (3.5 recommended)
* If using OpenMP:
  * A compiler supporting OpenMP and C++11 (GCC 4.8+ and Intel 17 tested)

## Getting the files

* Clone with git:

```bash
git clone git://github.com/goofit/goofit.git
cd goofit
```

You can either checkout a tagged version, or stay on the master for the latest and greatest. There are often development branches available, too.

## Building 

The build system now uses CMake. The procedure is standard for CMake builds:

```bash
mkdir build
cd build
cmake ..
make
```

> If you don't have a modern CMake, you can get a local copy on linux using:
>
> ```bash
> mkdir cmake && wget -qO- "https://cmake.org/files/v3.7/cmake-3.7.2-Linux-x86_64.tar.gz" | tar --strip-components=1 -xz -C cmake
> ```
>
> Now, before running cmake just use the local copy instead:
>
> ```bash
> export PATH=`pwd`/cmake/bin:$PATH
> ```

If you want to change compiler, set `CC` and `CXX` to appropriate defaults *before* you run cmake either inline or in your environment. If you want to set the host and device backends, you can set those options. The defaults are:
```
cmake .. -DGOOFIT_DEVICE=CUDA -DGOOFIT_HOST=OMP
```

Valid options are `CUDA` (device only), `OMP`, and CPP (host only). While the current backend is used, `CPP`, and `TBB` will probably remain unavailable for device calculations. The default is `Auto`, and will select `CUDA` if CUDA is found.

Other custom options supported along with the defaults:

* `-DGOOFIT_ARCH=Auto` (Auto, Common, All, valid number(s) or name(s)): sets the compute architecture. See [CUDA_SELECT_NVCC_ARCH_FLAGS].
* `-DGOOFIT_EXAMPLES=ON`: Build the examples
* `-DGOOFIT_PACKAGES=ON`: Build any packages found with the name `goofit*`
* `-DGOOFIT_SEPARATE_COMP=OFF`: Enable separable compilation of PDFs
* `-DGOOFIT_TESTS=OFF`: Build the goofit tests


A few standard cmake tricks:

* Use `VERBOSE=1` to see the commands used to build the files.
* Use `-L` to list the CMake options.
* Use `ccmake` if available to see a curses (terminal) gui, or `cmake-gui` for a completely graphical interface.
* Use `-G` and the name of a generator to use something other than `make`, like `Xcode` or `Ninja`.
* Open the `CMakeLists.txt` with QtCreator to generate for that IDE.
* Set the release type with `-DCMAKE_BUILD_TYPE=Release`, `RelWithDebInfo`, `Debug`, etc.
* Set up multiple build directories, like `build-omp` and `build-cuda`.
* CMake remembers your `-D` option selections in your build directory so you don't have to specify them again.
* CMake reruns when needed when you `make` unless you add a file that it globs for (like new `goofit_projects`).
* Use `-j12` to build with 12 cores (for example).
* Use `cmake --build .` to build without referring to your specific build tool, like `make` or `ninja`.


> ## Classic Makefile system (depreciated)
>   
> * You should set `CUDALOCATION` for your system
> 
> * You should have run source `thisroot.sh` to setup ROOT paths and other environment variables
> 
> * Set `TARGET_OMP=1` if you want to use OMP
>   * Checkout a copy of CUDA's thrust next to the goofit repository (there's nothing to compile)
>   * If you already have thrust in a different location, set `THRUSTLOCATION`
> 
> ```
> make all TARGET_OMP=1
> ```
> 
> Following the fairly standard makefile convention, all programs are built in-place instead of in a build directory.

## Running the Examples

* To run all the examples, with timing information, use:
```
./examples/RunAll.py
```

(This requires the [Plumbum] library, install with `pip install plumbum`, `pip install --user plumbum`, or `conda -c conda-forge plumbum`.)

If you want to run an individual example, those are in subdirectories in examples (built products are in your build directory, the source is in `goofit/examples`).


## Adding a new example:

The examples are designed to be easy to add to. Make a new directory, then add a new CMakeLists.txt in your directory with one or more of the following two lines:

```cmake
goofit_add_directory()
goofit_add_executible(MyNewExample MyNewExample.cu)
```

The first line adds your `.cu` file with goofit code as an executible, and the second one sets up a symbolic links to the source and datafiles in the build directory to the source directory. If you perfer to only have some files symbolically linked, use `goofit_add_link(filename.ext)` explicitly for each file. To get the example to build when you build goofit, add the name of your directory to `examples/CMakeLists.txt`.

If you are building with separable compilation, you can also use `goofit_add_pdf(mypdf.cu)` to add a PDF. This will also require that you include any directory that you need with `include_directory`, as usual.

> If you want to extend the Makefile system instead, copy a Makefile from a different directory, changing the relevent project name (only one program per directory supported), and make a new target in `examples/Makefile`. 

To add packages, use standard CMake tools. For example (CMake 3.5), to add [Boost][FindBoost] 1.49+ filesystem and `TTreeReader` from ROOT:

```cmake
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.49 REQUIRED COMPONENTS filesystem)

goofit_add_executable(K3Pi K3Pi.cu)
target_link_libraries(MyNewExample Boost::filesystem ROOT::TreePlayer)

```

## Adding a new project
  
If you'd like to make a separate goofit project, you can do so. Simply checkout your project inside goofit, with the name `work` or `goofit`+something. CMake will automatically pick up those directories and build them, and GooFit's git will ignore them. Otherwise, they act just like the example directory. If you add a new directory, you will need to explicitly rerun cmake, as that cannot be picked up by the makefile. The automatic search can be turned off with the `GOOFIT_PROJECTS` option.

## Converting from older GooFit code
 
The build system underwent a major upgrade in the move to CMake. The folders that were introduced to keep the includes structured require modifications of source code, converting lines like `#include "Variable.hh"` to `#include "goofit/Variable.h"`. This modification can be done for you by running the provided script, `scripts/ModernizeGooFit.py` on your source files (requires Python and Plumbum). You should remove your old Makefiles and use the new `CMakeFiles.txt` files provided in examples - this should require
writing two lines of code instead of the 50 or so previously needed.

## Acknowledgement

GooFit's development is supported by the National Science Foundation under grant number [1414736]
and was developed under grant number [1005530]. 
Any opinions, findings, and conclusions or recommendations expressed in this material are those of the developers
and do not necessarily reflect the views of the National Science Foundation.

[API documentation]: https://GooFit.github.io/GooFit
[travis-badge]:      https://travis-ci.org/GooFit/GooFit.svg?branch=master  
[travis-link]:       https://travis-ci.org/GooFit/GooFit
[codecov-badge]:     https://codecov.io/gh/GooFit/GooFit/branch/master/graph/badge.svg
[codecov-link]:      https://codecov.io/gh/GooFit/GooFit
[gitter-badge]:      https://badges.gitter.im/GooFit/GooFit.svg
[gitter-link]:       https://gitter.im/GooFit/Lobby
[license-badge]:     https://img.shields.io/badge/License-LGPL%20v3-blue.svg
[1005530]:           https://nsf.gov/awardsearch/showAward?AWD_ID=1005530
[1414736]:           https://nsf.gov/awardsearch/showAward?AWD_ID=1414736
[CUDA_SELECT_NVCC_ARCH_FLAGS]: https://cmake.org/cmake/help/v3.7/module/FindCUDA.html
[Plumbum]:           https://plumbum.readthedocs.io/en/latest/
[FindBoost]:         https://cmake.org/cmake/help/v3.7/module/FindBoost.html
