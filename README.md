[![Build Status](https://travis-ci.org/GooFit/GooFit.svg?branch=master)](https://travis-ci.org/GooFit/GooFit)
[![Gitter](https://badges.gitter.im/GooFit/GooFit.svg)](https://gitter.im/GooFit/Lobby)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](./LICENSE)

GooFit is a massively-parallel framework, written in CUDA, for
doing maximum-likelihood fits with a comfortable syntax.
It is also possible to build
GooFit using OpenMP.

* [Changelog](./CHANGELOG.md)
* [Contributing](./CONTRIBUTING.md)
* [API documentation](https://GooFit.github.io/GooFit)

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
> Now, just use the local copy instead:
>
> ```bash
> ./cmake/bin/cmake ..
> ```

If you want to change compiler, set `CC` and `CXX` to appropriate defaults *before* you run cmake either inline or in your environment. If you want to set the host and device backends, you can set those options. The defaults are:
```
cmake .. -DGOOFIT_DEVICE=CUDA -DGOOFIT_HOST=OMP
```

Valid options are `CUDA` (device only), `OMP`, and CPP (host only). While the current backend is used, `CPP`, and `TBB` will probably remain unavailable for device calculations. The default is `Auto`, and will select `CUDA` if CUDA is found.

Other custom options supported along with the defaults:

* `-DGOOFIT_ARCH=Auto`: (Auto, Common, All, valid number(s) or name(s)): sets the compute architecture. See [CUDA_SELECT_NVCC_ARCH_FLAGS](https://cmake.org/cmake/help/v3.7/module/FindCUDA.html).
* `-DGOOFIT_EXAMPLES=ON`: Build the examples
* `-DGOOFIT_PACKAGES=ON`: Build any packages found with the name `goofit*`
* `-DGOOFIT_SEPARATE_COMP=OFF`: Enable separable compilation of PDFs
* `-DGOOFIT_TESTS=OFF`: Build the goofit tests

Advanced Options:
* `-DGOOFIT_MPI=ON`: (OFF/ON.  With this feature on, GPU devices are selected automatically)
* `-DGOOFIT_CUDA_OR_GROUPSIZE:INT=128`: This sets the group size that thrust will use for distributing the problem.  This parameter can be thought of as 'Threads per block'.  These will be used after running 'find_optimal.py' to figure out the optimal size.
* `-DGOOFIT_CUDA_OR_GRAINSIZE:INT=7`: This is the grain size thrust uses for distributing the problem.  This parameter can be thought of as 'Items per thread'.


A few standard cmake tricks:

* Use `VERBOSE=1` to see the commands used to build the files
* Use `-L` to list the CMake options
* Use `ccmake` if available to see a curses (terminal) gui, or `cmake-gui` for a completely graphical interface
* Use `-G` and the name of a generator to use something other than `make`, like `Xcode` or `Ninja`
* Open the `CMakeLists.txt` with QtCreator to generate for that IDE
* Set the release type with `-DCMAKE_BUILD_TYPE=Release`, `RelWithDebInfo`, `Debug`, etc.
* Set up multiple build directories, like `build-omp`
* CMake remembers your `-D` option selections in your build directory so you don't have to specify them again
* CMake reruns when needed when you `make`
* Use `-j12` to build with 12 cores (for example)
* Use `cmake --build .` to build without refering to your specific build tool, like `make`


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

(This requires the [Plumbum](https://plumbum.readthedocs.io/en/latest/) library, install with `pip install plumbum`, `pip install --user plumbum`, or `conda -c conda-forge plumbum`.)

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

To add packages, use standard CMake tools. For example (CMake 3.5), to add [Boost](https://cmake.org/cmake/help/v3.7/module/FindBoost.html) 1.49+ filesystem and `TTreeReader` from ROOT:

```cmake
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.49 REQUIRED COMPONENTS filesystem)

goofit_add_executable(K3Pi K3Pi.cu)
target_link_libraries(MyNewExample Boost::filesystem ROOT::TreePlayer)

```

## Adding a new project
  
If you'd like to make a seperate goofit project, you can do so. Simply checkout your project inside goofit, with the name `work` or `goofit`+something. CMake will automatically pick up those directories and build them, and GooFit's git will ignore them. Otherwise, they act just like the example directory. If you add a new directory, you will need to explicitly rerun cmake, as that cannot be picked up by the makefile. The automatic search can be turned off with the `GOOFIT_PROJECTS` option.

## Converting from older GooFit code
 
The build system underwent a major upgrade in the move to CMake. The folders that were introduced to keep the includes structured require modifications of source code, converting lines like `#include "Variable.hh"` to `#include "goofit/Variable.h"`. This modification can be done for you by running the provided script, `scripts/ModernizeGooFit.py` on your source files (requires Python and Plumbum). You should remove your old Makefiles and use the new `CMakeFiles.txt` files provided in examples - this should require
writing two lines of code instead of the 50 or so previously needed.

## Improving Performance with MPI

Using the MPI verion with an appropriate environment setup will allow for multiple GPU's to be used, or allow for multiple nodes.  To use this feature simpley turn the flag on with cmake `-DGOOFIT_MPI=ON`.  This will divide the dataset by the number of processes involved.  For instance, if you have two nodes that will be involved in the calculation, the data will be split in half.  Currently, each node will load the entire buffer from disk, then load partitioned data it will work on.  

A few notes about using the MPI version:
* You will need to use the `CountingVariable` for any event numbers used or referenced within the code, or anything that counts with the events.
* Please call setDataSize after setData.  If you do not, setDataSize doesn't have `m_iEventsPerTask`, which will need to be recalculated. 

## Configuring Group Size & Grain Size

This advanced option is for GPU devices only.  What it will do is search the group and grain space in order to find the optimal configuration for the particular PDFs.  This should be run after an example has been developed and tested.  Please look at scripts/find_optimal.py to see how to formulate a particular script.  Depending on the searchable space, this can take hours to days to compute.  
The script will loop over the space and configure each parameter, then recompile and run the example a number of times.  The average is provided in the spreadsheet, and the fastest version is printed to the user.

## Acknowledgement

GooFit's development has been supported by the National Science Foundation under grant number NSF-1005530. 
Any opinions, findings, and conclusions or recommendations expressed in this material are those of the developers
and do not necessarily reflect the views of the National Science Foundation.
