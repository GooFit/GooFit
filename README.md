[![Actions Status][actions-badge]][actions-link]
[![Travis Status][travis-badge]][travis-link]
[![Join the chat at https://gitter.im/GooFit/Lobby][gitter-badge]][gitter-link]
[![License: BSD][license-badge]](./LICENSE)
[![Latest release][releases-badge]][releases-link]
[![PyPI Status][pypi-status]][pypi-link]
[![Conda-Forge Status][cf-status]][cf-link]
[![DOI][DOI-badge]][DOI-link]
[![Scikit-HEP][sk-badge]](https://scikit-hep.org/)

<p align="center">
  <img width="495" height="220" src="https://raw.githubusercontent.com/GooFit/GooFit/master/docs/GooFitLogo.png"/>
</p>

GooFit is a massively-parallel framework, written using Thrust for CUDA and OpenMP, for
doing maximum-likelihood fits with a familiar syntax.

[What's new](./docs/CHANGELOG.md)
• [Tutorials]
• [API documentation]
• [2.0 upgrade](./docs/CONVERTING20.md)
• [2.1 upgrade](./docs/CONVERTING21.md)
• [2.2 upgrade](./docs/CONVERTING22.md)
• [Build recipes](./docs/SYSTEM_INSTALL.md)
• [Python](https://pypi.python.org/pypi/goofit/)

## Known issues
https://github.com/GooFit/GooFit/labels/critical https://github.com/GooFit/GooFit/labels/amplitude%20analysis https://github.com/GooFit/GooFit/labels/cuda

## Requirements

* A recent version of CMake is required. The minimum is 3.9. CMake is incredibly easy to install, you can even use `pip` (see [the system install page](./docs/SYSTEM_INSTALL.md)). GooFit developers have supplied patches to CMake 3.12, so at least that is highly recommended. CMake 3.16 does not currently work with the Python bindings.
* A ROOT 6 build highly recommended -- GooFit will use the included Minuit2 submodule if ROOT is not found, and the Minuit1 based fitter will not be available. Supports 6.04-6.24 (6.10+ recommended).

<details><summary>If using CUDA: (click to expand)</summary><p>

* CMake 3.9+
* CUDA 8.0+ (with caveats below)
    * CUDA 8: Supported
    * CUDA 9.2, 10.0: Some warnings from Eigen, supported
    * CUDA 9.0, 9.1: Buggy, see [known issues](https://github.com/GooFit/GooFit/issues/173)
    * CUDA 10.1, 10.2: Not yet supported due to Thrust 1.8 incompatibility
* An nVidia GPU supporting compute capability at least 3.0 (3.5+ recommended)

</p></details>

<details><summary>If using OpenMP: (click to expand)</summary><p>

* A compiler supporting OpenMP and C++11 (GCC 4.8+, Clang, and Intel 17 tested, GCC 4.7 not supported)
* Note that TBB is also available as a backend, but it still requires OpenMP to be present.
* On macOS, this backend requires `brew install libomp` or a custom compiler.

</p></details>

<details><summary>If using CPP: (click to expand)</summary><p>

* Single threaded builds are available for debugging and development (such as on the default Clang on macOS)

</p></details>

<br/>

A list of exact commands required for several platforms is [available here](./docs/SYSTEM_INSTALL.md).


<details><summary>Python Bindings: (click to expand)</summary><p>

There are also Python Bindings. This requires Python (2 or 3), [NumPy](http://www.numpy.org), [SciKit-Build](http://scikit-build.readthedocs.io), and CMake. CUDA 8+ is required if using CUDA. If you want the most recent stable release, use `pip install -v goofit` (If you have pip 9 or less, you'll need scikit-build and cmake installed beforehand). Python 2 will be removed soon.

Repository method:

You can uses `pip install -v .` inside the repository. You can also directly force the bindings from a normal build with `-DGOOFIT_PYTHON=ON`. You can check your install with `python -m goofit`. You can debug a goofit file named `python_script.py` with gcc using `gdb -ex r --args python python_script.py`.

Other python requirements for the examples:

* numpy-1.11.1+
* pandas-0.15.1+
* uncertainties-3.0.2
* matplotlib
* plumbum

Optional:

* numba

</p></details>

<br/>

## Getting the files

* Clone with git:

```bash
git clone git://github.com/GooFit/GooFit.git --recursive
cd GooFit
```

You can either checkout a tagged version, or stay on the master for the latest and greatest. There are often development branches available, too. You can use `--jobs=N` or set git's `submodule.fetchJobs` configuration parameter to download the submodules in parallel with `N` threads.

## Building

If you just want to get started as fast as possible, running `make`, `make omp`, or `make cuda` in the main directory will make a build directory for you, and will run CMake and make. It is recommended that you instead directly use the CMake powered build system as described below, so that you will have a better understanding of what you are doing and more flexibility.

The build system uses CMake. The procedure is standard for CMake builds:

```bash
# Classic method
mkdir build
cd build
cmake ..
make -j4 # 4 threads, adjust as needed

# Newer method (CMake 3.13+)
cmake -S . -B build
cmake --build build -j4 # 4 threads, adjust as needed
```

If you don't have a modern CMake, Kitware provides installers for every OS. You can even get a copy using python: `pip install cmake` or locally with `pip install --user cmake`.
On a Mac, you can also use any package manager, such as Homebrew: `brew install cmake`.

If you want to change compiler, set `CC` and `CXX` to appropriate defaults *before* you run CMake either inline or in your environment. You can also set `CMAKE_C_COMPILER` and `CMAKE_CXX_COMPILER` directly on the command line with `-D`. If you want to set the host and device backends, you can set those options. The defaults are:

```bash
cmake .. -DGOOFIT_DEVICE=CUDA -DGOOFIT_HOST=CPP
```

Valid options are `CUDA` (device only), `OMP`, `TBB`, and `CPP`. The Thrust `TBB` backend requires the Intel compiler.  The default device is `Auto`, and will select `CUDA` if CUDA is found, `OMP` or `CPP` otherwise.

Other custom options supported along with the defaults:

* `-DGOOFIT_DEVICE=Auto`: The device to use for computation (`CUDA`, `OMP`, `TBB`, or `CPP`). Default setting of `Auto` looks for CUDA first, then OpenMP, then CPP.
* `-DGOOFIT_ARCH=Auto`: (`Auto`, `Common`, `All`, valid number(s) or name(s)): sets the compute architecture. See [CUDA_SELECT_NVCC_ARCH_FLAGS][]. Can be set to `OFF` to avoid adding any flags.
* `-DGOOFIT_EXAMPLES=ON`: Build the examples
* `-DGOOFIT_PACKAGES=ON`: Build any packages found with the name `goofit_*`
* `-DGOOFIT_DEBUG=ON` and `-DGOOFIT_TRACE=ON` will enable the matching printout macros
* `-DGOOFIT_PYTHON=ON`: Include the python bindings using [pybind11] if Python found (use `-DPYTHON_EXECUTABLE=$(which python3)` to use a specific interpreter).

<details><summary>Advanced Options: (click to expand)</summary><p>

* `-DGOOFIT_HOST=Auto`: This is CPP unless device is `OMP`, in which case it is also `OMP`. This changes `thrust::host_vector` calculations, and is not fully supported when set to a non-default setting.
* `-DGOOFIT_TESTS=ON`: Build the GooFit tests
* `-DGOOFIT_MPI=ON`: (OFF/ON.  With this feature on, GPU devices are selected automatically).  Tested with MVAPICH2/2.2 and OpenMPI.
* You can enable sanitizers on non-CUDA builds with `-DSANITIZE_ADDRESS=ON`, `-DSANITIZE_MEMORY=ON`, `-DSANITIZE_THREAD=ON` or `-DSANITIZE_UNDEFINED=ON`.
* If `clang-tidy` is available, it will automatically be used to check the source. If you set `-DGOOFIT_TIDY_FIX=ON`, fixes will be applied to the GooFit source.
* `-DGOOFIT_SPLASH=ON`: Controls the unicode splash at the beginning.
* `-DGOOFIT_CERNROOT=ON`: Allows you to disable the automatic search for ROOT (used by the PIP Python build)
* `-DCMAKE_UNITY_BUILD=OFF`: Turn on Unity builds in CMake 3.16+. Should be a bit faster (does not speed up CUDA portions of builds).

</p></details>

<details><summary>A few standard CMake tricks: (click to expand)</summary><p>

* Use `make VERBOSE=1` to see the commands used to build the files.
* Use `cmake .. -LH` to list the CMake options with help.
* Use `ccmake` if available to see a curses (terminal) gui, or `cmake-gui` for a completely graphical interface.
* Use `-G` and the name of a generator to use something other than `make`, like `Xcode` or `Ninja`.
* Open the `CMakeLists.txt` with QtCreator to generate for that IDE.
* Set the release type with `-DCMAKE_BUILD_TYPE=Release`, `RelWithDebInfo`, `Debug`, etc.
* Set up multiple build directories, like `build-omp` and `build-cuda`.
* CMake caches your `-D` option selections in your build directory so you don't have to specify them again.
* CMake reruns when needed when you `make` unless you add a file that it globs for (like new `goofit_projects`).
* Use `make -j12` to build with 12 cores (for example). You can set this as the `MAKEFLAGS` environment variable, too.
* Use `CMake --build .` to build without referring to your specific build tool, like `make` or `ninja`.
* If you are using the `llvm` tool-suite, you can use `-DCMAKE_EXPORT_COMPILE_COMMANDS=ON` to generate the .json file that the `clang-*` commands expect.

</p></details>



## Running the examples and tests

* To run all the examples, with timing information, use:

```bash
./examples/RunAll.py
./pyexamples/RunAll.sh # Python
```

(This requires the [Plumbum] library, install with `pip install plumbum`, `pip install --user plumbum`, or `conda -c conda-forge plumbum`.)

If you want to run an individual example, those are in subdirectories in examples (built products are in your build directory, the source is in `/examples`).

The tests can be run with `make test` or `ctest`. The python bindings, if built, can be tested with `pytest`, run from the main build directory. The python examples and tests folders are linked to the build directory with a `py` prefix.

## Other topics

<details><summary>Adding a new example: (click to expand)</summary><p>

The examples are designed to be easy to add to. Make a new directory, then add a new CMakeLists.txt in your directory with one or more of the following two lines:

```cmake
goofit_add_directory()
goofit_add_executible(MyNewExample MyNewExample.cu)
```

The first line adds your `.cu` file with GooFit code as an executable, and the second one sets up a symbolic links to the source and `dataFiles` in the build directory to the source directory. If you prefer to only have some files symbolically linked, use `goofit_add_link(filename.ext)` explicitly for each file. This happens at configure time. To get the example to build when you build GooFit, add the name of your directory to `examples/CMakeLists.txt`.

If you are building with separable compilation, you can also use `goofit_add_pdf(mypdf.cu)` to add a PDF. This will also require that you include any directory that you need with `include_directory`, as usual.

To add packages, use standard CMake tools. For example (CMake 3.5+), to add [Boost][FindBoost] 1.49+ filesystem and `TTreeReader` from ROOT:

```cmake
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.49 REQUIRED COMPONENTS filesystem)

goofit_add_executable(K3Pi K3Pi.cu)
target_link_libraries(MyNewExample Boost::filesystem ROOT::TreePlayer)
```

</p></details>


<details><summary>Adding a new project: (click to expand)</summary><p>

### External package (BETA)

GooFit now requires separable compilation, so it also now supports "external" packages, much like most other libraries. You can design your package with GooFit included as a subdirectory, and
it should just work. You'll also save time by not building examples, python bindings, and tests. The recommended procedure:

```bash
git add submodule <url to goofit> goofit
git submodule update --init --recursive
```

Then, you'll need a CMakeLists that looks something like this:

```bash
cmake_minimum_required(VERSION 3.9...3.16)

project(my_external_package LANGUAGES CXX)

add_subdirectory(goofit)
goofit_external_package()

goofit_add_executable(myapp myapp.cpp)
```

That's it! Just make a build directory and build. The `goofit_external_package()` command sets up optional CUDA, as well as links all reasonable files into your build directory. You can run `goofit_setup_std()`, `goofit_optional_cuda()` and `goofit_add_directory()` instead if you want.

### Classic method

If you'd like to make a separate GooFit project, you can do so. Simply checkout your project inside GooFit, with the name `work` or `goofit_`+something. CMake will automatically pick up those directories and build them, and GooFit's git will ignore them. Otherwise, they act just like the example directory. If you add a new directory, you will need to explicitly rerun CMake, as that cannot be picked up by the makefile. The automatic search can be turned off with the `GOOFIT_PROJECTS` option, or by using `GOOFIT_PROJECT_<name>` for a specific package.
GooFit packages should contain:

```cmake
goofit_add_package(MyPackageName)
```

After the package name, you can list `ROOT` to require that ROOT. The package will be disabled if ROOT is not found.

</p></details>




<details><summary style="font-size: 1.5em; margin-top: 24px; font-weight: 600; border-bottom: 1px solid #eaecef; line-height: 1.25">Using an IDE: (click to expand)</summary><p>

The following IDEs have been tested. Here `$SRC` refers to the source directory, and usually is `..` or `../GooFit`. You may want `-DCMAKE_BUILD_TYPE=Debug` and/or `-DGOOFIT_DEBUG=ON`.

| Name | Platform | Setup | Notes |
|------|----------|:------|:------|
| Xcode | macOS | `cmake $SRC -GXcode` | Only CPP version, works well though |
| Nsight-Eclipse | Linux | `cmake $SRC -G "Eclipse CDT4 - Unix Makefiles"` | Must be out-of-source, supports CUDA backend |
| QtCreator | All | Open from QtCreator dialog | Requires CMake extension (usually present). Might be able to use CMake 3.7+ Server |
| CLion | All | Open from CLion menu | Young but promising |

</p></details>


<details><summary>Converting from older GooFit Code: (click to expand)</summary><p>

The build system underwent a major upgrade in the move to CMake. The folders that were introduced to keep the includes structured require modifications of source code, converting lines like `#include "Variable.hh"` to `#include "GooFit/Variable.h"`. This modification can be done for you by running the provided script, `scripts/ModernizeGooFit.py` on your source files (requires Python and [Plumbum](https://github.com/tomerfiliba/plumbum)). You should remove your old Makefiles and use the new `CMakeFiles.txt` files provided in examples - this should require
writing two lines of code instead of the 50 or so previously needed. You should also add a GooFit Application to your code. (2 lines of CMake)

The new `GooFit::Application`, which is not required but provides GooFit options, like GPU selection and status, as well as MPI support and configurable command line options, is available by adding:

```cpp
#include "GooFit/Application.h"
using namespace GooFit;

// Place this at the beginning of main
Application app{"Optional description", argc, argv};

// Command line options can be added here.

GOOFIT_PARSE(app);
```

See [CLI11] for more details. The [pipipi0](./examples/pipipi0DPFit) example has an example of a complex set of options.

The other key differences in code are the addition of the `GooFit` namespace (`using namespace GooFit` allows fast conversion), and the removal of direct access to members of `Variable` (using getters/setters, or directly treat the variable like its value).

See Converting to [GooFit 2.0](./docs/CONVERTING20.md), [GooFit 2.1](./docs/CONVERTING21.md), and the [Changelog](./docs/CHANGELOG.md).

</p></details>


<details><summary>Improving performance with MPI: (click to expand)</summary><p>

Using the MPI version with an appropriate environment setup will allow for multiple GPU's to be used, and/or allow for multiple nodes.  To use this feature simply turn the flag on with CMake `-DGOOFIT_MPI=ON`.  This will divide the dataset by the number of processes involved.  For instance, if you have two nodes that will be involved in the calculation, the data will be split in half.  Currently, each node will load the entire buffer from disk, then load partitioned data it will work on.  It is highly recommended not to use more than one process per node for MPI+OpenMP versions.

A few notes about using the MPI version:

* You will need to use the `CountingVariable` for any event numbers used or referenced within the code, or anything that counts with the events.
* Please call `setDataSize` after `setData`.  If you do not, `setDataSize` doesn't have `m_iEventsPerTask`, which will need to be recalculated.

</p></details>


<details><summary>Configuring group size and grain size: (click to expand)</summary><p>

This advanced option is for GPU devices only. The script `scripts/find_optimal.py` will search a programmable group and grain space in order to find the optimal configuration for the particular PDFs.  This should be run after an example has been developed and tested.  Please look at `scripts/find_optimal.py` to see how to formulate a particular script.  Depending on the searchable space, this can take hours to days to compute.
The script will loop over the space and configure each parameter, then recompile and run the example a number of times.  A spreadsheet is calculated to help notice patterns, and the fastest version is printed to the user.

</p></details>


## Acknowledgement

GooFit's development is supported by the National Science Foundation under grant number [1414736]
and was developed under grant number [1005530].
Any opinions, findings, and conclusions or recommendations expressed in this material are those of the developers
and do not necessarily reflect the views of the National Science Foundation.
In addition, we thank the nVidia GPU Grant Program for donating hardware used in developing this framework.

GooFit is available under the BSD license, except for the Landau distribution & MINUIT code. You must remove these
to get a permissive version of GooFit.

[actions-badge]:     https://github.com/GooFit/GooFit/workflows/CI/badge.svg
[actions-link]:      https://github.com/GooFit/GooFit/actions
[DOI-badge]:         https://zenodo.org/badge/9017446.svg
[DOI-link]:          https://zenodo.org/badge/latestdoi/9017446
[API documentation]: https://GooFit.github.io/GooFit
[travis-badge]:      https://travis-ci.org/GooFit/GooFit.svg?branch=master
[travis-link]:       https://travis-ci.org/GooFit/GooFit
[gitter-badge]:      https://badges.gitter.im/GooFit/GooFit.svg
[gitter-link]:       https://gitter.im/GooFit/Lobby
[license-badge]:     https://img.shields.io/badge/License-BSD-blue.svg
[1005530]:           https://nsf.gov/awardsearch/showAward?AWD_ID=1005530
[1414736]:           https://nsf.gov/awardsearch/showAward?AWD_ID=1414736
[CUDA_SELECT_NVCC_ARCH_FLAGS]: https://cmake.org/cmake/help/v3.7/module/FindCUDA.html
[Plumbum]:           https://plumbum.readthedocs.io/en/latest/
[FindBoost]:         https://cmake.org/cmake/help/v3.7/module/FindBoost.html
[CLI11]:             https://github.com/CLIUtils/CLI11
[pybind11]:          http://pybind11.readthedocs.io/en/master
[ROOT]:              https://root.cern.ch
[Tutorials]:         https://goofit.gitlab.io/Goo2Torial
[pypi-status]:       https://img.shields.io/pypi/v/goofit.svg?logo=PyPI&logoColor=white
[pypi-link]:         https://pypi.python.org/pypi/goofit/
[cf-status]:         https://img.shields.io/conda/vn/conda-forge/goofit.svg?logo=Conda-Forge&logoColor=white
[cf-link]:           https://github.com/conda-forge/goofit-split-feedstock
[releases-badge]:    https://img.shields.io/github/release/GooFit/GooFit.svg
[releases-link]:     https://github.com/GooFit/GooFit/releases
[sk-badge]:          https://scikit-hep.org/assets/images/Scikit--HEP-Affiliated-blue.svg
