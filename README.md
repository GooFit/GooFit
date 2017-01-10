GooFit is a massively-parallel framework, written in CUDA, for
doing maximum-likelihood fits with a comfortable syntax.
It is also possible to build
GooFit using OpenMP.

## Requirements

* If using CUDA:
 * CUDA 7.5 and 8.0 tested, older versions should work for now
 * An nVidia GPU supporting compute capability at least 2.0 (3.5 recommended)
* If using OpenMP:
 * A compiler supporting OpenMP (GCC4.8+ and Intel 17 tested)

## Getting the files

* Clone with git:
```
git clone git://github.com/goofit/goofit.git
cd goofit
```

If you want to live on the bleeding edge, or are interested in contributing to GooFit, checkout the dev branch:
```
git checkout dev
```


## Building 

The build system now uses CMake. The procedure is standard for CMake builds:
```
mkdir build
cd build
cmake ..
make
```

If you want to change compiler, set `CC` and `CXX` to appropriate defaults *before* you run cmake either inline or in your environment. If you want to set the host and device backends, you can set those options. The defaults are:
```
cmake .. -DGOOFIT_DEVICE=CUDA -DGOOFIT_HOST=OMP
```

Valid options are `CUDA` (device only), `OMP`, `CPP`, and `TBB` (unavailable currently).

A few standard cmake tricks:

* Use `VERBOSE=1` to see the commands used to build the files
* Use `-L` to list the CMake options
* Use `ccmake` if available to see a curses (terminal) gui, or `cmake-gui` for a completely graphical interface
* Use `-G` and the name of a generator to use something other than `make`, like `Xcode` or `ninja`
* Open the `CMakeLists.txt` with QtCreator to generate for that IDE
* Set the release type with `-DCMAKE_BUILD_TYPE=Release`, `RelWithDebInfo`, `Debug`, etc.
* Set up multiple build directories, like `build-omp`
* CMake remembers your `-D` option selections in your build directory so you don't have to specify them again
* CMake reruns when needed when you `make`



## Classic Makefile system (depreciated)

* You should set `CUDALOCATION` for your system

* You should have run source thisroot.sh to setup ROOT paths and other environment variables

* Set `TARGET_OMP=1` if you want to use OMP
  * Checkout a copy of CUDA's thrust next to the goofit repository (there's nothing to compile)
  * If you already have thrust in a different location, set `THRUSTLOCATION`

```
make all TARGET_OMP=1
```

Following the fairly standard makefile convention, all programs are built in-place instead of in a build directory.

## Running the Examples

* To run all the examples, with timing information, use:
```
./examples/int_testing.py
```

(This requires the Plumbum library, install with `pip install plumbum`, `pip install --user plumbum`, or `conda -c conda-forge plumbum`)

If you want to run an individual example, those are in subdirectories in examples (built products are in your build directory, the source is in `goofit/examples`).


## Adding a new example:

The examples are designed to be easy to add to. Make a new directory, then add a new CMakeLists.txt in your directory with one or more of the following two lines:

```
goofit_add_executible(MyNewExample MyNewExample.cu)
goofit_add_link(SomeNeededDataFile.txt)
```

The first line adds your `.cu` file with goofit code as an executible, and the second one sets up a symbolic link to the data file (or any file) in the build directory to the source directory. To get the example to build when you build goofit, add the name of your directory to `examples/CMakeLists.txt`.

> ## Note:
> If you want to extend the Makefile system instead, copy a Makefile from a different directory, changing the relevent project name (only one program per directory supported), and make a new target in `examples/Makefile`. 

## Acknowledgement

GooFit's development has been supported by the National Science Foundation under grant number NSF-1005530. 
Any opinions, findings, and conclusions or recommendations expressed in this material are those of the developers
and do not necessarily reflect the views of the National Science Foundation.
