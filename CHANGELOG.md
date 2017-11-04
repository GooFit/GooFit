
## v2.0.0: C++11 and CMake

GooFit is now easy to build on a wide variety of Unix systems, and supports debuggers and IDEs. GooFit is faster, has unit tests, and working examples. More PDFs and examples have been added, as well as newly released example datasets that are downloaded automatically. GooFit now has built in support for MPI, and can use that to deploy to multiple graphics cards on the same machine. A new command line parser ([CLI11]) and drastically improved logging and errors have made code easier to write and debug. Usage of GooFit specific terminology is now reduced, using standard Thrust or CUDA terms when possible, lowering the barrier for new developers. A new Python script has been added to assist users converting from pre 2.0 code.

The file structure of GooFit and the build system have been completely revamped. The fake `nvcc` features have been removed, as have the `rootstuff` copies of ROOT classes. PDFs are now organized by type and compile and link separately. Multiple PDF caching support has improved.  The build system now uses CMake and manages external libraries.

A new feature of the CMake build system is GooFit Packages, which are complete packages that can be added to GooFit and built, allowing analysis code to live in a separate location from GooFit, rather than the old method of simply forking GooFit and adding your analysis manually. A GooFit Package can be made into an example trivially. See [this package](https://github.com/maddocbf/goofit_KKPiPi) for an example.

GooFit 2.0 will receive continuing support while development on GooFit 2.1 presses on.  

#### Key features:

* New 4-body Dalitz plot support and example
* New 4 body signal generation example and time-dependent Dalitz plot generation example powered by [MCBooster](https://github.com/MultithreadCorner/MCBooster)
* Restructured files with script to aid in renaming includes, available for assisting existing projects in converting to 2.0
* CMake build system: See [Issue 22](https://github.com/GooFit/GooFit/issues/22) and [PR 23](https://github.com/GooFit/GooFit/pull/23).
	* Auto compute capability detection based on detected GPUs
	* Auto CUDA/OMP selection 
	* Added CPP single threaded backend, support for macOS and IDEs/debuggers
	* Separable compilation for PDFs
	* Support for more compilers, such as Clang and Intel
	* Macros for `CMakeLists.txt` for adding a new package in 2-3 lines
	* Auto linking for build directory
	* Auto download of dependencies through git submodules and CMake
	* Data for examples automatically downloaded
* CUDA 7.0+ required, C++11 compliant compiler required. Large portions of the code have been moved to cleaner C++11 syntax, both by hand and with the `clang-tidy` tool ([PR 86](https://github.com/GooFit/GooFit/pull/88) and [PR 88](https://github.com/GooFit/GooFit/pull/88)).
* ROOT 6 recommended, but no longer required
* Fixes for building examples, nicer warnings with incorrect command line parameters
* Rootstuff, fakecuda, and other hacks removed ([PR 56](https://github.com/GooFit/GooFit/pull/56))
* Examples have a script that run all of them with timing info
* Travis CI builds ([PR 32](https://github.com/GooFit/GooFit/pull/32))
* Improved documentation, automatically builds on changes to master
* Added optional `GooFit::Application`, based on [CLI11](https://github.com/CLIUtils/CLI11), with standard GooFit options and logging output, fully user extendable for new options. See [PR 26](https://github.com/GooFit/GooFit/pull/36) and [Issue 33](https://github.com/GooFit/GooFit/issues/33).
* Better naming to match CUDA ([PR 61](https://github.com/GooFit/GooFit/pull/61))
* Added the GooFit namespace to all GooFit classes and variables.
* Better Variable based caching with multi-pdf support ([PR 65](https://github.com/GooFit/GooFit/pull/65) and [PR 68](https://github.com/GooFit/GooFit/pull/68))
* Logging and formatting support, cleanup of old commented code ([PR 66](/GooFit/GooFit/pull/66))
* Support for Minuit2 (default and available without ROOT) joins Minuit1, rebuilt fitters provide better output and automatic timing information
* CountingVariable added for MPI ready event numbers
* Added MPI support in [PR 51](https://github.com/GooFit/GooFit/pull/51), supporting multiple GPUs per node and multiple nodes
* Added preliminary Python bindings using [PyBind11](http://pybind11.readthedocs.io/en/master/)
* Started a new tutorial series, [GooFit 2torial](https://henryiii.gitbooks.io/goofit/content/), to replace [GooTorial](https://github.com/GooFit/GooTorial)
* Added a changelog and version information

[CLI11]: https://github.com/CLIUtils/CLI11

## Special tag: Final Makefile release
#### March 31, 2017

The Makefile system was partially maintained and adapted to the new file structure, but was deprecated after version 1.0, and received one special tag before it was removed. It is not possible to have an in-source CMake build.

## v1.0.0: Final Classic Makefile Release
#### December 17, 2016

This is the final release before the massive reworking of the build system. This was the "master" version of GooFit for some time. This release improved support for OMP and building with ROOT 6, but most of the work was done on a per example basis, so some examples still require ROOT 5.

The GitHub release includes the data needed for `pipipi0`.

This release includes speed improvements for compute arch 3.5+ boards, see [PR 21](https://github.com/GooFit/GooFit/pull/21). 

## v0.4: OMP Support
#### October 3, 2013

This release supports parallelising using OMP. To target OMP, do

```
gmake TARGET_OMP=yes
```

when compiling GooFit, and add the options

```
-fopenmp -DTHRUST_DEVICE_BACKEND=THRUST_DEVICE_BACKEND_OMP -lgomp
```

when compiling your own source code. (In the case of `nvcc`, you must have `-Xcompiler` in front of the `-fopenmp`). The example Makefiles show how to do this. 



## v0.3: The Great Rename
#### September 11, 2013

This changes the names of the core classes to better explain what they do. In particular `PdfBuilder` becomes `FitManager`; `ThrustPdfFunctor` becomes `GooPdf`; `FunctorBase` becomes `PdfBase`; the `FPOINTER` directory becomes PDFs; and `FooThrustFunctor` becomes `FooPdf`. 

NB! This breaks old code! To update your programs you can do this:

```bash
sed -i 's/ThrustPdfFunctor/GooPdf/g' `grep -l ThrustPdfFunctor *`
sed -i 's/FPOINTER/PDFs/g' `grep -l FPOINTER *`
sed -i 's/PdfFunctor/FitManager/g' `grep -l PdfFunctor *`
sed -i 's/FunctorBase/PdfBase/g' `grep -l FunctorBase *`
sed -i 's:\([a-zA-Z][a-zA-Z]*\)ThrustFunctor:\1Pdf:g' `grep -l "[a-zA-Z][a-zA-Z]*ThrustFunctor" *`
```

in that order. 

## v0.2: Docs and cleanup
#### July 17, 2013

## v0.1: First release
#### April 11, 2013

