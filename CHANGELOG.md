
## v2.0.0: CMake
#### In progress

The structure of GooFit and the build system have been completely revamped. The fake `nvcc` features have been removed, and the non-deterministic build failures from globbing have been eliminated (as has the globbing itself). The makefile system is much cleaner and nicer, but is also deprecated in favor of the new CMake builds. Many of the new features, like GooFit Packages, are not available unless using CMake. Several new examples and new PDFs have been added.

A new feature of the CMake build system is GooFit Packages, which are complete packages that can be added to GooFit and built, allowing analysis code to live in a separate location from GooFit, rather than the old method of simply forking GooFit and adding your analysis manually. A GooFit Package can be made into an example trivially. See [this package](https://github.com/maddocbf/goofit_KKPiPi) for an example.

#### Key features:

* Restructured files with script to aid in renaming includes
* CMake build system: See [Issue 22](https://github.com/GooFit/GooFit/issues/22) and [PR 23](https://github.com/GooFit/GooFit/pull/23).
  * Auto compute capability detection
  * Auto Cuda/OMP selection
  * Added CPP single threaded backend, support for MacOS
  * Optional separable compilation for PDFs, automatic for non-CUDA builds
  * (Almost) supports Intel compilers
  * Macros for `CMakeLists.txt` for adding a new package in 2-3 lines
  * Auto linking for build directory
* Fixes for building examples, nicer warnings with incorrect command line parameters.
* Examples have a script that run all of them with timing info
* Travis builds [PR 32](https://github.com/GooFit/GooFit/pull/32)
* Improved documentation, automatically builds on changes to master
* `GooFit::Application`, based on [CLI11](https://github.com/CLIUtils/CLI11). See [PR](https://github.com/GooFit/GooFit/pull/36) and [Issue](https://github.com/GooFit/GooFit/issues/33).
* Added MPI support in [PR 51](https://github.com/GooFit/GooFit/pull/36)
* Added PyGooFit: preliminary Python bindings using [PyBind11](http://pybind11.readthedocs.io/en/master/)
* Added (this) changelog

## Special tag: Final Makefile release
#### March 31, 2017

The Makefile system was partially maintained and adapted to the new file structure, but was deprecated after version 1.0, and received one special tag before it was removed. It was not possible to have an in-source CMake build while the makefile system was in place (and still not recommended).

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

