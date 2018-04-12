## v2.1.3: Preparing for new indexing
#### April 12, 2018

This release is preparing for the major change in 2.2. The changes are mostly small, but help with several key pieces of needed functionality from the NIPS project. GooFit testing is significanly more stable.

User facing:

* Adding a few 3-body utilities: `get_amp_real`/`imag` and `getDecayInfo` ([#142])
* Generation no longer changes the number of events requested ([#142])
* Add support for fit fraction calculation (3 body only ATM) ([#143])
* kMatrix is now optional, and off by default for CUDA (slow build) ([#143])
* Can directly `pip install goofit` with Pip 10
* Packages can now be turned on and off individually (`goofit_add_package`)

Backend:

* Using Application no longer requires NVCC compile ([#139]), better C++ vs. CUDA separation
* Travis CI is much cleaner ([#140], [#141]), 
* Examples are part of tests now ([#144])
* Allow changing the bin number for normalization without restarting GooFit ([#143])
* Dropped symmDP where not supported ([#143])
* Dropped code coverage (test timeout issues when adding more tests) ([#144])
* Removed POSIX redefine warnings for Python, other warning fixes
* Packages do not trigger main repo git in binary dir
* Added feature summary

[#139]: https://github.com/GooFit/GooFit/pull/139
[#140]: https://github.com/GooFit/GooFit/pull/140
[#141]: https://github.com/GooFit/GooFit/pull/141
[#142]: https://github.com/GooFit/GooFit/pull/142
[#143]: https://github.com/GooFit/GooFit/pull/143
[#144]: https://github.com/GooFit/GooFit/pull/144


## v2.1.2: Live Python
#### March 12, 2018

GooFit has received even more Pythonization with the conversion of two more examples, ChiSquare and SimpleFit. It also has interactive Jupyter Notebook examples, where the plots update as you change the sliders.

* Grids can now be produced from DataSets with event numbers  ([#133])
* Far more of Minuit2 available from Python ([#138])
* MC generator for 1D functions (useful for testing in future) ([#138])
* Correlation matrices can now be accessed from Python ([#134])
* Some examples use numba for high-speed calculation of CPU code ([#134])
* Binned dataset access in Python improved ([#134])
* Travis timeouts reduced by separating codecov run
* Better error message if `Python.h` missing in Python build ([#136])
* OpenMP and CUDA updates to build system

[#133]: https://github.com/GooFit/GooFit/pull/133
[#134]: https://github.com/GooFit/GooFit/pull/134
[#136]: https://github.com/GooFit/GooFit/pull/136
[#138]: https://github.com/GooFit/GooFit/pull/138

## v2.1.1: Better Python
#### February 12, 2018

GooFit has received a few bugfixes that needed to be made available before the 2.2 merge. These include:

* More complete Python examples ([2437af9], [#131])
* DalitzPlotter helpers ([#131])
* Regression fix: FeatureDetector now respects current compilation options ([#131])

[#131]: https://github.com/GooFit/GooFit/pull/131
[2437af9]: https://github.com/GooFit/GooFit/commit/2437af9958e8534aa299b754a4e4361ae18d5301

## v2.1.0: Python bindings
#### December 7, 2017


GooFit now has working, fully supported Python bindings using PyBind11. All PDFs have been bound ([#92]). Python bindings now build by default if Python development files are found ([#93]). Pythonization additions, such as supporting the same shortcuts for Variables as C++, from/to numpy converters, and more, was added in ([#99], [#109]). Pip install is now supported using SciKit-Build, and source releases are being made on PyPI ([#107]). The build will use CUDA if found, and OpenMP otherwise.
Many examples converted ([#118], [#120], [#126])


Other Python additions:

* Live printout in Python Notebooks ([#114])
* Minuit2 wrapper started ([#115], [#128])
* `print_goofit_info` (and `print_splash`) added to Python  ([#126])
* `pyexamples/RunAll.sh` added to run all examples, run by Travis ([#126])


Major changes:

* `Observable`s are now their own class, and `CountingVariable` is now `EventNumber` ([#123])
* Variables are now passed by copy everywhere, handling smart pointers internally ([#124])
* DecayInfo renamed and split ([#124])
* FitControl is now explicitly a `std::shared_ptr` ([#126])
* Resonances and Lineshapes are now classes instead of using enums or ordering ([#119])
* OpenMP now supported on macOS Apple Clang on High Sierra with Homebrew, using `brew install cliutils/apple/libomp` ([#126])


Other changes include:

* New `ResonancePDF` types ([#114])
* Spline, KMatrix, FOCUS Lineshapes added, untested ([#119])
* TravisCI now uses Trusty ([#98]), performs style checks ([#117]), and runs all examples ([#114]).
* Minuit2 now can be missing from ROOT and GooFit will use its own copy ([#102], [#113]).
* Eigen is now included ([#104]), helper functions added ([#119])
* Large updates to CMake, including CUDA as a language ([#122]), CMake 3.9 FindCUDA backport, fully target-based build ([#119])
* Initial CMake IPO support ([#114])
* Better folder structure for PDF source files ([#114])
* `fpcomplex` shortcut ([#114])
* IDE support improvements ([#114]))
* Splash screen ([#114])
* Macros to help setup PDFs ([#119])
* CCache support
* `argc, argv` no longer required in `GOOFIT_PARSE` ([#126])
* The include order is now checked by the clang-format run ([#126])
* The ModernizeGooFit script supports 2.1 and more ([#126])
* `.cpp` is used everywhere, instead of `.cc` and `.cpp` mix
* CLI11 has been updated to 1.3 ([#130])

[#92]: https://github.com/GooFit/GooFit/pull/92
[#93]: https://github.com/GooFit/GooFit/pull/93
[#98]: https://github.com/GooFit/GooFit/pull/98
[#99]: https://github.com/GooFit/GooFit/pull/99
[#102]: https://github.com/GooFit/GooFit/pull/102
[#104]: https://github.com/GooFit/GooFit/pull/104
[#109]: https://github.com/GooFit/GooFit/pull/106
[#107]: https://github.com/GooFit/GooFit/pull/107
[#113]: https://github.com/GooFit/GooFit/pull/113
[#114]: https://github.com/GooFit/GooFit/pull/114
[#115]: https://github.com/GooFit/GooFit/pull/115
[#117]: https://github.com/GooFit/GooFit/pull/117
[#118]: https://github.com/GooFit/GooFit/pull/118
[#119]: https://github.com/GooFit/GooFit/pull/119
[#120]: https://github.com/GooFit/GooFit/pull/120
[#122]: https://github.com/GooFit/GooFit/pull/122
[#123]: https://github.com/GooFit/GooFit/pull/123
[#124]: https://github.com/GooFit/GooFit/pull/124
[#126]: https://github.com/GooFit/GooFit/pull/126
[#128]: https://github.com/GooFit/GooFit/pull/128
[#130]: https://github.com/GooFit/GooFit/pull/130


## v2.0.0: C++11 and CMake
#### June 9, 2017

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

