## v2.3.0: Lots of polish
## May 23, 2022

Several years of fixes, updates, and polish. GooFit's core is now
BSD licensed (still uses LGPL Landau & Minuit code).

* Minimum supported CUDA version: 8.0
* Minimum required CMake version: 3.9 (matches ROOT) [#240][]
* Lots of dependent version updates (pybind11, CLI11, etc)
* All files formatted with pre-commit [#241][]
* 3-body preliminary kMatrix [#209][]
* 3 gauss resolution updates (now has one extra param) [#218][]
* LASS parameters now settable [#218][]
* Changed to Level to GlobalLevel for ROOT 6.24 [#266][]
* GeneateSig available on Amp3Body [#254][]
* Added rho omega lineshape [#278][]
* New `PdfBase::status` for debugging [#298][]
* Nicer particle array access from Python [#307][]

[#218]: https://github.com/GooFit/GooFit/pull/218
[#209]: https://github.com/GooFit/GooFit/pull/209
[#240]: https://github.com/GooFit/GooFit/pull/240
[#241]: https://github.com/GooFit/GooFit/pull/241
[#254]: https://github.com/GooFit/GooFit/pull/254
[#266]: https://github.com/GooFit/GooFit/pull/266
[#278]: https://github.com/GooFit/GooFit/pull/278
[#278]: https://github.com/GooFit/GooFit/pull/278
[#298]: https://github.com/GooFit/GooFit/pull/298
[#307]: https://github.com/GooFit/GooFit/pull/307


## v2.2.3: Python builds
#### January 31, 2020

This release sees a huge increase in platforms built on CI, and lots
of fixes to assist the new conda-forge package.

* GooFit version is now available as a Python tuple
* Check CUDA 8.0, 9.2, and 10.0 builds on GHA [#224][]
* Allow `goofit -m python` when no GPU present [#229][]
* Support `GOOFIT_ARCH=OFF`, fix `FORCE_LOCAL_THRUST` for CMake 3.12+ [#224][]
* Bump vendored library versions [#223][], [#225][]
* Square Dalitz added [#207][]
* Added `setDataSize`, `getName`, `getBinSize` to Python bindings
* CrystalBall fixes [#200][]
* 2D MC bound [#196][]
* Make standalone source package available on GitHub [#232][]

[#196]: https://github.com/GooFit/GooFit/pull/196
[#200]: https://github.com/GooFit/GooFit/pull/200
[#207]: https://github.com/GooFit/GooFit/pull/207
[#223]: https://github.com/GooFit/GooFit/pull/223
[#224]: https://github.com/GooFit/GooFit/pull/224
[#225]: https://github.com/GooFit/GooFit/pull/225
[#229]: https://github.com/GooFit/GooFit/pull/229
[#232]: https://github.com/GooFit/GooFit/pull/232


## v2.2.2: Python access
#### February 1, 2019

This adds a few minor fixes and tools for easier debugging.

* Added M12 and M23 access to DalitzPlotter [#190]
* Fix non-CUDA and older CMake builds [#191]
* Allow environment variables to set Python build [#192]
* Better variable formatting, fixed `print_goofit_info` [#193], [#194], more information on CUDA and ROOT [#195]
* `app.set_floating_exceptions()` will cause floating point errors to throw for debugging
* `app.get_filename` now supports GooFit as a package
* Fix for Python segfaulting when calling `evaluatePdf` multiple times
* GSpline fixes [#194]
* Support ROOT 6.04 [#195]
* Uses CLI11 1.7.1 now [#195]

[#190]: https://github.com/GooFit/GooFit/pull/190
[#191]: https://github.com/GooFit/GooFit/pull/191
[#192]: https://github.com/GooFit/GooFit/pull/192
[#193]: https://github.com/GooFit/GooFit/pull/193
[#194]: https://github.com/GooFit/GooFit/pull/194
[#195]: https://github.com/GooFit/GooFit/pull/195



## v2.2.1: Bugfixes
#### October 23, 2018

This release contains some performance improvements [#176], and
fixes some of the bugs present in the rewrite of 2.2.0:

* Fixed bug with smooth histogram [#178]
* Restored spline, derivatives are now properly calculated too [#181], limits corrected [#188]
* MCBooster seed now available [#177]

It also fixes several building problems found by users:

* Fixed builds with CMake 3.5 (ROOT Docker uses this by default)
* Better Pip 10+ support, limiting SciKit-Build version (0.7+ buggy on macOS) [#188]
* Docker is used to format if clang-format is not installed and docker is
* Forking now supported again [#179], git without https access easier [#183]
* Better CUDA 9+ support [#180], followed by full CUDA 9+ support [#189]
* LCG on SLC 6 supported [#185], [#184]
* OpenMP supported again on latest macOS (not on Anaconda)
* Python binding bugfixes [#188]

[#176]: https://github.com/GooFit/GooFit/pull/176
[#177]: https://github.com/GooFit/GooFit/pull/177
[#178]: https://github.com/GooFit/GooFit/pull/178
[#179]: https://github.com/GooFit/GooFit/pull/179
[#180]: https://github.com/GooFit/GooFit/pull/180
[#181]: https://github.com/GooFit/GooFit/pull/181
[#183]: https://github.com/GooFit/GooFit/pull/183
[#184]: https://github.com/GooFit/GooFit/pull/184
[#185]: https://github.com/GooFit/GooFit/pull/185
[#188]: https://github.com/GooFit/GooFit/pull/188
[#189]: https://github.com/GooFit/GooFit/pull/189

## v2.2.0: New indexing
#### July 31, 2018

The internals of GooFit have been updated for GPU performance and much simpler PDF authoring.  The new improvements will make PDF generation and debugging much simpler. [#125] Lots of Python improvements are part of this change, as well, including the removal of several compile time limits. Much better Python help display. The Physics PDFs have undergone a rename, with a coherent naming scheme (old names still work for now).

* Indexing improves performance and readability
* A number of new tests demonstrating simple usage of each PDF. [#148]
* Travis CI uses the additional tests to ensure changes work well.
* `registerFunction` simplifies PDF construction [#151]
* The Catch2 test framework makes tests easier to write [#151]
* 1D Monte Carlo is easy to generate, and test with [#151]
* Lots of improvements to error checking, `GOOFIT_MAXFUNC` added [#151]
* FitControl is always enum based now [#151]
* Packages can now request changes on the current configure run (bugfix) [#151]
* Improved help, also with Jupyter notebook customization [#154], [#160]
* Parameters can now be registered in the GooPdf Constructor [#154]
* Separable compilation is now required, faster compile times [#159]
* Nicer Makefile, fixed normalize (internal) spelling [#164]
* ParameterContainer is checked when `GOOFIT_TRACE` is on and using CPU [#164]
* Normalization improvements, PDFs store cache directly [#164]
* Removed defunct code [#165]
* Better file separation, better functional separation [#166]
    * `sumOfNll` removed
    * Added `get_event_size`, `get_bin_grid_size`, `reduce_with_metric`, `reduce_with_bins`, `evaluate_with_metric`
* No longer limit the number of parameters and functions at compile time [#163]
    * Deleting a PDF or an Application will clear the GPU parameter/function memory
* Internal rename, PDF structure, and breaking up of files [#167]
    * All the Dalitz PDFs are now in AmpNBody (Old names and includes still provided for compatibility)
    * Resonances and Lineshapes are now one per file
    * All helper classes are split out and available in `PDFs/physics/detail`
    * Inheritance structure added
* Better testing for MPI [#169]
* Partial support for CUDA 9.x [#172]
* External package support (BETA) [#174]

[#125]: https://github.com/GooFit/GooFit/pull/125
[#148]: https://github.com/GooFit/GooFit/pull/148
[#151]: https://github.com/GooFit/GooFit/pull/151
[#154]: https://github.com/GooFit/GooFit/pull/154
[#159]: https://github.com/GooFit/GooFit/pull/159
[#160]: https://github.com/GooFit/GooFit/pull/160
[#163]: https://github.com/GooFit/GooFit/pull/163
[#164]: https://github.com/GooFit/GooFit/pull/164
[#165]: https://github.com/GooFit/GooFit/pull/165
[#166]: https://github.com/GooFit/GooFit/pull/166
[#167]: https://github.com/GooFit/GooFit/pull/167
[#169]: https://github.com/GooFit/GooFit/pull/169
[#172]: https://github.com/GooFit/GooFit/pull/172
[#174]: https://github.com/GooFit/GooFit/pull/174

## v2.1.3: Preparing for new indexing
#### April 21, 2018

This release is the last feature release for the 2.1 series. The changes are mostly focused on stability and support, but there are also several key pieces of needed functionality added from the NIPS project. GooFit testing is significantly more stable, C++/CUDA are better separated, one more C++ example is now available in Python, and 6.10 and 6.12 versions of ROOT are finally supported.

User facing:

* Adding a few 3-body utilities: `get_amp_real`/`imag` and `getDecayInfo` ([#142])
* Generation no longer changes the number of events requested ([#142])
* Add support for fit fraction calculation (3 body only ATM) ([#143])
* kMatrix is now optional, and off by default for CUDA (slow build) ([#143])
* Can directly `pip install goofit` with Pip 10
* Packages can now be turned on and off individually (`goofit_add_package`)
* Moved files users might use to `goofit/utilities` from `goofit/detail` ([#145])
* Adds symmetrized support for RBW (no runtime impact if not requested!) ([#145])
* Python example available for Zachfit ([#149])


Backend:

* Using Application no longer requires NVCC compile ([#139]), better C++ vs. CUDA separation
* Travis CI is much cleaner ([#140], [#141]),
* Examples are part of tests now ([#144])
* Allow changing the bin number for normalization without restarting GooFit ([#143])
* Dropped symmDP where not supported ([#143])
* Dropped code coverage (test timeout issues when adding more tests) ([#144])
* Removed POSIX redefine warnings for Python (except CUDA), other warning fixes
* Build most examples and Python bindings with C++
* ROOT versions newer than 6.08 are now supported
* Packages do not trigger main repo git in binary dir
* Added feature summary
* Cleaned up GlobalCudaDefines, dropped GooError in favor of CUDA naming ([#145])
* Added some repairs for CUDA 7 ([#146])
* Added some minor changes to C++ examples ([#149])
* `setIntegrationconstant` Python bindings added ([#149])


[#139]: https://github.com/GooFit/GooFit/pull/139
[#140]: https://github.com/GooFit/GooFit/pull/140
[#141]: https://github.com/GooFit/GooFit/pull/141
[#142]: https://github.com/GooFit/GooFit/pull/142
[#143]: https://github.com/GooFit/GooFit/pull/143
[#144]: https://github.com/GooFit/GooFit/pull/144
[#145]: https://github.com/GooFit/GooFit/pull/145
[#146]: https://github.com/GooFit/GooFit/pull/146
[#149]: https://github.com/GooFit/GooFit/pull/149

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


GooFit now has working, fully supported Python bindings using pybind11. All PDFs have been bound ([#92]). Python bindings now build by default if Python development files are found ([#93]). Pythonization additions, such as supporting the same shortcuts for Variables as C++, from/to NumPy converters, and more, was added in ([#99], [#109]). Pip install is now supported using SciKit-Build, and source releases are being made on PyPI ([#107]). The build will use CUDA if found, and OpenMP otherwise.
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
* Added preliminary Python bindings using [pybind11](http://pybind11.readthedocs.io/en/master/)
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
