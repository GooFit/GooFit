


## v2.0.0: CMake
#### In progress

See [Issue 22](https://github.com/GooFit/GooFit/issues/22) and [PR 23](https://github.com/GooFit/GooFit/pull/23).


## v1.0.0: Final Classic Makefile Release
#### December 17, 2016

This is the final release before the massive reworking of the build system. This was the "master" version of GooFit for some time. This release improved support for OMP and building with ROOT 6, but most of the work was done on a per example basis, so some examples still require ROOT 5.

The github release includes the data needed for `pipipi0`.

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

