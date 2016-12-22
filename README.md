GooFit is a massively-parallel framework, written in CUDA, for
doing maximum-likelihood fits.   It is also possible to build
GooFit using OpenMP.

## Requirements

* If using CUDA:
 * CUDA 7.5 and 8.0 tested, older versions should work for now
 * An nVidia GPU supporting compute capability at least 2.0
 * A compiler supporting OpenMP

## Installation

* Clone with git:
```
git clone git://github.com/goofit/goofit.git
```

* You should set `CUDALOCATION` for your system

* Set `TARGET_OMP=1` if you want to use OMP
  * Checkout a copy of CUDA's thrust next to the goofit repository (in `../thrust`) (there's nothing to compile)
  * If you already have thrust, set `THRUSTLOCATION`

```
make TARGET_OMP=1
```

## Building the Examples Fits

* Run `source thisroot.sh` from the bin directory of your root install to set up your environment for root
* Change to the examples directory and build the program(s): 
```
cd examples && make
```
or
```
cd examples && make TARGET_OMP=1
```

## Running the Example Fit

* To run the simple fit example type:
```
./simpleFitExample
```

Acknowledgement

GooFit's development has been supported by the National Science Foundation under grant number NSF-1005530. 
Any opinions, findings, and conclusions or recommendations expressed in this material are those of the developers
and do not necessarily reflect the views of the National Science Foundation.
