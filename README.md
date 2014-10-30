GooFit is a massively-parallel framework, written in CUDA, for
doing maximum-likelihood fits.   It is also possible to build
GooFit using OpenMP.

Requirements

* If using CUDA:
 * CUDA 4.2 or 5.x
 * An nVidia GPU supporting compute capability at least 2.0
* If building for OpenMP:
 * a compiler supporting OpenMP

Installation

* Clone with git:

        git clone git://github.com/GooFit/GooFit.git

* check that the path setting in Makefile for CUDALOCATION is correct for your system
 * If using OpenMP, CUDALOCATION should point at ./fakecuda
  * Set TARGET_OMP = 1 in the Makefile
  * Put ./fakcuda on your PATH, e.g. export PATH=$PATH:$HOME/GooFit/fakecuda
    (n.b. the use of g++ is hardwired in fakecuda/nvcc; this is bad)
  * Install a copy of CUDA's thrust in ./thrust (there's nothing to compile)
* then compile with gmake.

As an alternative to modifying the Makefile, you can pass the values as e.g.
  gmake CUDALOCATION=$HOME/Src/GooFit/fakecuda TARGET_OMP=1

Building the Examples Fits

* If you have Root on your system,
 * set the path to full root installation (bash syntax)
        export ROOTSYS=<path on your system>
* If you don't have Root, do NOT set ROOTSYS

* check that the path setting in examples/simpleFit/Makefile for CUDALOCATION is correct for your system
 * If using OpenMP, CUDALOCATION should point at ./fakecuda
  * Set TARGET_OMP = 1 in the Makefile
* change to the examples/simpleFit directory and buid the program(s): 
        cd examples/simpleFit ; gmake all
(you'll recall that you can set CUDALOCATION and TARGET_OMP on the command line instead of
modifying the Makefile)

Running the Example Fit

* Setup the bin and library paths for root

        export PATH=$PATH:$ROOTSYS/bin:.
        # On Linux:
        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib:<GooFit path>/rootstuff
        # On Mac:
        export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$ROOTSYS/lib:<GooFit path>/rootstuff
(If you don't have Root, don't modify the PATH and omit $ROOTSYS/lib:)
	
* To run the simple fit example type:
        ./simpleFitExample
(If you didn't define ROOTSYS the data and fits will be written to stdout)	

Acknowledgement

GooFit's development has been supported by the National Science Foundation under grant number NSF-1005530. 
Any opinions, findings, and conclusions or recommendations expressed in this material are those of the developers
and do not necessarily reflect the views of the National Science Foundation.
