GooFit is a massively-parallel framework, written in CUDA, for
doing maximum-likelihood fits. 

Requirements

* CUDA 4.2 or 5.0
* An nVidia GPU supporting compute capability at least 2.0

Installation

* Clone with git:

        git clone git://github.com/GooFit/GooFit.git

* check that the path setting in Makefile for CUDALOCATION is correct for your system
* then compile with gmake. 

Building the Examples Fits

* set the path to full root installation (bash syntax)

        export ROOTSYS=<path on your system>

* check that the path setting in examples/Makefile for CUDALOCATION is correct for your system
* change to the example directory and buid the program(s): 

        cd examples ; gmake all

Running the Example Fit

* Setup the bin and library paths for root

        export PATH=${PATH}:${ROOTSYS}/bin/:./
        # On Linux:
        export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ROOTSYS}/lib:<GooFit path>/rootstuff
        # On Mac:
        export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${ROOTSYS}/lib:<GooFit path>/rootstuff

* To run the simple fit example type:

        ./simpleFitExample

Acknowledgement

GooFit's development has been supported by the National Science Foundation under grant number NSF-1005530. 
Any opinions, findings, and conclusions or recommendations expressed in this material are those of the developers
and do not necessarily reflect the views of the National Science Foundation.
