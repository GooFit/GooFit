// the configured options
#define MCBooster_VERSION_MAJOR 1
#define MCBooster_VERSION_MINOR 1
#define MCBooster_VERSION_PATCH 0
/**
\mainpage Documentation

Table of Contents
=================

  * [MCBooster](#mcbooster)

    * [What is it?](#what-is-it)

    * [Main features](#main-features)

    * [The Latest Version](#the-latest-version)

    * [Documentation](#documentation)

    * [Installation and requirements ](#installation-and-requirements-)

    * [Examples](#examples)

    * [Licensing](#licensing)

    * [Contact the developers](#contact-the-developers)

    * [Author](#author)

    * [Acknowledgement](#acknowledgement)

What is it?
-----------

MCBooster is an header only library designed for the fast generation of
phase space events. The library makes use of Thrust and can deploy OpenMP
threads, CUDA and Xeon Phi cores. It is focused on performance and precision. 

The libray core algorithm follow as close as it is possible the implementation of the class [TGenPhaseSpace](https://root.cern.ch/doc/master/TGenPhaseSpace_8cxx.html)
from [ROOT framwork](https://root.cern.ch/),
which is based on the [GENBOD function (W515 from CERNLIB)](http://cernlib.web.cern.ch/cernlib/mc/genbod.html)
using the Raubold and Lynch method as described in 

[_F. James, Monte Carlo Phase Space, CERN 68-15 (1968)_](https://cds.cern.ch/record/275743/files/CERN-68-15.pdf).

Main features
-------------

Generates phase space Monte Carlo samples with up to nine particles in the final state, using very simple
and intuitive interface. Example:

```
	//generating 10M events of B0 -> J/psi K pi
	#include <mcbooster/GTypes.h>
	#include <mcbooster/Vector4R.h>
	#include <mcbooster/Generate.h>
	 ...
	//setting the mother particle
	Vector4R B0(5.2795, 0.0, 0.0, 0.0);
	 
	//setting the masses of the daughter particles
	vector<GReal_t> masses;
	masses.push_back(3.096916); // J/psi
	masses.push_back(0.493677); // K
	masses.push_back(0.13957018); // pi
	 
	//generator ctor for 10M events
	PhaseSpace phsp(B0.mass(), massesB0, 10000000);
	 
	//run the generator
	phsp.Generate(B0);
	 
	//Unweight the events flags the accepted and rejected events
	phsp.Unweight();
	 
	//export events to the host (in case it is necessary)
	Events *GenEvents = new Events(masses.size(), 10000000);
	phsp.Export(GenEvents);
	...
```

Other key features are:

1. Decays can be generated with mother particles are rest or with a definite four-vector.
this feature allows the generation of sequential decays.

2. Generates weighted and "unweighted" samples simultaneously. 

3. Allows the fast parallel evaluation of arbitrary functions taking as 
argument up to ten particles as input. 

4. Allows the fast parallel evaluation of arrays of variables simultaneously.

MCBooster also provides a bunch of custom types, containers and an increasing number of algorithms
to maximaze performance, avoid unecessary usage of memory and grant flexibility and protability between 
host and device calculations and deployment scenarios. 

Just changing .cu to .cpp in any source code writen only using the provided cosntructs is enough
to compile your application for OpenMP using GCC in a machine without a NVIDIA GPU installed.  

Many other possibilities and functionaties, bounded only by the creativity of the users. 

The Latest Version
------------------

The latest version can be found on the 
[project relases page](https://github.com/MultithreadCorner/MCBooster/releases).

Documentation
-------------

The complete and updated [Doxygen](http://www.doxygen.org/) source code documentation of this release is available in HTML format on the
[reference manual](http://multithreadcorner.github.io/MCBooster/) webpage.
Users can also browse the documentation by class, file or name using the following links:

1.[classes](http://multithreadcorner.github.io/MCBooster/classes.html)

2.[files](http://multithreadcorner.github.io/MCBooster/files.html)

3.[names](http://multithreadcorner.github.io/MCBooster/namespacemembers.html)

Installation and requirements 
-----------------------------

MCBooster is a header only library, so no build process is necessary to install it. 
Just place the `mcbooter` folder and its contents where your system can find it.
The library run on Linux systems and requires C++11 and the [Thrust library](https://thrust.github.io/). 
Some examples demonstrating the basic features of the library are included in the `src` folder. 
These code samples require [ROOT](https://root.cern.ch/) and [TCLAP](http://tclap.sourceforge.net/) library. 
CUDA based projects will require a local installation of [CUDA Tookit](https://developer.nvidia.com/cuda-toolkit) with version 6.5 or higher.  
Alternatively, projects the targeting [OpenMP](http://openmp.org/wp/) backend can be compiled with either nvcc or gcc. 
The CUDA runtime is not required to use OpemMP with gcc. 

Examples
--------

Some example code samples demonstrating the basic usage of the library are stored in the src directory, in the project source tree. 
These samples can be built using [CMAKE](https://cmake.org/) according the following instructions:

1. clone the git repository: `git clone https://github.com/MultithreadCorner/MCBooster.git`
2. go to MCBooster directory: `cd MCBooster`
3. create a build directory: `mkdir build` 
4. go to build directory: `cd build`
4. `cmake ../`
5. `make`

Users with root privilegies can do `make install` and get the targets installed into system-install-dir/bin 
(usually /usr/local. __Notice the project installation path is printed out in the setp 4__). Users without root privileges can point the installation path to a different location cmake `-DCMAKE_INSTALL_PREFIX=<user-path>/ ../`.

The examples are named according to the convention `MCBooster_Example_<BACKEND AND COMPILER>_<EXAMPLE NAME>`. To run an example do `./example-name`.
The examples are described below:

1. __B2KPiJpsi__ : Generates a sample of B0 -> Jpsi K pi, with J/psi -> mu+ mu- and calculates in parallel, for each event, the variables: 
  * M(K,pi), the Kpi invariant mass.
  * M(J/psi,pi), the J/psipi invariant mass.
  * cos theta(K), the helicity angle of the Kpi.
  * cos theta(mu), the helicity angle of the J/psi
  * phi, the angle between the decay planes 
  
The program print some events and timing information to sdtout and plotsthe distributions of the above variables plus the B0 -> J/psiK pi Dalitz plot.

2. __GenerateSample__ : Takes arguments from the command line, generates a sample and save it into a ROOT TTree. 

3. __PerformanceTest__: Meausure the time to generate and export samples in function of the number of events and number of particles.

Licensing
---------

MCBooster is released under the [GNU General Public License version 3](http://www.gnu.org/licenses/gpl-3.0.en.html). Please see the file called [COPYING](https://github.com/MultithreadCorner/MCBooster/blob/master/COPYING).

Contact the developers
----------------------
Here’s what you should do if you need help or would like to contribute:

1. If you need help or would like to ask a general question, subscribe and use https://groups.google.com/d/forum/mcbooster.
2. If you found a bug, use GitHub issues.
3. If you have an idea, use GitHub issues.
4. If you want to contribute, submit a pull request.

Author
--------

MCBooster was created and is mantained by [Antonio Augusto Alves Jr](https://github.com/AAAlvesJr).

Acknowledgement
---------------

MCBooster's development has been supported by the [National Science Foundation](http://nsf.gov/index.jsp) under grant number [PHY-1414736](http://nsf.gov/awardsearch/showAward?AWD_ID=1414736). Any opinions, findings, and conclusions or recommendations expressed in this material are those of the developers and do not necessarily reflect the views of the National Science Foundation.

 */
