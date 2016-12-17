#------------------------------------------------------------------------------
ifeq ($(TARGET_MIC),)
CXX=nvcc
LD=g++  
else
# Intel Xeon Phi/MIC requires using Intel C++ Compiler (ICC)
CXX=icpc
LD=icpc
CXXFLAGS=-mmic -x c++
LDFLAGS=-mmic
endif

CXXFLAGS     += -O3 -std=c++11
CUDALIBDIR    = lib64

# ARMv7 32-bit (Nvidia Jetson TK1)
ifeq ($(shell uname -m),armv7l)
CUDALIBDIR    = lib
endif

ifeq ($(CUDALOCATION), )
CUDALOCATION = /usr/local/cuda/
endif
CUDAHEADERS = $(CUDALOCATION)/include/

GOOFITDIR = ${PWD}/../../
SRCDIR = $(GOOFITDIR)/PDFs
WRKDIR = $(GOOFITDIR)/wrkdir/

ifeq ($(TARGET_OMP),)
# nvcc (CUDA)
CXXFLAGS    += -DTARGET_SM35 -arch=sm_35
#CXXFLAGS    += -arch=sm_20
LIBS        += -L$(CUDALOCATION)/$(CUDALIBDIR) -lcudart
else
# OpenMP common flags
CXXFLAGS    += -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_BACKEND_OMP
LDFLAGS     += -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_BACKEND_OMP
ifeq ($(TARGET_MIC),)
# GCC/Clang
CXFFLAGS += -fopenmp
LDFLAGS += -fopenmp -lgomp
else
# Intel C++ Compiler (ICC)
CXXFLAGS += -openmp
LDFLAGS += -openmp
endif
endif 

INCLUDES += -I$(CUDAHEADERS) -I$(SRCDIR) -I$(GOOFITDIR) -I$(GOOFITDIR)/rootstuff 
LIBS += -L$(GOOFITDIR)/rootstuff -lRootUtils

# These are for user-level programs that want access to the ROOT plotting stuff, 
# not just the fitting stuff included in the GooFit-local ripped library. 
ROOT_INCLUDES = -I$(shell root-config --incdir)
ROOT_LIBS     = $(shell root-config --libs)

THRUSTO		= $(WRKDIR)/Variable.o $(WRKDIR)/FitManager.o $(WRKDIR)/GooPdfCUDA.o $(WRKDIR)/Faddeeva.o $(WRKDIR)/FitControl.o $(WRKDIR)/PdfBase.o $(WRKDIR)/DataSet.o $(WRKDIR)/BinnedDataSet.o $(WRKDIR)/UnbinnedDataSet.o $(WRKDIR)/FunctorWriter.o 
ROOTRIPDIR	= $(GOOFITDIR)/rootstuff
ROOTUTILLIB	= $(ROOTRIPDIR)/libRootUtils.so 

.SUFFIXES: 

all:	zach

%.o:	%.cu
	$(CXX) $(INCLUDES) $(ROOT_INCLUDES) $(DEFINEFLAGS) $(CXXFLAGS) -c -o $@ $<

zach:		ZachFit.o $(THRUSTO) $(ROOTUTILLIB) 
		$(LD) $(LDFLAGS) $^ $(LIBS) $(ROOT_LIBS) -o $@
		@echo "$@ done"

clean:
		@rm -f *.o core mixfit 
