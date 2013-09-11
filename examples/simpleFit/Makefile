#------------------------------------------------------------------------------
CXX=nvcc
LD=g++  
OutPutOpt = -o

CXXFLAGS     = -O3 -arch=sm_20 -g 
CUDALIBDIR=lib64

UNAME=$(shell uname)
ifeq ($(UNAME), Darwin)
CUDALIBDIR=lib
CXXFLAGS+=-m64
else
# Adds library symbols to debug info 
# - not strictly required, nice to have.
# Turn it off if it gives you trouble. 
LDFLAGS += -rdynamic 
endif

ifeq ($(CUDALOCATION), )
CUDALOCATION = /usr/local/cuda/
endif
CUDAHEADERS = $(CUDALOCATION)/include/
GOOFITDIR = $(PWD)/../../

INCLUDES += -I$(CUDAHEADERS) -I$(GOOFITDIR) -I$(GOOFITDIR)/rootstuff -I$(GOOFITDIR)/PDFs/ 
LIBS += -L$(CUDALOCATION)/$(CUDALIBDIR) -lcudart -L$(GOOFITDIR)/rootstuff -lRootUtils 

# These are for user-level programs that want access to the ROOT plotting stuff, 
# not just the fitting stuff included in the GooFit-local ripped library. 
ROOT_INCLUDES = -I$(ROOTSYS)/include/
ROOT_LIBS     = -L$(ROOTSYS)/lib/ -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lMatrix -lPhysics -lMathCore -pthread -lThread -lMinuit2 -lMinuit -lFoam 
WRKDIR = $(GOOFITDIR)/wrkdir/

THRUSTO		= $(WRKDIR)/Variable.o $(WRKDIR)/FitManager.o $(WRKDIR)/GooPdfCUDA.o $(WRKDIR)/Faddeeva.o $(WRKDIR)/FitControl.o $(WRKDIR)/PdfBase.o $(WRKDIR)/DataSet.o $(WRKDIR)/BinnedDataSet.o $(WRKDIR)/UnbinnedDataSet.o $(WRKDIR)/FunctorWriter.o 
ROOTRIPDIR	= $(GOOFITDIR)/rootstuff
ROOTUTILLIB	= $(ROOTRIPDIR)/libRootUtils.so 
PROGRAMS	= simpleFitExample

.SUFFIXES: 

all:	$(PROGRAMS)

%.o:	%.cu
	$(CXX) $(INCLUDES) $(ROOT_INCLUDES) $(DEFINEFLAGS) $(CXXFLAGS) -c -o $@ $<

simpleFitExample:	simpleFitExample.o 
			$(LD) $(LDFLAGS) $^ $(THRUSTO) $(LIBS) $(ROOT_LIBS) $(OutPutOpt) $@
			@echo "$@ done"

clean:
		@rm -f *.o core $(PROGRAMS) 
