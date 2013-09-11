#------------------------------------------------------------------------------
CXX=nvcc
LD=g++  
OutPutOpt = -o

CXXFLAGS     = -O3 -arch=sm_20 
ifeq ($(CUDALOCATION), )
CUDALOCATION = /usr/local/cuda/
endif
CUDAHEADERS = $(CUDALOCATION)/include/

GOOFITDIR = ${PWD}/../../
SRCDIR = $(GOOFITDIR)/PDFs
WRKDIR = $(GOOFITDIR)/wrkdir/

INCLUDES += -I$(CUDAHEADERS) -I$(SRCDIR) -I$(GOOFITDIR) -I$(GOOFITDIR)/rootstuff 
LIBS += -L$(CUDALOCATION)/lib64 -lcudart -L$(GOOFITDIR)/rootstuff -lRootUtils 

# These are for user-level programs that want access to the ROOT plotting stuff, 
# not just the fitting stuff included in the GooFit-local ripped library. 
ROOT_INCLUDES = -I$(ROOTSYS)/include/
ROOT_LIBS     = -L$(ROOTSYS)/lib/ -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lMatrix -lPhysics -lMathCore -pthread -lThread -lMinuit -rdynamic -lFoam 

THRUSTO		= $(WRKDIR)/Variable.o $(WRKDIR)/FitManager.o $(WRKDIR)/GooPdfCUDA.o $(WRKDIR)/Faddeeva.o $(WRKDIR)/FitControl.o $(WRKDIR)/PdfBase.o $(WRKDIR)/DataSet.o $(WRKDIR)/BinnedDataSet.o $(WRKDIR)/UnbinnedDataSet.o $(WRKDIR)/FunctorWriter.o 
ROOTRIPDIR	= $(GOOFITDIR)/rootstuff
ROOTUTILLIB	= $(ROOTRIPDIR)/libRootUtils.so 

.SUFFIXES: 

all:	dalitz

%.o:	%.cu
	$(CXX) $(INCLUDES) $(ROOT_INCLUDES) $(DEFINEFLAGS) $(CXXFLAGS) -c -o $@ $<

dalitz:		DalitzFit.o $(THRUSTO) $(ROOTUTILLIB) 
		$(LD) $(LDFLAGS) $^ $(LIBS) $(ROOT_LIBS) $(OutPutOpt) $@
		@echo "$@ done"

clean:
		@rm -f *.o core dalitz dalitz.o 
