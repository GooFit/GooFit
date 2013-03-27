#------------------------------------------------------------------------------
CXX=nvcc
LD=g++  
OutPutOpt = -o

CXXFLAGS     = -O3 -m64 -arch=sm_20

CUDALOCATION = /usr/local/cuda/
CUDAHEADERS = $(CUDALOCATION)/include/
GOOFITDIR = $(PWD)/../

INCLUDES += -I$(CUDAHEADERS) -I$(GOOFITDIR) -I$(GOOFITDIR)/rootstuff -I$(GOOFITDIR)/FPOINTER/ 
LIBS += -L$(CUDALOCATION)/lib -lcudart -L$(GOOFITDIR)/rootstuff -lRootUtils

# These are for user-level programs that want access to the ROOT plotting stuff, 
# not just the fitting stuff included in the GooFit-local ripped library. 
ROOT_INCLUDES = -I$(ROOTSYS)/include/
ROOT_LIBS     = -L$(ROOTSYS)/lib/ -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lMatrix -lPhysics -lMathCore -pthread -lThread -lMinuit2 -lMinuit -lFoam
BUILD_DIR = $(GOOFITDIR)/build/

THRUSTO		= $(BUILD_DIR)/Variable.o $(BUILD_DIR)/PdfBuilder.o $(BUILD_DIR)/ThrustPdfFunctorCUDA.o $(BUILD_DIR)/Faddeeva.o $(BUILD_DIR)/FitControl.o $(BUILD_DIR)/FunctorBase.o $(BUILD_DIR)/DataSet.o $(BUILD_DIR)/BinnedDataSet.o $(BUILD_DIR)/UnbinnedDataSet.o $(BUILD_DIR)/FunctorWriter.o 
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
