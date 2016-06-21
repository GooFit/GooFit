#------------------------------------------------------------------------------
CXX=nvcc
LD=g++  
OutPutOpt = -o

CXXFLAGS     = -w -arch=sm_20 -std=c++11 

ifeq ($(CUDALOCATION), )
CUDALOCATION = /usr/local/cuda/
endif
CUDAHEADERS = $(CUDALOCATION)/include/

GOOFITDIR = ${PWD}/../../
SRCDIR = $(GOOFITDIR)/PDFs
WRKDIR = $(GOOFITDIR)/wrkdir/

ifneq ($(TARGET_OMP),)
CXXFLAGS    += -fopenmp -DTHRUST_DEVICE_BACKEND=THRUST_DEVICE_BACKEND_OMP -lgomp -DMCBOOSTER_BACKEND=OMP -DTHRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_OMP
LDFLAGS     += -fopenmp -DTHRUST_DEVICE_BACKEND=THRUST_DEVICE_BACKEND_OMP -lgomp -DMCBOOSTER_BACKEND=OMP -DTHRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_OMP
else
LDFLAGS     += -DMCBOOSTER_BACKEND=CUDA
LIBS 		+= -L$(CUDALOCATION)/lib64 -lcudart
endif 

INCLUDES += -I$(CUDAHEADERS) -I$(SRCDIR) -I$(GOOFITDIR) -I$(GOOFITDIR)/rootstuff -I$(GOOFITDIR)/MCBooster
LIBS += -L$(GOOFITDIR)/rootstuff -lRootUtils 

# These are for user-level programs that want access to the ROOT plotting stuff, 
# not just the fitting stuff included in the GooFit-local ripped library. 
ROOT_INCLUDES = -I$(ROOTSYS)/include/
ROOT_LIBS     = -L$(ROOTSYS)/lib/ -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lMatrix -lPhysics -lMathCore -pthread -lThread -lMinuit -rdynamic -lFoam -lTreePlayer

THRUSTO		= $(WRKDIR)/Variable.o $(WRKDIR)/FitManager.o $(WRKDIR)/GooPdfCUDA.o $(WRKDIR)/Faddeeva.o $(WRKDIR)/FitControl.o $(WRKDIR)/PdfBase.o $(WRKDIR)/DataSet.o $(WRKDIR)/BinnedDataSet.o $(WRKDIR)/UnbinnedDataSet.o $(WRKDIR)/FunctorWriter.o 
ROOTRIPDIR	= $(GOOFITDIR)/rootstuff
ROOTUTILLIB	= $(ROOTRIPDIR)/libRootUtils.so 

.SUFFIXES: 

all:	tddp4

%.o:	%.cu
	$(CXX) $(INCLUDES) $(ROOT_INCLUDES) $(DEFINEFLAGS) $(CXXFLAGS) -c -o $@ $<

tddp4:		tddp4.o $(THRUSTO) $(ROOTUTILLIB) 
		$(LD) $(LDFLAGS) $^ $(LIBS) $(ROOT_LIBS) $(OutPutOpt) $@
		@echo "$@ done"

clean:
		@rm -f *.o tddp4 tddp4.o *.png *.root
		
