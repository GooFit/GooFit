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
CXXFLAGS    += -fopenmp -DTHRUST_DEVICE_BACKEND=THRUST_DEVICE_BACKEND_OMP -lgomp
LDFLAGS     += -fopenmp -DTHRUST_DEVICE_BACKEND=THRUST_DEVICE_BACKEND_OMP -lgomp
else
LIBS 		+= -L$(CUDALOCATION)/lib64 -lcudart
endif 

INCLUDES += -I$(CUDAHEADERS) -I$(SRCDIR) -I$(GOOFITDIR) -I$(GOOFITDIR)/rootstuff -I$(GOOFITDIR)/MCBooster
LIBS += -L$(CUDALOCATION)/lib64 -lcudart -L$(GOOFITDIR)/rootstuff -lRootUtils 

# These are for user-level programs that want access to the ROOT plotting stuff, 
# not just the fitting stuff included in the GooFit-local ripped library. 
ROOT_INCLUDES = -I$(ROOTSYS)/include/
ROOT_LIBS     = -L$(ROOTSYS)/lib/ -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lMatrix -lPhysics -lMathCore -pthread -lThread -lMinuit -rdynamic -lFoam -lTreePlayer

THRUSTO		= $(WRKDIR)/Variable.o $(WRKDIR)/FitManager.o $(WRKDIR)/GooPdfCUDA.o $(WRKDIR)/Faddeeva.o $(WRKDIR)/FitControl.o $(WRKDIR)/PdfBase.o $(WRKDIR)/DataSet.o $(WRKDIR)/BinnedDataSet.o $(WRKDIR)/UnbinnedDataSet.o $(WRKDIR)/FunctorWriter.o 
ROOTRIPDIR	= $(GOOFITDIR)/rootstuff
ROOTUTILLIB	= $(ROOTRIPDIR)/libRootUtils.so 

.SUFFIXES: 

all:	plots

%.o:	%.cu
	$(CXX) $(INCLUDES) $(ROOT_INCLUDES) $(DEFINEFLAGS) $(CXXFLAGS) -c -o $@ $<

%.o:	%.cpp
	g++ $(ROOT_INCLUDES) -w -g -std=c++11 -c -o $@ $<

mcevents:		mcevents.o 
		$(LD) $^ $(ROOT_LIBS) $(OutPutOpt) $@

		@echo "$@ done"
DP4:		DP4.o $(THRUSTO) $(ROOTUTILLIB) 
		$(LD) $(LDFLAGS) $^ $(LIBS) $(ROOT_LIBS) $(OutPutOpt) $@
		@echo "$@ done"

phspMC.root: mcevents
			./mcevents

DP4.txt: DP4 phspMC.root
		./DP4 > DP4.txt

plots: DP4.txt
	chmod 775 plot.py
	./plot.py
clean:
		@rm -f *.o DP4 DP4.o DP4.txt *.png mcevents phspMC.root
