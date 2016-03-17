#------------------------------------------------------------------------------
CXX=nvcc
LD=g++  
OutPutOpt = -o

CXXFLAGS     = -w -arch=sm_20 -std=c++11 -O3 
#-O3
ifeq ($(CUDALOCATION), )
CUDALOCATION = /usr/local/cuda/
endif
CUDAHEADERS = $(CUDALOCATION)/include/

GOOFITDIR = ${PWD}/../../../
SRCDIR = $(GOOFITDIR)/PDFs
WRKDIR = $(GOOFITDIR)/wrkdir/

ifneq ($(TARGET_OMP),)
CXXFLAGS    += -Xcompiler -fopenmp -DTHRUST_DEVICE_BACKEND=THRUST_DEVICE_BACKEND_OMP -lgomp
LDFLAGS     += -fopenmp -DTHRUST_DEVICE_BACKEND=THRUST_DEVICE_BACKEND_OMP -lgomp
endif 

INCLUDES += -I$(CUDAHEADERS) -I$(SRCDIR) -I$(GOOFITDIR) -I$(GOOFITDIR)/rootstuff 
LIBS += -L$(CUDALOCATION)/lib64 -lcudart -L$(GOOFITDIR)/rootstuff -lRootUtils 

# These are for user-level programs that want access to the ROOT plotting stuff, 
# not just the fitting stuff included in the GooFit-local ripped library. 
ROOT_INCLUDES = -I$(ROOTSYS)/include/
ROOT_LIBS     = -L$(ROOTSYS)/lib/ -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lMatrix -lPhysics -lMathCore -pthread -lThread -lMinuit -rdynamic -lFoam 

THRUSTO		= $(WRKDIR)/Variable.o $(WRKDIR)/FitManager.o $(WRKDIR)/GooPdfCUDA.o $(WRKDIR)/Faddeeva.o $(WRKDIR)/FitControl.o $(WRKDIR)/PdfBase.o $(WRKDIR)/DataSet.o $(WRKDIR)/BinnedDataSet.o $(WRKDIR)/UnbinnedDataSet.o $(WRKDIR)/FunctorWriter.o 
ROOTRIPDIR	= $(GOOFITDIR)/rootstuff
ROOTUTILLIB	= $(ROOTRIPDIR)/libRootUtils.so 

.SUFFIXES: 

all:	dalitz

%.o:	%.cu
	$(CXX) $(INCLUDES) $(ROOT_INCLUDES) $(DEFINEFLAGS) $(CXXFLAGS) -c -o $@ $<

dalitz:		DalitzFit2.o $(THRUSTO) $(ROOTUTILLIB) 
		$(LD) $(LDFLAGS) $^ $(LIBS) $(ROOT_LIBS) $(OutPutOpt) $@
		@echo "$@ done"

dalitz_plot:		DalitzFit_plot.o $(THRUSTO) $(ROOTUTILLIB) 
		$(LD) $(LDFLAGS) $^ $(LIBS) $(ROOT_LIBS) $(OutPutOpt) $@
		@echo "$@ done"

dalitz2:		DalitzFit2.o $(THRUSTO) $(ROOTUTILLIB) 
		$(LD) $(LDFLAGS) $^ $(LIBS) $(ROOT_LIBS) $(OutPutOpt) $@
		@echo "$@ done"

dalitz_D3pi_sig:		DalitzFit_D3pi_sig.o $(THRUSTO) $(ROOTUTILLIB) 
		$(LD) $(LDFLAGS) $^ $(LIBS) $(ROOT_LIBS) $(OutPutOpt) $@
		@echo "$@ done"

dalitz_D3K_sig:		DalitzFit_D3K_sig.o $(THRUSTO) $(ROOTUTILLIB) 
		$(LD) $(LDFLAGS) $^ $(LIBS) $(ROOT_LIBS) $(OutPutOpt) $@
		@echo "$@ done"

dalitz2_gen:		DalitzFit2_gen.o $(THRUSTO) $(ROOTUTILLIB) 
		$(LD) $(LDFLAGS) $^ $(LIBS) $(ROOT_LIBS) $(OutPutOpt) $@
		@echo "$@ done"

dalitz3:		DalitzFit3.o $(THRUSTO) $(ROOTUTILLIB) 
		$(LD) $(LDFLAGS) $^ $(LIBS) $(ROOT_LIBS) $(OutPutOpt) $@
		@echo "$@ done"

dalitz_symm:		DalitzFit_symm.o $(THRUSTO) $(ROOTUTILLIB) 
		$(LD) $(LDFLAGS) $^ $(LIBS) $(ROOT_LIBS) $(OutPutOpt) $@
		@echo "$@ done"

dp_plot:	dp_plot.o $(THRUSTO) $(ROOTUTILLIB) 
		$(LD) $(LDFLAGS) $^ $(LIBS) $(ROOT_LIBS) $(OutPutOpt) $@
		@echo "$@ done"

clean:
		@rm -f *.o core dalitz dalitz.o 
