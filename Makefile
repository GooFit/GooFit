#------------------------------------------------------------------------------
CXX=nvcc
LD=g++  
OutPutOpt = -o

CXXFLAGS     = -O3 -arch=sm_20 
DEFINEFLAGS=-DDUMMY=dummy 
CUDALIBDIR=lib64

UNAME=$(shell uname)
ifeq ($(UNAME), Darwin)
CUDALIBDIR=lib
CXXFLAGS+=-m64
endif

ifneq ($(CUDAPRINT),)
DEFINEFLAGS += -DCUDAPRINT=yes
endif 

ifneq ($(PRINTCALLS),)
DEFINEFLAGS += -DPRINTCALLS=yes
endif 

ifneq ($(PROFILE),)
DEFINEFLAGS += -DPROFILING=yes
endif 

CUDALOCATION = /usr/local/cuda/
CUDAHEADERS = $(CUDALOCATION)/include/

SRCDIR = $(PWD)/PDFs

INCLUDES += -I$(CUDAHEADERS) -I$(SRCDIR) -I$(PWD) -I$(PWD)/rootstuff 
LIBS += -L$(CUDALOCATION)/$(CUDALIBDIR) -lcudart -L$(PWD)/rootstuff -lRootUtils 

FUNCTORLIST = $(SRCDIR)/EngineCore.cu 
FUNCTORLIST += $(wildcard $(SRCDIR)/*Pdf.cu)
FUNCTORLIST += $(wildcard $(SRCDIR)/*Aux.cu)
HEADERLIST = $(patsubst %.cu,%.hh,$(FUNCTORLIST))
WRKFUNCTORLIST = $(patsubst $(SRCDIR)/%.cu,wrkdir/%.cu,$(FUNCTORLIST))
#NB, the above are used in the SRCDIR Makefile.

THRUSTO		= wrkdir/Variable.o wrkdir/FitManager.o wrkdir/EngineCoreCUDA.o wrkdir/Faddeeva.o wrkdir/FitControl.o wrkdir/PdfBase.o wrkdir/DataSet.o wrkdir/BinnedDataSet.o wrkdir/UnbinnedDataSet.o wrkdir/FunctorWriter.o 
ROOTRIPDIR	= $(PWD)/rootstuff
ROOTRIPOBJS	= $(ROOTRIPDIR)/TMinuit.o $(ROOTRIPDIR)/TRandom.o $(ROOTRIPDIR)/TRandom3.o 
ROOTUTILLIB	= $(ROOTRIPDIR)/libRootUtils.so 

.SUFFIXES: 
.PHONY:		goofit clean 

goofit:		$(THRUSTO) $(ROOTUTILLIB) 
		@echo "Built GooFit objects" 

# One rule for GooFit objects.
wrkdir/%.o:	%.cc %.hh 
		@mkdir -p wrkdir 
		$(CXX) $(INCLUDES) $(CXXFLAGS) $(DEFINEFLAGS) -c -o $@ $<

# A different rule for user-level objects. Notice ROOT_INCLUDES. 
%.o:	%.cu
	$(CXX) $(INCLUDES) $(ROOT_INCLUDES) $(DEFINEFLAGS) $(CXXFLAGS) -c -o $@ $<

# Still a third rule for the ROOT objects - these have their own Makefile. 
$(ROOTRIPDIR)/%.o:	$(ROOTRIPDIR)/%.cc 
			rm -f $@ 
			@echo "Postponing $@ for separate Makefile" 

$(ROOTUTILLIB):	$(ROOTRIPOBJS)
		@cd rootstuff; $(MAKE) 

include $(SRCDIR)/Makefile 

FitManager.o:		FitManager.cc FitManager.hh wrkdir/ThrustFitManagerCUDA.o Variable.o 
			$(CXX) $(DEFINEFLAGS) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

wrkdir/EngineCoreCUDA.o:	wrkdir/CUDAglob.cu PdfBase.cu 
				nvcc $(CXXFLAGS) $(INCLUDES) -I. $(DEFINEFLAGS) -c $< -o $@ 
				@echo "$@ done"

clean:
		@rm -f *.o wrkdir/*
		cd rootstuff; $(MAKE) clean 
