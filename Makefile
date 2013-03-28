#------------------------------------------------------------------------------
CXX=nvcc
LD=g++  
OutPutOpt = -o

CXXFLAGS     = -O3 -arch=sm_20 
DEFINEFLAGS=-DDUMMY=dummy 

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

SRCDIR = $(PWD)/FPOINTER

INCLUDES += -I$(CUDAHEADERS) -I$(SRCDIR) -I$(PWD) -I$(PWD)/rootstuff 
LIBS += -L$(CUDALOCATION)/lib64 -lcudart -L$(PWD)/rootstuff -lRootUtils 

# These are for user-level programs that want access to the ROOT plotting stuff, 
# not just the fitting stuff included in the GooFit-local ripped library. 
ROOT_INCLUDES = -I$(ROOTSYS)/include/
ROOT_LIBS     = -L$(ROOTSYS)/lib/ -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lMatrix -lPhysics -lMathCore -pthread -lThread -lMinuit2 -lMinuit -rdynamic -lFoam 

FUNCTORLIST = $(SRCDIR)/ThrustPdfFunctor.cu 
FUNCTORLIST += $(wildcard $(SRCDIR)/*ThrustFunctor.cu)
FUNCTORLIST += $(wildcard $(SRCDIR)/*Aux.cu)
HEADERLIST = $(patsubst %.cu,%.hh,$(FUNCTORLIST))
WRKFUNCTORLIST = $(patsubst $(SRCDIR)/%.cu,wrkdir/%.cu,$(FUNCTORLIST))

THRUSTO		= wrkdir/Variable.o wrkdir/PdfBuilder.o wrkdir/ThrustPdfFunctorCUDA.o wrkdir/Faddeeva.o wrkdir/FitControl.o wrkdir/FunctorBase.o wrkdir/DataSet.o wrkdir/BinnedDataSet.o wrkdir/UnbinnedDataSet.o wrkdir/FunctorWriter.o 
ROOTRIPDIR	= $(PWD)/rootstuff
ROOTRIPOBJS	= $(ROOTRIPDIR)/TMinuit.o $(ROOTRIPDIR)/TRandom.o $(ROOTRIPDIR)/TRandom3.o 
ROOTUTILLIB	= $(ROOTRIPDIR)/libRootUtils.so 

.SUFFIXES: 
.PHONY:		goofit clean 

goofit:		$(THRUSTO)
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

PdfBuilder.o:		PdfBuilder.cc PdfBuilder.hh wrkdir/ThrustPdfFunctorCUDA.o Variable.o 
			$(CXX) $(DEFINEFLAGS) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

wrkdir/ThrustPdfFunctorCUDA.o:	wrkdir/CUDAglob.cu FunctorBase.cu 
				nvcc $(CXXFLAGS) $(INCLUDES) -I. $(DEFINEFLAGS) -c $< -o $@ 
				@echo "$@ done"

clean:
		@rm -f *.o wrkdir/*
		cd rootstuff; $(MAKE) clean 
