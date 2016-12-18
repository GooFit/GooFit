.PHONY:		goofit clean 

GOOFITDIR = $(shell /bin/pwd)
include Makefile.common

# GooPdf must be first in CUDAglob, as it defines global variables.
FUNCTORLIST    = $(SRCDIR)/GooPdf.cu 
FUNCTORLIST    += $(filter-out $(SRCDIR)/GooPdf.cu, $(wildcard $(SRCDIR)/*Pdf.cu))
FUNCTORLIST	   += $(SRCDIR)/DalitzPlotHelpers.cu $(SRCDIR)/SpinFactors.cu
FUNCTORLIST   += $(wildcard $(SRCDIR)/*Aux.cu)
HEADERLIST     = $(patsubst %.cu,%.hh,$(SRCFILES))
WRKFUNCTORLIST = $(patsubst $(SRCDIR)/%.cu,wrkdir/%.cu,$(FUNCTORLIST))
#NB, the above are used in the SRCDIR Makefile.

ROOTRIPOBJS	= $(ROOTRIPDIR)/TMinuit.o $(ROOTRIPDIR)/TRandom.o $(ROOTRIPDIR)/TRandom3.o 

goofit:		$(THRUSTO) $(ROOTUTILLIB) 
		@echo "Built GooFit objects" 

# One rule for GooFit objects.
wrkdir/%.o:	%.cc %.hh 
		@mkdir -p wrkdir 
		$(CXX) $(INCLUDES) $(CXXFLAGS) -c -o $@ $<

# A different rule for user-level objects. Notice ROOT_INCLUDES. 
%.o:	%.cu
	$(CXX) $(INCLUDES) $(ROOT_INCLUDES)  $(CXXFLAGS) -c -o $@ $<

# Still a third rule for the ROOT objects - these have their own Makefile. 
$(ROOTRIPDIR)/%.o:	$(ROOTRIPDIR)/%.cc 
			rm -f $@ 
			@echo "Postponing $@ for separate Makefile" 

$(ROOTUTILLIB):	$(ROOTRIPOBJS)
		@cd rootstuff && $(MAKE) 

include $(SRCDIR)/Makefile 

FitManager.o:		FitManager.cc FitManager.hh wrkdir/ThrustFitManagerCUDA.o Variable.o 
			$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

wrkdir/GooPdfCUDA.o:	wrkdir/CUDAglob.cu PdfBase.cu 
			$(CXX) $(CXXFLAGS) $(INCLUDES) -I. -c $< -o $@ 
			@echo "$@ done"

clean:
		@rm -f *.o wrkdir/*
		cd rootstuff; $(MAKE) clean 
