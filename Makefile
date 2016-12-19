.PHONY:		goofit clean examples

GOOFITDIR = $(shell /bin/pwd)
include Makefile.common

FUNCTORLIST   = $(wildcard $(SRCDIR)/*.cu)
FUNCTORLIST   += $(wildcard $(SRCDIR)/*.hh)

# Handy print for makefile variables (debugging)
print-%  : ; @echo $* = $($*)

ROOTRIPOBJS	= $(ROOTRIPDIR)/TMinuit.o $(ROOTRIPDIR)/TRandom.o $(ROOTRIPDIR)/TRandom3.o 

goofit:		$(THRUSTO_B) $(ROOTUTILLIB) 
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

FitManager.o:		FitManager.cc FitManager.hh wrkdir/ThrustFitManagerCUDA.o Variable.o 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

wrkdir/AllPdfs.o:	$(FUNCTORLIST)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -I. -c PDFs/AllPdfs.cu -o $@ 
	@echo "$@ done"

examples:
	cd examples && $(MAKE)

clean:
	@rm -f *.o wrkdir/*
	cd rootstuff; $(MAKE) clean 

