.PHONY:		goofit clean examples

GOOFITDIR = $(shell /bin/pwd)
include Makefile.common

FUNCTORLIST   = $(wildcard src/PDFs/*.cu)
FUNCTORLIST   += $(wildcard include/goofit/PDFs/*.h)

warning:
	@echo "Depreciated: please use the CMake build. A classic makefile build is still available with make all or make goofit."
	@echo "The standard procedure:"
	@echo "  mkdir build"
	@echo "  cd build"
	@echo "  cmake .."
	@echo "  make"
	@echo ""
	@echo "Remember: you can make multiple build directories to have different options, like build-cuda and build-omp." 
	@echo "You can use 'make omp cuda' to setup or rebuild the build-cuda and build-omp directories."

omp:
	@mkdir -p build-omp
	cd build-omp && cmake .. -DGOOFIT_DEVICE=OMP && $(MAKE) --no-print-directory

cuda:
	@mkdir -p build-cuda
	cd build-cuda && cmake .. -DGOOFIT_DEVICE=CUDA && $(MAKE) --no-print-directory

all: goofit examples

goofit: $(THRUSTO_B) wrkdir/libRootUtils.so
	@echo "Built GooFit objects" 

# Handy print for makefile variables (debugging)
print-%  : ; @echo $* = $($*)

# One rule for GooFit objects.
wrkdir/%.o:	src/goofit/%.cc include/goofit/%.h 
	@mkdir -p wrkdir
	$(CXX) $(INCLUDES) $(CXXFLAGS) -c -o $@ $<


# A different rule for user-level objects. Notice ROOT_INCLUDES. 
%.o:	%.cu
	$(CXX) $(INCLUDES) $(ROOT_INCLUDES)  $(CXXFLAGS) -c -o $@ $<

# Still a third rule for the ROOT objects - these have their own Makefile. 
$(ROOTRIPDIR)/%.o:	$(ROOTRIPDIR)/%.cc 
	rm -f $@ 
	@echo "Postponing $@ for separate Makefile" 

wrkdir/libRootUtils.so:
	$(MAKE) -f Makefile.rootstuff 

FitManager.o:		FitManager.cc FitManager.h wrkdir/ThrustFitManagerCUDA.o Variable.o 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

wrkdir/AllPdfs.o:	$(FUNCTORLIST)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -I$(SRCDIR)/PDFs  -c src/PDFs/AllPdfs.cu -o $@ 
	@echo "$@ done"

examples:
	cd examples && $(MAKE)

clean:
	@rm -f  wrkdir/*
	cd examples && $(MAKE) clean

