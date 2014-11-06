ifeq ($(TARGET_MIC),)
CXX=g++
LD=g++
else
# Intel Xeon Phi/MIC requires using Intel C++ Compiler (ICC)
CXX=icpc
LD=icpc
CXXFLAGS=-mmic -x c++
LDFLAGS=-mmic
endif

CXXFLAGS     += -O2 -fPIC -pthread -g
LDFLAGS      += -fPIC -pthread -g

UNAME=$(shell uname)
ifeq ($(UNAME), Darwin)
CXXFLAGS+=-m64
endif

OBJECTS		= TMinuit.o TRandom.o TRandom3.o 

.SUFFIXES: 

all:	libRootUtils.so 

%.o:	%.cc %.hh 
	$(CXX) -I. $(CXXFLAGS) -c -o $@ $<

libRootUtils.so:	$(OBJECTS)
			$(LD) $(LDFLAGS) -shared -o $@ $^ 

clean:
		@rm -f *.o *.so 

