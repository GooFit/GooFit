CXX=g++ 
LD=g++  

CXXFLAGS     = -O2 -fPIC -pthread -g

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
			$(LD) -shared $(CXXFLAGS) -o $@ $^ 

clean:
		@rm -f *.o *.so 

