CXX=g++ 
LD=g++  

CXXFLAGS     = -O2 -m64 -fPIC -pthread -g

OBJECTS		= TMinuit.o TRandom.o TRandom3.o 

.SUFFIXES: 

all:	libRootUtils.so 

%.o:	%.cc %.hh 
	$(CXX) -I. $(CXXFLAGS) -c -o $@ $<

libRootUtils.so:	$(OBJECTS)
			$(LD) -shared $(CXXFLAGS) -o $@ $^ 
#			$(LD) -shared -Wl,-soname,$@ $(CXXFLAGS) -o $@ $^ 

clean:
		@rm -f *.o *.so 

