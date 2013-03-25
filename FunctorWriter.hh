#ifndef FUNCTOR_WRITER_HH
#define FUNCTOR_WRITER_HH

class FunctorBase; 

#include <thrust/host_vector.h>
#include "GlobalCudaDefines.hh" // Need this for 'fptype' 

void writeToFile (FunctorBase* pdf, const char* fname);
void readFromFile (FunctorBase* pdf, const char* fname);

void writeListOfNumbers (thrust::host_vector<fptype>& target, const char* fname);
void readListOfNumbers (thrust::host_vector<fptype>& target, const char* fname); 

#endif
