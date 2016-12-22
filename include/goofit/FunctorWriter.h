#ifndef FUNCTOR_WRITER_HH
#define FUNCTOR_WRITER_HH

class PdfBase; 

#include <thrust/host_vector.h>
#include "goofit/GlobalCudaDefines.h" // Need this for 'fptype' 

void writeToFile (PdfBase* pdf, const char* fname);
void readFromFile (PdfBase* pdf, const char* fname);

void writeListOfNumbers (thrust::host_vector<fptype>& target, const char* fname);
void readListOfNumbers (thrust::host_vector<fptype>& target, const char* fname); 

#endif
