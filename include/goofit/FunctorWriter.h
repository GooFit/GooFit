#pragma once

#include <thrust/host_vector.h>

#include "goofit/GlobalCudaDefines.h" // Need this for 'fptype'

namespace GooFit {

class PdfBase;

void writeToFile(PdfBase *pdf, const char *fname);
void readFromFile(PdfBase *pdf, const char *fname);

void writeListOfNumbers(thrust::host_vector<fptype> &target, const char *fname);
void readListOfNumbers(thrust::host_vector<fptype> &target, const char *fname);

} // namespace GooFit
