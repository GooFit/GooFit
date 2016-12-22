#ifndef MAPPED_PDF_HH
#define MAPPED_PDF_HH

#include "goofit/PDFs/GooPdf.h" 

class MappedPdf : public GooPdf {
public:
  MappedPdf (std::string n, GooPdf* m, vector<GooPdf*>& t); 
  // Map function m must be custom written to correspond to order of function list t. 
  __host__ fptype normalise () const;
private:

};

#endif
