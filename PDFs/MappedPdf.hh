#ifndef MAPPED_PDF_HH
#define MAPPED_PDF_HH

#include "EngineCore.hh" 

class MappedPdf : public EngineCore {
public:
  MappedPdf (std::string n, EngineCore* m, vector<EngineCore*>& t); 
  // Map function m must be custom written to correspond to order of function list t. 
  __host__ fptype normalise () const;
private:

};

#endif
