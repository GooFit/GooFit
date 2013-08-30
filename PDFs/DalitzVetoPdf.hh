#ifndef DALITZVETO_PDF_HH
#define DALITZVETO_PDF_HH

#include "EngineCore.hh" 
#include "TddpPdf.hh"

struct VetoInfo {
  DaughterPair cyclic_index; 
  Variable* minimum;
  Variable* maximum; 
};

class DalitzVetoPdf : public EngineCore {
public:
  __host__ DalitzVetoPdf (std::string n,  Variable* _x, Variable* _y, Variable* motherM, Variable* d1m, Variable* d2m, Variable* d3m, vector<VetoInfo*> vetos);

private:

};

#endif
