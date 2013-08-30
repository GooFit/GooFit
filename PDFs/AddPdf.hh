#ifndef ADD_PDF_HH
#define ADD_PDF_HH

#include "EngineCore.hh" 

class AddPdf : public EngineCore {
public:

  AddPdf (std::string n, std::vector<Variable*> weights, std::vector<PdfBase*> comps); 
  AddPdf (std::string n, Variable* frac1, PdfBase* func1, PdfBase* func2); 
  __host__ virtual fptype normalise () const;
  __host__ virtual bool hasAnalyticIntegral () const {return false;}

protected:
  __host__ virtual double sumOfNll (int numVars) const;

private:
  bool extended; 
};

#endif
