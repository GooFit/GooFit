#ifndef PROD_PDF_HH
#define PROD_PDF_HH

#include "goofit/PDFs/GooPdf.h" 

class ProdPdf : public GooPdf {
public:

  ProdPdf (std::string n, std::vector<PdfBase*> comps); 
  __host__ virtual fptype normalise () const;
  __host__ virtual bool hasAnalyticIntegral () const {return false;}

private:
  bool varOverlaps; // True if any components share an observable. 
};

#endif
