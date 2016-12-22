#ifndef LANDAU_PDF_HH
#define LANDAU_PDF_HH

#include "goofit/PDFs/GooPdf.h" 

class LandauPdf : public GooPdf {
public:
  LandauPdf (std::string n, Variable* _x, Variable* mpv, Variable* sigma); 

private:

};

#endif
