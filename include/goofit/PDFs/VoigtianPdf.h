#ifndef VOIGTIAN_PDF_HH
#define VOIGTIAN_PDF_HH

#include "goofit/PDFs/GooPdf.h" 

class VoigtianPdf : public GooPdf {
public:
  VoigtianPdf (std::string n, Variable* _x, Variable* m, Variable* s, Variable* w); 

private:

};

#endif
