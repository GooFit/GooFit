#ifndef BW_PDF_HH
#define BW_PDF_HH

#include "goofit/PDFs/GooPdf.h" 

class BWPdf : public GooPdf {

public:
  BWPdf (std::string n, Variable* _x, Variable* m, Variable* s); 
private:

};

#endif
