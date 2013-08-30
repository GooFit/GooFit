#ifndef VOIGTIAN_PDF_HH
#define VOIGTIAN_PDF_HH

#include "EngineCore.hh" 

class VoigtianPdf : public EngineCore {
public:
  VoigtianPdf (std::string n, Variable* _x, Variable* m, Variable* s, Variable* w); 

private:

};

#endif
