#ifndef LANDAU_PDF_HH
#define LANDAU_PDF_HH

#include "EngineCore.hh" 

class LandauPdf : public EngineCore {
public:
  LandauPdf (std::string n, Variable* _x, Variable* mpv, Variable* sigma); 

private:

};

#endif
