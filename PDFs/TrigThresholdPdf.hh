#ifndef TRIGTHRESHOLD_PDF_HH
#define TRIGTHRESHOLD_PDF_HH

#include "EngineCore.hh" 

class TrigThresholdPdf : public EngineCore {
public:
  TrigThresholdPdf (std::string n, Variable* _x, Variable* thresh, Variable* trigConst, Variable* linConst, bool upper = true); 
  TrigThresholdPdf (std::string n, Variable* _x, Variable* _y, Variable* thresh, Variable* trigConst, Variable* linConst, Variable* massConstant, bool upper);

private:

};

#endif
