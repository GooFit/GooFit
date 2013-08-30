#ifndef EXPGAUS_PDF_HH
#define EXPGAUS_PDF_HH

#include "EngineCore.hh" 

class ExpGausPdf : public EngineCore {
public:
  ExpGausPdf (std::string n, Variable* _x, Variable* m, Variable* s, Variable* t); 

private:

};

#endif
