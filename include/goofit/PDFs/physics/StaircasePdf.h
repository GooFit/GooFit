#ifndef STAIRCASE_THRUST_FUNCTOR_HH
#define STAIRCASE_THRUST_FUNCTOR_HH

#include <goofit/PDFs/GooPdf.h> 
#include <goofit/PdfBase.h>
#include <goofit/Variable.h>

#include <vector>
namespace GooFit {

class StaircasePdf : public GooPdf {
 public:
  StaircasePdf (std::string n, Observable _x, std::vector<Variable> x0list);
  //__host__ fptype integrate (fptype lo, fptype hi) const;
  //__host__ virtual bool hasAnalyticIntegral () const {return true;}
 private:
};
}
#endif
