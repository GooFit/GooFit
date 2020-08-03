#include <goofit/PDFs/physics/StaircasePdf.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/Variable.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/PdfBase.h>

namespace GooFit {

__device__ fptype device_Staircase (fptype* evt, ParameterContainer &pc)
{

  int id = pc.getObservable(0);
  fptype x   = RO_CACHE(evt[id]);

  fptype ret = 0;

  for(unsigned int i = 0; i < 2; i++)
  {
    fptype param = pc.getParameter(i);
    if(x < param)
    {
      ret = i;
      break;
    }
  }

  pc.incrementIndex(1, 2, 0, 0, 1);

  return ret;
}

__device__ device_function_ptr ptr_to_Staircase = device_Staircase;
device_function_ptr hptr_to_Staircase = device_Staircase;

__host__ StaircasePdf::StaircasePdf(std::string n, Observable _x, std::vector<Variable> x0list)
  : GooPdf("StaircasePdf", n, _x) 
{

  for(Variable &x0: x0list)
    registerParameter(x0);

  registerFunction("ptr_to_Staircase", ptr_to_Staircase);

  initialize();

}


} //namespace GooFit