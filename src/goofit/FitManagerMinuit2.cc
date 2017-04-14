#include "goofit/PdfBase.h"
#include "goofit/FitManagerMinuit2.h"
#include "goofit/PDFs/GooPdf.h"
#include "goofit/Variable.h"
#include "goofit/FCN.h"

#include <Minuit2/MnUserParameters.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/FunctionMinimum.h>
#include "Minuit2/MnPrint.h"

namespace GooFit {
    
FitManagerMinuit2::FitManagerMinuit2(PdfBase* dat) {
    fcn_.reset(new FCN(dat));
}

Minuit2::FunctionMinimum FitManagerMinuit2::fit() {
    auto val = Minuit2::MnPrint::Level();
    Minuit2::MnPrint::SetLevel(3);
    
    // Setting global call number to 0
    host_callnumber = 0;
    Minuit2::MnUserParameters* upar = fcn_->get_params();
    
    Minuit2::MnMigrad migrad{*fcn_, *upar};
    Minuit2::FunctionMinimum min = migrad();
    std::cout << min << std::endl;
    
    Minuit2::MnPrint::SetLevel(val);
    return min;
}
    
}
  

