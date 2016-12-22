#include "FitControl.hh"
#include "PdfBase.hh" 

FitControl::FitControl (bool bin, std::string mn) 
  : binned(bin) 
  , metricName(mn)
  , owner(0)
  , errorsOnBins(false)
{}

FitControl::~FitControl () {} 

void FitControl::setOwner (PdfBase* dat) {
  assert(!owner); 
  owner = dat;
} 

UnbinnedNllFit::UnbinnedNllFit () 
  : FitControl(false, "ptr_to_NLL")
{}

BinnedNllFit::BinnedNllFit () 
  : FitControl(true, "ptr_to_BinAvg")
{}

BinnedErrorFit::BinnedErrorFit ()
  : FitControl(true, "ptr_to_BinWithError")
{
  errorsOnBins = true; 
}

BinnedChisqFit::BinnedChisqFit () 
  : FitControl(true, "ptr_to_Chisq")
{}


