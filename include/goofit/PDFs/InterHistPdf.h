#pragma once

#include "goofit/BinnedDataSet.h"
#include "goofit/PDFs/GooPdf.h"
#include <thrust/device_vector.h>

namespace GooFit {

class InterHistPdf : public GooPdf {
  public:
    InterHistPdf(std::string n, BinnedDataSet *x, std::vector<Variable *> params, std::vector<Variable *> obses);
    //__host__ virtual fptype normalize () const;

  private:
    thrust::device_vector<fptype> *dev_base_histogram;
    fptype totalEvents;
    fptype *host_constants;
    int numVars;
};

} // namespace GooFit
