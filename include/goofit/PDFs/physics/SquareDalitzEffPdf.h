
#pragma once

#include <vector>

#include <goofit/GlobalCudaDefines.h> // Need this for 'fptype'
#include <goofit/PDFs/GooPdf.h>
#include <goofit/Variable.h>
#include <goofit/detail/Complex.h>

namespace GooFit {

class SquareDalitzEffPdf : public GooPdf {
  public:
    // Very specific efficiency parametrisation for semileptonically-tagged D0->KSPiPi decays as determined from data
    // Uses variables of square Dalitz plot - m' and theta'
    SquareDalitzEffPdf(std::string n,
                       Observable m12,
                       Observable m13,
                       Variable c0,
                       Variable c1,
                       Variable c2,
                       Variable c3,
                       Variable c4,
                       Variable c5,
                       Variable c6,
                       Variable c7);

  private:
};

} // namespace GooFit
