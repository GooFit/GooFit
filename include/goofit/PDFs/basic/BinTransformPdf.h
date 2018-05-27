#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

/**
Returns the global bin of its argument; in one
dimension:

\f[
    P(x;l,s) = \mathrm{floor}\left(\frac{x-l}{s}\right)
\f]

where \f$l\f$ is the lower limit and \f$s\f$ is the bin size. The utility of
this is perhaps not immediately obvious; one application is as an
intermediate step in a `MappedPdf`. For example, suppose I want to
model a \f$y\f$ distribution with a different set of parameters in five
slices of \f$x\f$; then I would use a `BinTransformPdf` to calculate
which slice each event is in.

The constructor takes `vector`s of the observables \f$\vec x\f$, lower
bounds \f$\vec l\f$, bin sizes \f$\vec b\f$, and number of bins \f$\vec n\f$.
The last is used for converting local (i.e. one-dimensional) bins
into global bins in the case of multiple dimensions.
**/
class BinTransformPdf : public GooPdf {
  public:
    // Transforms ND coordinates into a single bin number.
    BinTransformPdf(std::string n,
                    std::vector<Observable> obses,
                    std::vector<fptype> limits,
                    std::vector<fptype> binSizes,
                    std::vector<int> numBins);
};
} // namespace GooFit
