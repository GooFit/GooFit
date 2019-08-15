#include <mcbooster/Evaluate.h>
#include <mcbooster/GContainers.h>
#include <mcbooster/GFunctional.h>
#include <mcbooster/GTypes.h>
#include <mcbooster/functors/FlagAcceptReject.h>

#include <thrust/extrema.h>

#include <goofit/PDFs/physics/Amp4BodyBase.h>
#include <goofit/PDFs/physics/Amp4BodyGlobals.h>

namespace GooFit {
/// Start with mcbooster::BoolVector_d flags(numEvents);
void Amp4BodyBase::fillMCFlags(mcbooster::BoolVector_d &flags,
                               const mcbooster::RealVector_d &weights,
                               unsigned int numEvents) {
    thrust::counting_iterator<mcbooster::GLong_t> first(0);
    thrust::counting_iterator<mcbooster::GLong_t> last = first + numEvents;

    auto max = thrust::max_element(weights.begin(), weights.end());
    thrust::transform(first, last, weights.begin(), flags.begin(), mcbooster::FlagAcceptReject((fptype)*max));
}

} // namespace GooFit
