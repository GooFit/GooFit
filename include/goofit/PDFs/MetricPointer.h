#include <goofit/FitControl.h>
#include <goofit/GlobalCudaDefines.h>

namespace GooFit {

/// Pass event, parameters, index into parameters.
typedef fptype (*device_metric_ptr)(fptype, fptype *, fptype);

auto getMetricPointer(EvalFunc val) -> void *;

} // namespace GooFit
