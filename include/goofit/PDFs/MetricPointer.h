#include <goofit/FitControl.h>
#include <goofit/GlobalCudaDefines.h>

namespace GooFit {

/// Pass event, parameters, index into parameters.
typedef fptype (*device_metric_ptr)(fptype, fptype *, fptype);

void *getMetricPointer(EvalFunc val);

} // namespace GooFit
