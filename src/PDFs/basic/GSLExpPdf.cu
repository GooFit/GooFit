#include <goofit/PDFs/basic/GSLExpPdf.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <gsl/gsl_sf_exp.h>

using namespace std;

namespace GooFit {

__device__ fptype device_GSLExpPdf(fptype *evt, ParameterContainer &pc) {
    int id       = pc.getObservable(0);
    fptype alpha = pc.getParameter(0);
    fptype x     = RO_CACHE(evt[id]);

    pc.incrementIndex(1, 1, 0, 1, 1);
    return gsl_sf_exp_mult(x, alpha);
}

__device__ device_function_ptr ptr_to_GSLExpPdf = device_GSLExpPdf;

__host__ GSLExpPdf::GSLExpPdf(std::string n, Observable _x, Variable alpha)
    : GooPdf("GSLExpPdf", n, _x, alpha) {
    registerFunction("ptr_to_GSLExp", ptr_to_GSLExpPdf);

    initialize();
}

} // namespace GooFit
