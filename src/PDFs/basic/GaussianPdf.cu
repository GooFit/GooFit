#include <goofit/Log.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/basic/GaussianPdf.h>

namespace GooFit {

__device__ auto device_Gaussian(fptype *evt, ParameterContainer &pc) -> fptype {
    int id       = pc.getObservable(0);
    fptype x     = RO_CACHE(evt[id]);
    fptype mean  = pc.getParameter(0);
    fptype sigma = pc.getParameter(1);
    pc.incrementIndex(1, 2, 0, 1, 1);

// mds fptype ret = exp(-0.5 * (x - mean) * (x - mean) / (sigma * sigma));
// to avoid potential problems with very, very small values,
// set all values for arg more than 10 sigma from center to be 2x10e-9
// [which is just below exp(-20.)]
    double arg = -0.5 * (x - mean) * (x - mean) / (sigma * sigma);
    fptype ret = (arg < -20.)?0.000000002:exp(arg);
// mds     std::cout << "in device_Gaussian: x, mean, sigma, ret = " <<
// mds                   x << " " << mean << "  " << sigma << "  " << ret << "\n";

    return ret;
}

__device__ device_function_ptr ptr_to_Gaussian = device_Gaussian;

__host__ GaussianPdf::GaussianPdf(std::string n, Observable _x, Variable mean, Variable sigma)
    : GooPdf("GaussianPdf", n, _x, mean, sigma) {

// mds     std::cout << "entered GaussianPdf::GaussianPdf \n";
// mds     PdfBase::status("    entered GaussianPdf::GaussianPdf");

    registerFunction("ptr_to_Gaussian", ptr_to_Gaussian);

    initialize();
// mds     std::cout << "  about to leave GaussianPdf::GaussianPdf \n";
// mds     PdfBase::status("      about to leave GaussianPdf::GaussianPdf");
}

__host__ auto GaussianPdf::integrate(fptype lo, fptype hi) const -> fptype {
    static const fptype rootPi = sqrt(atan2(0.0, -1.0));

    // Integral over all R.
    fptype sigma = host_parameters[parametersIdx + 2];
    sigma *= root2 * rootPi;
    return sigma;
}

} // namespace GooFit
