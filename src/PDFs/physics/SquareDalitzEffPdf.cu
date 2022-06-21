#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/SquareDalitzEffPdf.h>
#include <vector>

namespace GooFit {

__device__ auto thetaprime(fptype m12, fptype m13, fptype mD, fptype mKS0, fptype mh1, fptype mh2) -> fptype {
    // Helper function to calculate theta'
    fptype m23 = mD * mD + mKS0 * mKS0 + mh1 * mh1 + mh2 * mh2 - m12 - m13;
    if(m23 < 0)
        return 0;

    fptype num     = m23 * (m12 - m13);
    fptype lambda1 = ((m23 - mh1 * mh1 + mh2 * mh2) * (m23 - mh1 * mh1 + mh2 * mh2) - 4 * m23 * mh2 * mh2);
    fptype lambda2 = ((mD * mD - mKS0 * mKS0 - m23) * (mD * mD - mKS0 * mKS0 - m23) - 4 * m23 * mKS0 * mKS0);
    fptype denum   = sqrt(lambda1) * sqrt(lambda2);
    if(lambda1 < 0. || lambda2 < 0.) {
        GOOFIT_TRACE("WARNING NAN.");
        return -99;
    }
    if(denum != 0.) {
        return num / denum;
    } else {
        GOOFIT_TRACE("WARNING deunum is zero.");
        return -99;
    }
}

__device__ auto device_SquareDalitzEff(fptype *evt, ParameterContainer &pc) -> fptype {
    // Define observables
    int idx = pc.getObservable(0);
    int idy = pc.getObservable(1);

    // don't use RO_CACHE as this is used as efficiency for Amp3Body
    fptype x = evt[idx];
    fptype y = evt[idy];

    // Define coefficients
    fptype c0 = pc.getParameter(0);
    fptype c1 = pc.getParameter(1);
    fptype c2 = pc.getParameter(2);
    fptype c3 = pc.getParameter(3);
    fptype c4 = pc.getParameter(4);
    fptype c5 = pc.getParameter(5);
    fptype c6 = pc.getParameter(6);
    fptype c7 = pc.getParameter(7);

    fptype mD   = 1.86483;
    fptype mKS0 = 0.497611;
    fptype mh1  = 0.13957;
    fptype mh2  = 0.13957;

    pc.incrementIndex(1, 8, 0, 2, 1);

    // Check phase space
    if(!inDalitz(x, y, mD, mKS0, mh1, mh2))
        return 0;

    // Call helper functions
    fptype thetap = thetaprime(x, y, mD, mKS0, mh1, mh2);
    if(thetap > 1. || thetap < -1.)
        return 0;

    fptype m23 = mD * mD + mKS0 * mKS0 + mh1 * mh1 + mh2 * mh2 - x - y;
    if(m23 < 0 || m23 > 2)
        return 0;

    fptype ret = c0 * m23 * m23 + c1 * m23 + c2 * m23 * thetap * thetap + c3 * thetap * thetap + c4 * thetap + c5
                 + c6 * m23 * m23 * m23 * m23 + c7 * m23 * m23 * m23;

    return ret;
}

__device__ device_function_ptr ptr_to_SquareDalitzEff = device_SquareDalitzEff;

__host__ SquareDalitzEffPdf::SquareDalitzEffPdf(std::string n,
                                                Observable m12,
                                                Observable m13,
                                                Variable c0,
                                                Variable c1,
                                                Variable c2,
                                                Variable c3,
                                                Variable c4,
                                                Variable c5,
                                                Variable c6,
                                                Variable c7)

    : GooPdf("SquareDalitzEffPdf", n, m12, m13, c0, c1, c2, c3, c4, c5, c6, c7) {
    registerFunction("ptr_to_SquareDalitzEff", ptr_to_SquareDalitzEff);

    initialize();
}

// __host__ fptype SquareDalitzEffPdf::normalize() { return 1; }

} // namespace GooFit
