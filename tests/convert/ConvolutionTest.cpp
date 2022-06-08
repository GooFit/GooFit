#include <goofit/Catch.h>

#include "testhelpers.h"
#include <goofit/PDFs/basic/GaussianPdf.h>
#include <goofit/PDFs/combine/ConvolutionPdf.h>
#include <goofit/UnbinnedDataSet.h>

#include <goofit/Variable.h>

#include <random>

using namespace GooFit;

TEST_CASE("Convolution", "[convert][fit]") {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Observable xvar{"xvar", 0, 10};

    // Data set
    UnbinnedDataSet data(xvar);

    // Generate toy events.
    for(int i = 0; i < 10000; ++i) {
        double val = d(gen);
        if(val < 10 && val > 0) {
            xvar.setValue(val);
            data.addEvent();
        }
    }

    // Fit parameter
    Variable a1{"alpha1", 1, 0.1, 1, 10};
    Variable s1{"sigma1", 1, 0.1, 3};
    Variable a2{"alpha2", 1, 0.1, 1, 10};
    Variable s2{"sigma2", 1, 0.1, 3};

    GaussianPdf m{"gauss1", xvar, a1, s1};
    GaussianPdf r{"gauss2", xvar, a2, s2};

    // GooPdf object
    ConvolutionPdf pdf{"convolutionpdf", xvar, &m, &r};
    pdf.setData(&data);

    bool fitter = test_fitter(&pdf);

    CHECK(fitter);
    // CHECK(a1.getError() ==  .1);
    // CHECK(a1.getValue() == Approx(0.665178392).margin(a1.getError() * 3));
    // CHECK(a2.getError() ==  .1);
    // CHECK(a2.getValue() == Approx(0.665178392).margin(a2.getError() * 3));
}
