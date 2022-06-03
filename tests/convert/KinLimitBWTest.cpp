#include <goofit/Catch.h>

#include "testhelpers.h"
#include <goofit/PDFs/basic/KinLimitBWPdf.h>
#include <goofit/UnbinnedDataSet.h>

#include <goofit/Variable.h>

#include <random>

using namespace GooFit;

TEST_CASE("Gaussian", "[convert][fit]") {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Observable xvar{"xvar", -10, 10};

    // Data set
    UnbinnedDataSet data(xvar);

    // Generate toy events.
    for(int i = 0; i < 1000; ++i) {
        double val = d(gen);
        if(val < 10) {
            xvar.setValue(val);
            data.addEvent();
        }
    }

    // Fit parameter
    Variable alpha{"alpha", 1, 0.1, 1, 10};
    Variable sigma{"sigma", 1, 1, 3};

    // GooPdf object
    KinLimitBWPdf kinlimitbwpdf{"kinlimitbwpdf", xvar, alpha, sigma};
    kinlimitbwpdf.setData(&data);

    bool fitter = test_fitter(&kinlimitbwpdf);

    CHECK(fitter);
    CHECK(alpha.getError() < .1);
    CHECK(alpha.getValue() == Approx(1).margin(alpha.getError() * 3));
}
