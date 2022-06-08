#include <goofit/Catch.h>

#include "testhelpers.h"
#include <goofit/PDFs/basic/LandauPdf.h>
#include <goofit/UnbinnedDataSet.h>

#include <goofit/Variable.h>

#include <random>

using namespace GooFit;

TEST_CASE("Landau", "[convert][fit]") {
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
    Variable alpha{"alpha", 1, 0.1, -10, 10};

    Variable sigma{"sigma", 1, 0, 3};

    // GooPdf object
    LandauPdf landaupdf{"landaupdf", xvar, alpha, sigma};
    landaupdf.setData(&data);

    bool fitter = test_fitter(&landaupdf);

    CHECK(fitter);
    CHECK(alpha.getError() < .1);
    CHECK(alpha.getValue() == Approx(0.259).margin(alpha.getError() * 3));
}
