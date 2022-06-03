#include <goofit/Catch.h>

#include "testhelpers.h"
#include <goofit/PDFs/basic/ArgusPdf.h>
#include <goofit/UnbinnedDataSet.h>

#include <goofit/Variable.h>

#include <random>

using namespace GooFit;

TEST_CASE("ArgusUpper", "[convert][fit]") {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Observable xvar{"xvar", 0, 10};

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

    Variable beta{"beta", 1, 0.1, -10, 10};

    // GooPdf object
    ArgusPdf arguspdf{"Arguspdf", xvar, alpha, beta, false};
    arguspdf.setData(&data);

    bool fitter = test_fitter_minuit1(&arguspdf);

    CHECK(fitter);
    CHECK(alpha.getError() < .1);
    CHECK(alpha.getValue() == Approx(0.998854).margin(alpha.getError() * 3));
}

TEST_CASE("ArgusLower", "[convert][fit]") {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Observable xvar{"xvar", 0, 10};

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

    Variable beta{"beta", 1, 0.1, -10, 10};

    // GooPdf object
    ArgusPdf arguspdf{"Arguspdf", xvar, alpha, beta, true};
    arguspdf.setData(&data);

    bool fitter = test_fitter_minuit1(&arguspdf);

    CHECK(fitter);
    CHECK(alpha.getError() < .1);
    CHECK(alpha.getValue() == Approx(1.0012696).margin(alpha.getError() * 3));
}
