#include <goofit/Catch.h>

#include "testhelpers.h"
#include <goofit/PDFs/basic/CorrGaussianPdf.h>
#include <goofit/UnbinnedDataSet.h>

#include <goofit/Variable.h>

#include <random>

using namespace GooFit;

TEST_CASE("CorrGaussianPdf", "[convert][fit]") {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Observable xvar{"xvar", -10, 10};
    Observable yvar{"yvar", -10, 10};

    // Data set
    UnbinnedDataSet data{{xvar, yvar}};

    // Generate toy events.
    for(int i = 0; i < 10000; ++i) {
        double v = d(gen);
        // double vy = d(gen);

        if(v > 0 && v < 10) {
            xvar.setValue(v);
            yvar.setValue(v);
            data.addEvent();
        }
    }

    // Fit parameter
    Variable mean1{"mean1", 4, 0.01, -10, 10};
    Variable sigma1{"sigma1", 5, 1, 10};
    Variable mean2{"mean2", 4, 0.01, -10, 10};
    Variable sigma2{"sigma2", 5, 1, 10};

    Variable correlation{"correlation", 1, 0.01, 1, 3};

    // GooPdf object
    CorrGaussianPdf corrgausspdf{"corrgausspdf", xvar, yvar, mean1, sigma1, mean2, sigma2, correlation};
    corrgausspdf.setData(&data);

    bool fitter = test_fitter(&corrgausspdf);

    CHECK(fitter);
    // CHECK(mean1.getError() ==  .1);
    // CHECK(mean1.getValue() == Approx(0.665178392).margin(mean1.getError() * 3));
    // CHECK(mean2.getError() ==  .1);
    // CHECK(mean2.getValue() == Approx(0.6837).margin(mean2.getError() * 3));
}
