#include <gtest/gtest.h>

#include "testhelpers.h"
#include <goofit/PDFs/basic/ExpGausPdf.h>
#include <goofit/UnbinnedDataSet.h>

#include <goofit/Variable.h>

#include <random>

using namespace std;
using namespace GooFit;

TEST(ExpGausPdf, SimpleFit) {
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

    Variable beta{"beta", 1, 0.1, -10, 10};

    // GooPdf object
    ExpGausPdf expgauspdf{"expgausspdf", xvar, alpha, sigma, beta};
    expgauspdf.setData(&data);

    bool fitter = test_fitter_minuit1(&expgauspdf);

    EXPECT_TRUE(fitter);
    EXPECT_LT(alpha.getError(), 12.84);
    EXPECT_NEAR(12.8379, alpha.getValue(), alpha.getError() * 3);
}
