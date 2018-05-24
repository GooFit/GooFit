#include <gtest/gtest.h>

#include "testhelpers.h"
#include <goofit/PDFs/basic/ScaledGaussianPdf.h>
#include <goofit/UnbinnedDataSet.h>

#include <goofit/Variable.h>

#include <random>

using namespace std;
using namespace GooFit;

TEST(ScaledGaussian, SimpleFit) {
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

    Variable sigma{"sigma", 1, 1, 3};

    Variable scale1{"scale1", 1, 1, 5};
    Variable scale2{"scale2", 1, 1, 5};

    // GooPdf object
    ScaledGaussianPdf scaledgausspdf{"scaledgausspdf", xvar, alpha, sigma, scale1, scale2};
    scaledgausspdf.setData(&data);

    bool fitter = test_fitter(&scaledgausspdf);

    EXPECT_TRUE(fitter);
    EXPECT_LT(alpha.getError(), .1);
    EXPECT_NEAR(-0.33, alpha.getValue(), alpha.getError() * 3);
}
