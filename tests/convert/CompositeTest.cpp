#include <gtest/gtest.h>

#include "testhelpers.h"
#include <goofit/PDFs/basic/GaussianPdf.h>
#include <goofit/PDFs/combine/CompositePdf.h>
#include <goofit/UnbinnedDataSet.h>

#include <goofit/Variable.h>

#include <random>

using namespace std;
using namespace GooFit;

TEST(Composite, SimpleFit) {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Observable xvar{"xvar", 1, 10};

    // Data set
    UnbinnedDataSet data(xvar);

    // Generate toy events.
    for(int i = 0; i < 1000; ++i) {
        double val = d(gen);
        if(val < 10 && val > 1) {
            xvar.setValue(val);
            data.addEvent();
        }
    }

    // Fit parameter
    Variable alpha1{"alpha1", 2, 0.1, 1, 10};
    Variable sigma1{"sigma1", 2, 1, 3};
    GaussianPdf gauss1{"gauss1", xvar, alpha1, sigma1};

    Variable alpha2{"alpha2", 2, 0.1, 1, 10};
    Variable sigma2{"sigma2", 2, 1, 3};
    GaussianPdf gauss2{"gauss2", xvar, alpha1, sigma1};

    // GooPdf object
    CompositePdf compositepdf{"compositepdf", &gauss1, &gauss2};
    compositepdf.setData(&data);

    bool fitter = test_fitter(&compositepdf);

    EXPECT_TRUE(fitter);
    EXPECT_LT(alpha1.getError(), .1);
    EXPECT_NEAR(1.8, alpha1.getValue(), alpha1.getError() * 3);
}
