#include <gtest/gtest.h>

#include "testhelpers.h"
#include <goofit/PDFs/basic/ArgusPdf.h>
#include <goofit/UnbinnedDataSet.h>

#include <goofit/Variable.h>

#include <random>

using namespace std;
using namespace GooFit;

TEST(ArgusUpper, SimpleFit) {
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

    EXPECT_TRUE(fitter);
    EXPECT_LT(alpha.getError(), .1);
    EXPECT_NEAR(0.998854, alpha.getValue(), alpha.getError() * 3);
}

TEST(ArgusLower, SimpleFit) {
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

    EXPECT_TRUE(fitter);
    EXPECT_LT(alpha.getError(), .1);
    EXPECT_NEAR(1.0012696, alpha.getValue(), alpha.getError() * 3);
}
