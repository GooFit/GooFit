#include <gtest/gtest.h>

#include <random>

#include "goofit/Variable.h"
#include "goofit/UnbinnedDataSet.h"
#include "goofit/PDFs/basic/GaussianPdf.h"
#include "goofit/FitManager.h"

#ifdef ROOT_FOUND
#include "goofit/fitting/FitManagerMinuit1.h"
#endif

using namespace GooFit;

TEST(Simple, NoBlind) {
    // Independent variable.
    Variable xvar{"xvar", 0, 10};

    // Data set
    UnbinnedDataSet data{&xvar};

    // Random number generation
    std::mt19937 gen(137);
    std::normal_distribution<> d(1.5, .3);

    // Generate toy events.
    for(int i = 0; i < 100000; ++i) {
        double val = d(gen);
        if(std::fabs(val) < 10) {
            xvar.setValue(val);
            data.addEvent();
        }
    }

    // Fit parameter
    Variable alpha{"alpha", -2, -10, 10};
    Variable sigma{"sigma", 1, 0, 3};

    // GooPdf object
    GaussianPdf gausspdf{"gausspdf", &xvar, &alpha, &sigma};
    gausspdf.fitTo(&data);

    EXPECT_LT(alpha.getError(), .01);
    EXPECT_NEAR(1.5, alpha.getValue(), alpha.getError() * 3);
}

TEST(Simple, WithBlind) {
    // Independent variable.
    Variable xvar{"xvar", 0, 10};

    // Data set
    UnbinnedDataSet data{&xvar};

    // Random number generation
    std::mt19937 gen(137);
    std::normal_distribution<> d(1.5, .3);

    // Generate toy events.
    for(int i = 0; i < 100000; ++i) {
        double val = d(gen);
        if(std::fabs(val) < 10) {
            xvar.setValue(val);
            data.addEvent();
        }
    }

    // Fit parameter
    Variable alpha{"alpha", -2, -10, 10};
    Variable sigma{"sigma", 1, 0, 3};

    // Blinding (normally would be randomly generated)
    alpha.setBlind(1);

    // GooPdf object
    GaussianPdf gausspdf{"gausspdf", &xvar, &alpha, &sigma};
    gausspdf.fitTo(&data);

    EXPECT_LT(alpha.getError(), .01);
    EXPECT_NEAR(2.5, alpha.getValue(), alpha.getError() * 3);
}
#ifdef ROOT_FOUND
TEST(Simple, Min1Blind) {
    // Independent variable.
    Variable xvar{"xvar", 0, 10};

    // Data set
    UnbinnedDataSet data{&xvar};

    // Random number generation
    std::mt19937 gen(137);
    std::normal_distribution<> d(1.5, .3);

    // Generate toy events.
    for(int i = 0; i < 100000; ++i) {
        double val = d(gen);
        if(std::fabs(val) < 10) {
            xvar.setValue(val);
            data.addEvent();
        }
    }

    // Fit parameter
    Variable alpha{"alpha", -2, -10, 10};
    Variable sigma{"sigma", 1, 0, 3};

    // Blinding (normally would be randomly generated)
    alpha.setBlind(1);

    // GooPdf object
    GaussianPdf gausspdf{"gausspdf", &xvar, &alpha, &sigma};
    gausspdf.setData(&data);

    GooFit::FitManagerMinuit1 fitman{&gausspdf};
    fitman.fit();

    EXPECT_LT(alpha.getError(), .01);
    EXPECT_NEAR(2.5, alpha.getValue(), alpha.getError() * 3);
}
#endif
