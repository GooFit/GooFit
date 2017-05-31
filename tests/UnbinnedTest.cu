#include "gtest/gtest.h"

#include "goofit/FitManager.h"
#include "goofit/BinnedDataSet.h"
#include "goofit/PDFs/basic/ExpPdf.h"
#include "goofit/PDFs/basic/GaussianPdf.h"
#include "goofit/PDFs/combine/ProdPdf.h"

#include "goofit/Variable.h"

#include <sys/time.h>
#include <sys/times.h>
#include <iostream>

#include <random>

using namespace GooFit;

TEST(UnbinnedFit, SimpleFit) {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Variable xvar{"xvar", 0, 10};

    // Data set
    UnbinnedDataSet data(&xvar);

    // Generate toy events.
    for(int i = 0; i < 1000; ++i) {
        double val = d(gen);
        if(val < 10) {
            xvar.setValue(val);
            data.addEvent();
        }
    }

    // Fit parameter
    Variable alpha{"alpha", -2, 0.1, -10, 10};

    // GooPdf object
    ExpPdf exppdf{"exppdf", &xvar, &alpha};
    exppdf.setData(&data);

    FitManager fitter{&exppdf};
    fitter.fit();

    EXPECT_TRUE(fitter);
    EXPECT_LT(alpha.getError(), .1);
    EXPECT_NEAR(-1.5, alpha.getValue(), alpha.getError() * 3);
}

TEST(UnbinnedFit, DualFit) {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> dx(1.5);
    std::exponential_distribution<> dy(.75);

    // Independent variable.
    Variable xvar{"xvar", 0, 10};
    Variable yvar{"yvar", 0, 10};

    // Data set
    UnbinnedDataSet data{{&xvar, &yvar}};

    // Generate toy events.
    for(int i = 0; i < 20000; ++i) {
        double xval = dx(gen);
        double yval = dy(gen);
        if(xval < 10 && yval < 10) {
            xvar.setValue(xval);
            yvar.setValue(yval);
            data.addEvent();
        }
    }

    // Fit parameter
    Variable xalpha{"xalpha", -2, 0.1, -10, 10};
    // Fit parameter
    Variable yalpha{"yalpha", -2, 0.1, -10, 10};

    // GooPdf object
    ExpPdf xpdf{"xpdf", &xvar, &xalpha};
    ExpPdf ypdf{"ypdf", &yvar, &yalpha};
    ProdPdf totalpdf{"totalpdf", {&xpdf, &ypdf}};
    totalpdf.setData(&data);

    FitManager fitter{&totalpdf};
    fitter.fit();

    EXPECT_TRUE(fitter);
    EXPECT_LT(xalpha.getError(), .1);
    EXPECT_LT(yalpha.getError(), .1);
    EXPECT_NEAR(-1.5, xalpha.getValue(), xalpha.getError() * 3);
    EXPECT_NEAR(-.75, yalpha.getValue(), yalpha.getError() * 3);
}

TEST(UnbinnedFit, DifferentFitterVariable) {
    // Random number generation
    std::mt19937 gen(137);
    std::normal_distribution<> dx(1.5, .3);
    std::normal_distribution<> dy(.75, .5);

    // Independent variable.
    Variable xvar{"xvar", -10, 10};
    Variable yvar{"yvar", -10, 10};

    // Data set
    UnbinnedDataSet data{{&xvar, &yvar}, "Some name"};

    // Generate toy events.
    for(int i = 0; i < 20000; ++i) {
        xvar = dx(gen);
        yvar = dy(gen);
        if(xvar && yvar)
            data.addEvent();
    }

    // Fit parameter
    Variable xalpha{"xalpha", 0, 0.1, -10, 10};
    Variable xsigma{"xsigma", .1, 0.1, 0, 3};
    Variable unused{"unused", 1, .1, 0, 2};
    // Fit parameter
    Variable ysigma{"ysigma", .1, 0.1, 0, 3};
    Variable yalpha{"yalpha", 0, 0.1, -10, 10};

    // GooPdf object
    GaussianPdf ypdf{"ypdf", &yvar, &yalpha, &ysigma};
    GaussianPdf xpdf{"xpdf", &xvar, &xalpha, &xsigma};
    ProdPdf totalpdf{"totalpdf", {&xpdf, &ypdf}};
    totalpdf.setData(&data);

    FitManager fitter{&totalpdf};
    fitter.fit();

    EXPECT_TRUE(fitter);
    EXPECT_LT(xalpha.getError(), .1);
    EXPECT_LT(yalpha.getError(), .1);
    EXPECT_LT(xsigma.getError(), .1);
    EXPECT_LT(ysigma.getError(), .1);

    EXPECT_NEAR(1.5, xalpha.getValue(), xalpha.getError() * 3);
    EXPECT_NEAR(.75, yalpha.getValue(), yalpha.getError() * 3);
    EXPECT_NEAR(.3, xsigma.getValue(), xsigma.getError() * 3);
    EXPECT_NEAR(.5, ysigma.getValue(), ysigma.getError() * 3);
}
