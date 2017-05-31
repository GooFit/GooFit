#include "gtest/gtest.h"

#include "goofit/fitting/FitManagerMinuit1.h"
#include "goofit/UnbinnedDataSet.h"
#include "goofit/PDFs/basic/ExpPdf.h"
#include "goofit/PDFs/basic/GaussianPdf.h"
#include "goofit/PDFs/combine/ProdPdf.h"

#include "goofit/Variable.h"

#include <sys/time.h>
#include <sys/times.h>
#include <iostream>

#include <random>

using namespace std;
using namespace GooFit;

TEST(Minuit1, SimpleFit) {
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

    GooFit::FitManagerMinuit1 fitter{&exppdf};
    fitter.setVerbosity(2);
    fitter.fit();

    EXPECT_TRUE(fitter);
    EXPECT_LT(alpha.getError(), .1);
    EXPECT_NEAR(-1.5, alpha.getValue(), alpha.getError() * 3);
}

TEST(Minuit1, DualFit) {
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

    GooFit::FitManagerMinuit1 fitter{&totalpdf};
    fitter.setVerbosity(0);
    fitter.fit();

    EXPECT_TRUE(fitter);
    EXPECT_LT(xalpha.getError(), .1);
    EXPECT_LT(yalpha.getError(), .1);
    EXPECT_NEAR(-1.5, xalpha.getValue(), xalpha.getError() * 3);
    EXPECT_NEAR(-.75, yalpha.getValue(), yalpha.getError() * 3);
}

TEST(Minuit1, DifferentFitterVariable) {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> dx(1.5);
    std::exponential_distribution<> dy(.75);

    // Independent variable.
    Variable xvar{"xvar", 0, 10};
    Variable yvar{"yvar", 0, 10};

    // Data set
    UnbinnedDataSet data{{&xvar, &yvar}, "Some name"};

    // Generate toy events.
    for(int i = 0; i < 20000; ++i) {
        xvar = dx(gen);
        yvar = dy(gen);
        if(xvar && yvar) {
            data.addEvent();
        }
    }

    // Fit parameter
    Variable xalpha{"xalpha", -2, 0.1, -10, 10};
    // Fit parameter
    Variable yalpha{"yalpha", -2, 0.1, -10, 10};

    // GooPdf object
    ExpPdf ypdf{"ypdf", &yvar, &yalpha};
    ExpPdf xpdf{"xpdf", &xvar, &xalpha};
    ProdPdf totalpdf{"totalpdf", {&xpdf, &ypdf}};
    totalpdf.setData(&data);

    GooFit::FitManagerMinuit1 fitter{&totalpdf};
    fitter.setVerbosity(2);
    fitter.fit();

    EXPECT_TRUE(fitter);
    EXPECT_LT(xalpha.getError(), .1);
    EXPECT_LT(yalpha.getError(), .1);
    EXPECT_NEAR(-1.5, xalpha.getValue(), xalpha.getError() * 3);
    EXPECT_NEAR(-.75, yalpha.getValue(), yalpha.getError() * 3);
}

TEST(Minuit1, FitterConstants) {
    // Random number generation
    std::mt19937 gen(137);
    std::normal_distribution<> dx(1.5, .3);
    std::normal_distribution<> dy(-.75, .2);

    // Independent variable.
    Variable xvar{"xvar", -5, 5};
    Variable yvar{"yvar", -5, 5};

    // Data set
    UnbinnedDataSet data{{&xvar, &yvar}, "Some name"};

    // Generate toy events.
    for(int i = 0; i < 20000; ++i) {
        double xval = dx(gen);
        double yval = dy(gen);
        try {
            xvar.setValue(xval);
            yvar.setValue(yval);
            data.addEvent();
        } catch(const GooFit::OutOfRange &) {
        }
    }

    // Fit parameter
    Variable xalpha{"xalpha", 2, 0.1, -10, 10};
    Variable xsigma{"xsigma", .2};
    Variable yalpha{"yalpha", -1, 0.1, -10, 10};
    Variable ysigma{"ysigma", .3};

    // GooPdf object
    GaussianPdf xpdf{"xpdf", &xvar, &xalpha, &xsigma};
    GaussianPdf ypdf{"ypdf", &yvar, &yalpha, &ysigma};
    ProdPdf totalpdf{"totalpdf", {&xpdf, &ypdf}};
    totalpdf.setData(&data);

    GooFit::FitManagerMinuit1 fitter{&totalpdf};
    fitter.setVerbosity(2);
    fitter.fit();

    EXPECT_TRUE(fitter);
    EXPECT_LT(xalpha.getError(), .1);
    EXPECT_LT(yalpha.getError(), .1);
    EXPECT_NEAR(1.5, xalpha.getValue(), xalpha.getError() * 3);
    EXPECT_NEAR(-.75, yalpha.getValue(), yalpha.getError() * 3);
    EXPECT_EQ(.2, xsigma.getValue());
    EXPECT_EQ(.3, ysigma.getValue());
}
