#include <gtest/gtest.h>

#include <goofit/PDFs/basic/KinLimitBWPdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/fitting/FitManagerMinuit1.h>

#include <goofit/Variable.h>

#include <iostream>
#include <sys/time.h>
#include <sys/times.h>

#include <random>

using namespace std;
using namespace GooFit;

TEST(Gaussian, SimpleFit) {
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
    Variable alpha{"alpha", 1, 0.1, 1, 10};
    Variable sigma{"sigma", 1, 1, 3};

    // GooPdf object
    KinLimitBWPdf kinlimitbwpdf{"kinlimitbwpdf", xvar, alpha, sigma};
    kinlimitbwpdf.setData(&data);

    GooFit::FitManagerMinuit1 fitter{&kinlimitbwpdf};
    fitter.setVerbosity(2);
    fitter.fit();

    EXPECT_TRUE(fitter);
    EXPECT_LT(alpha.getError(), .1);
    EXPECT_NEAR(1, alpha.getValue(), alpha.getError() * 3);
}
