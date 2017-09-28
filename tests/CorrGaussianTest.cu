#include "gtest/gtest.h"

#include "goofit/fitting/FitManagerMinuit1.h"
#include "goofit/UnbinnedDataSet.h"
#include "goofit/PDFs/basic/CorrGaussianPdf.h"

#include "goofit/Variable.h"

#include <sys/time.h>
#include <sys/times.h>
#include <iostream>

#include <random>

using namespace std;
using namespace GooFit;

TEST(CorrGaussianPdf, SimpleFit) {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Variable xvar{"xvar", -10, 10};
    Variable yvar{"yvar", -10, 10};

    std::vector<Variable *> list;
    list.push_back(&xvar);
    list.push_back(&yvar);

    // Data set
    UnbinnedDataSet data(list);

    // Generate toy events.
    for(int i = 0; i < 1000; ++i) {
        double vx = d(gen);
        double vy = d(gen);
        
        //if(val < 10) {
            xvar.setValue(vx);
            yvar.setValue(vy);
            data.addEvent();
        //}
    }

    // Fit parameter
    Variable mean1{"mean1", 1, 0.1, -10, 10};
    Variable sigma1{"sigma1", 1, 0, 3};
    Variable mean2{"mean2", 1, 0.1, -10, 10};
    Variable sigma2{"sigma2", 1, 0, 3};

    Variable correlation {"correlation", 1, 0, 3};

    // GooPdf object
    CorrGaussianPdf corrgausspdf{"corrgausspdf", &xvar, &yvar, &mean1, &sigma1, &mean2, &sigma2, &correlation};
    corrgausspdf.setData(&data);

    GooFit::FitManagerMinuit1 fitter{&corrgausspdf};
    fitter.setVerbosity(2);
    fitter.fit();

    EXPECT_TRUE(fitter);
    EXPECT_LT(mean1.getError(), .1);
    EXPECT_NEAR(0.665178392, mean1.getValue(), mean1.getError() * 3);
    EXPECT_LT(mean2.getError(), .1);
    EXPECT_NEAR(0.6837, mean2.getValue(), mean2.getError() * 3);
}

