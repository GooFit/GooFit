#include "gtest/gtest.h"

#include "goofit/fitting/FitManagerMinuit1.h"
#include "goofit/UnbinnedDataSet.h"
#include "goofit/PDFs/ExpPdf.h"
#include "goofit/PDFs/ProdPdf.h"

#include "goofit/Variable.h"

#include <sys/time.h>
#include <sys/times.h>
#include <iostream>

#include <random>

using namespace std;

TEST(Minuit1, SimpleFit) {
    
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Variable xvar{"xvar", 0, 10};
    
    // Data set
    UnbinnedDataSet data(&xvar);
    
    // Generate toy events.
    for(int i=0; i<1000; ++i) {
        double val = d(gen);
        if(val < 10) {
            xvar.value = val;
            data.addEvent();
        }
    }
    
    // Fit parameter
    Variable alpha{"alpha", -2, 0.1, -10, 10};
    
    // GooPdf object
    ExpPdf exppdf{"exppdf", &xvar, &alpha};
    exppdf.setData(&data);
    
    GooFit::FitManagerMinuit1 fitter{&exppdf};
    fitter.setVerbosity(0);
    fitter.fit();

    EXPECT_TRUE(fitter);
    EXPECT_LT(alpha.error, .1);
    EXPECT_NEAR(-1.5, alpha.value, alpha.error*3);
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
    UnbinnedDataSet data {{&xvar, &yvar}, "Some name"};
    
    // Generate toy events.
    for(int i=0; i<20000; ++i) {
        double xval = dx(gen);
        double yval = dy(gen);
        if(xval < 10 && yval < 10) {
            xvar.value = xval;
            yvar.value = yval;
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
    ProdPdf totalpdf {"totalpdf", {&xpdf, &ypdf}};
    totalpdf.setData(&data);
    
    GooFit::FitManagerMinuit1 fitter{&totalpdf};
    fitter.setVerbosity(0);
    fitter.fit();
    
    EXPECT_TRUE(fitter);
    EXPECT_LT(xalpha.error, .1);
    EXPECT_LT(yalpha.error, .1);
    EXPECT_NEAR(-1.5, xalpha.value, xalpha.error*3);
    EXPECT_NEAR(-.75, yalpha.value, yalpha.error*3);
}

