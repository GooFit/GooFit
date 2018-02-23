#include <gtest/gtest.h>

#include <goofit/PDFs/basic/ExpPdf.h>
#include <goofit/UnbinnedDataSet.h>
#include "testhelpers.h"

#include <goofit/Variable.h>





#include <random>

using namespace std;
using namespace GooFit;

TEST(Exp, SimpleFit) {
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

    // GooPdf object
    ExpPdf exppdf{"exppdf", xvar, alpha, sigma};
    exppdf.setData(&data);

    bool fitter =(&exppdf);
    
    

    EXPECT_TRUE(fitter);
    // EXPECT_LT(alpha.getError(), .1);
    // EXPECT_NEAR(0.665178392, alpha.getValue(), alpha.getError() * 3);
}

TEST(ExpPoly, SimpleFit) {
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
    std::vector<Variable> weights;
    Variable alpha{"alpha", 1, 0.1, -10, 10};
    Variable sigma{"sigma", 1, 0, 3};

    weights.push_back(alpha);
    weights.push_back(sigma);

    // Variable offset{"sigma", 1, 0, 3};

    // GooPdf object
    ExpPdf exppdf{"exppdf", xvar, weights};
    exppdf.setData(&data);

    bool fitter = test_fitter(&exppdf);
    
    

    EXPECT_TRUE(fitter);
    // EXPECT_LT(alpha.getError(), .1);
}

TEST(ExpOffset, SimpleFit) {
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

    Variable offset{"offset", 1, 0, 3};

    // GooPdf object
    ExpPdf exppdf{"exppdf", xvar, alpha, offset};
    exppdf.setData(&data);

    bool fitter = test_fitter(&exppdf);
    
    

    EXPECT_TRUE(fitter);
}

TEST(ExpPolyOffset, SimpleFit) {
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
    std::vector<Variable> weights;
    Variable alpha{"alpha", 1, 0.1, -10, 10};
    Variable sigma{"sigma", 1, 0.1, -10, 10};

    weights.push_back(alpha);
    weights.push_back(sigma);

    Variable offset{"offset", 1, 0, 3};

    // GooPdf object
    ExpPdf exppdf{"exppdf", xvar, weights, offset};
    exppdf.setData(&data);

    bool fitter = test_fitter(&exppdf);
    
    

    EXPECT_TRUE(fitter);
}
