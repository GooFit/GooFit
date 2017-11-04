#include "gtest/gtest.h"

#include "goofit/fitting/FitManagerMinuit1.h"
#include "goofit/UnbinnedDataSet.h"
#include "goofit/PDFs/combine/EventWeightedAddPdf.h"
#include "goofit/PDFs/basic/GaussianPdf.h"

#include "goofit/Variable.h"

#include <sys/time.h>
#include <sys/times.h>
#include <iostream>

#include <random>

using namespace std;
using namespace GooFit;

TEST(EventWeightedAdd, SimpleFit) {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Variable w1{"w1", 1, 0.1, -10, 10};
    Variable w2{"w2", 1, 0.1, -10, 10};

    std::vector<Variable *> weights;
    weights.push_back(&w1);
    weights.push_back(&w2);

    // Data set
    UnbinnedDataSet data(weights);

    // Generate toy events.
    for(int i = 0; i < 1000; ++i) {
        double val = d(gen);
        if(val < 10) {
            w1.setValue(val);
            w2.setValue(val);
            data.addEvent();
        }
    }

    // Fit parameter
    Variable a1{"a1", 1, 0.1, -10, 10};
    Variable s1{"s1", 1, 0, 3};
    GaussianPdf gauss1{"gauss1", &w1, &a1, &s1};

    Variable a2{"a2", 1, 0.1, -10, 10};
    Variable s2{"s2", 1, 0, 3};
    GaussianPdf gauss2{"gauss2", &w2, &a2, &s2};

    std::vector<PdfBase *> pdfs;
    pdfs.push_back(&gauss1);
    pdfs.push_back(&gauss2);

    // GooPdf object
    EventWeightedAddPdf eventpdf{"eventweightedaddpdf", weights, pdfs};
    eventpdf.setData(&data);

    GooFit::FitManagerMinuit1 fitter{&eventpdf};
    fitter.setVerbosity(2);
    fitter.fit();

    EXPECT_TRUE(fitter);
    EXPECT_LT(w1.getError(), .11);
    EXPECT_NEAR(0.1, w1.getValue(), w1.getError() * 3);
    EXPECT_LT(w2.getError(), .11);
    EXPECT_NEAR(0.1, w2.getValue(), w2.getError() * 3);
}
