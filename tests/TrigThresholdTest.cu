#include "gtest/gtest.h"

#include "goofit/fitting/FitManagerMinuit1.h"
#include "goofit/UnbinnedDataSet.h"
#include "goofit/PDFs/basic/TrigThresholdPdf.h"

#include "goofit/Variable.h"

#include <sys/time.h>
#include <sys/times.h>
#include <iostream>

#include <random>

using namespace std;
using namespace GooFit;

TEST(TrigThresholdUpper, SimpleFit) {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Variable xvar{"xvar", -10, 10};

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
    Variable thres{"m", 1, 0.1, -10, 10};
    Variable trig{"s", 1, 0, 3};
    Variable linCons{"w", 1, 0, 3};

    // GooPdf object
    TrigThresholdPdf pdf{"trigthresholdpdf", &xvar, &thres, &trig, &linCons, true};
    pdf.setData(&data);

    GooFit::FitManagerMinuit1 fitter{&pdf};
    fitter.setVerbosity(2);
    fitter.fit();

    EXPECT_TRUE(fitter);
}

TEST(TrigThresholdLower, SimpleFit) {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Variable xvar{"xvar", -10, 10};

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
    Variable thres{"m", 1, 0.1, -10, 10};
    Variable trig{"s", 1, 0, 3};
    Variable linCons{"w", 1, 0, 3};

    // GooPdf object
    TrigThresholdPdf pdf{"trigthresholdpdf", &xvar, &thres, &trig, &linCons, false};
    pdf.setData(&data);

    GooFit::FitManagerMinuit1 fitter{&pdf};
    fitter.setVerbosity(2);
    fitter.fit();

    EXPECT_TRUE(fitter);
}

TEST(TrigThresholdVeryUpper, SimpleFit) {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Variable xvar{"xvar", -10, 10};
    Variable yvar{"yvar", -10, 10};

    std::vector<Variable*> list;
    list.push_back(&xvar);
    list.push_back(&yvar);

    // Data set
    UnbinnedDataSet data(list);

    // Generate toy events.
    for(int i = 0; i < 100000; ++i) {
        double val = d(gen);
        if(val < 10) {
            xvar.setValue(val);
            yvar.setValue(val);
            data.addEvent();
        }
    }

    // Fit parameter
    Variable thres{"m", 1, 0.1, 1, 10};
    Variable trig{"s", 1, 0, 3};
    Variable linCons{"w", 1, 0, 3};
    Variable mass{"mass", 1, 0, 3};

    // GooPdf object
    TrigThresholdPdf pdf{"trigthresholdpdf", &xvar, &yvar, &thres, &trig, &linCons, &mass, true};
    pdf.setData(&data);

    GooFit::FitManagerMinuit1 fitter{&pdf};
    fitter.setVerbosity(2);
    fitter.fit();

    EXPECT_TRUE(fitter);
}

TEST(TrigThresholdVeryLower, SimpleFit) {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Variable xvar{"xvar", -10, 10};
    Variable yvar{"yvar", -10, 10};

    std::vector<Variable*> list;
    list.push_back (&xvar);
    list.push_back (&yvar);

    // Data set
    UnbinnedDataSet data(list);

    // Generate toy events.
    for(int i = 0; i < 100000; ++i) {
        double val = d(gen);
        if(val < 10) {
            xvar.setValue(val);
            yvar.setValue(val);
            data.addEvent();
        }
    }

    // Fit parameter
    Variable thres{"m", 1, 0.1, -10, 10};
    Variable trig{"s", 1, 0, 3};
    Variable linCons{"w", 1, 0, 3};
    Variable mass{"mass", 1, 0, 3};

    // GooPdf object
    TrigThresholdPdf pdf{"trigthresholdpdf", &xvar, &yvar, &thres, &trig, &linCons, &mass, false};
    pdf.setData(&data);

    GooFit::FitManagerMinuit1 fitter{&pdf};
    fitter.setVerbosity(2);
    fitter.fit();

    EXPECT_TRUE(fitter);
}

