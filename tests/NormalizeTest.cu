#include "gtest/gtest.h"

#include "goofit/FitManager.h"
#include "goofit/UnbinnedDataSet.h"
#include "goofit/PDFs/basic/ExpPdf.h"
#include "goofit/PDFs/combine/ProdPdf.h"

#include "goofit/Variable.h"

#include <sys/time.h>
#include <sys/times.h>
#include <iostream>

#include <random>

using namespace GooFit;

TEST(Normalize, Dual) {
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

    std::cout << exppdf.normalize() << std::endl;
    std::cout << exppdf.normalize() << std::endl;

    FitManager fitter{&exppdf};
    fitter.setVerbosity(0);
    fitter.fit();

    std::cout << exppdf.normalize() << std::endl;
    std::cout << exppdf.normalise() << std::endl;

    std::cout << exppdf.integrate(0, 10) << std::endl;
}
