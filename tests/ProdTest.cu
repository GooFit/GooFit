#include <gtest/gtest.h>

#include <goofit/fitting/FitManagerMinuit1.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/PDFs/basic/GaussianPdf.h>
#include <goofit/PDFs/combine/ProdPdf.h>

#include <goofit/Variable.h>

#include <sys/time.h>
#include <sys/times.h>
#include <iostream>

#include <random>

using namespace std;
using namespace GooFit;

TEST(Prod, SimpleFit) {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Observable xvar{"xvar", 1, 10};

    // Data set
    UnbinnedDataSet data(xvar);

    // Generate toy events.
    for(int i = 0; i < 10000; ++i) {
        double val = d(gen);
        if(val < 10 && val > 1) {
            xvar.setValue(val);
            data.addEvent();
        }
    }

    // Fit parameter
    Variable a1{"alpha1", 1, 0.1, -10, 10};
    Variable s1{"sigma1", 1, 0, 3};
    Variable a2{"alpha2", 1, 0.1, -10, 10};
    Variable s2{"sigma2", 1, 0, 3};

    GaussianPdf gauss1{"gauss1", xvar, a1, s1};
    GaussianPdf gauss2{"gauss1", xvar, a2, s2};

    std::vector<PdfBase *> list;
    list.push_back(&gauss1);
    list.push_back(&gauss2);

    // GooPdf object
    ProdPdf pdf{"prodpdf", list};
    pdf.setData(&data);

    GooFit::FitManagerMinuit1 fitter{&pdf};
    fitter.setVerbosity(2);
    fitter.fit();

    EXPECT_TRUE(fitter);
}
