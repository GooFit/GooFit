#include <goofit/Catch.h>

#include "testhelpers.h"
#include <goofit/PDFs/basic/GaussianPdf.h>
#include <goofit/PDFs/combine/ProdPdf.h>
#include <goofit/UnbinnedDataSet.h>

#include <goofit/Variable.h>

#include <random>

using namespace GooFit;

TEST_CASE("Prod", "[convert][fit]") {
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

    bool fitter = test_fitter(&pdf);

    CHECK(fitter);
}
