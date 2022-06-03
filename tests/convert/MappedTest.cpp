#include <goofit/Catch.h>

#include "testhelpers.h"
#include <goofit/PDFs/basic/GaussianPdf.h>
#include <goofit/PDFs/combine/MappedPdf.h>
#include <goofit/UnbinnedDataSet.h>

#include <goofit/Variable.h>

#include <random>

using namespace GooFit;

TEST_CASE("Mapped", "[convert][fit]") {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Observable xvar{"xvar", 0, 10};

    // Data set
    UnbinnedDataSet data(xvar);

    // Generate toy events.
    for(int i = 0; i < 100000; ++i) {
        double val = d(gen);
        if(val < 10 && val > 0) {
            xvar.setValue(val);
            data.addEvent();
        }
    }

    // Fit parameter
    Variable a1{"a1", 1, 0.1, 1, 10};
    Variable s1{"s1", 1, 0, 3};
    Variable a2{"a2", 1, 0.1, 1, 10};
    Variable s2{"s2", 1, 0, 3};
    Variable a3{"a3", 1, 0.1, 1, 10};
    Variable s3{"s3", 1, 0, 3};

    GaussianPdf g1{"gauss1", xvar, a1, s1};
    GaussianPdf g2{"gauss2", xvar, a2, s2};
    GaussianPdf g3{"gauss3", xvar, a3, s3};

    std::vector<GooPdf *> list;
    list.push_back(&g2);
    list.push_back(&g3);

    // GooPdf object
    MappedPdf pdf{"mappedpdf", &g1, list};
    pdf.setData(&data);

    bool fitter = test_fitter(&pdf);

    CHECK(fitter);
    // CHECK(alpha.getError() ==  .1);
    // CHECK(alpha.getValue() == Approx(0.665178392).margin(alpha.getError() * 3));
}
