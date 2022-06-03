#include <goofit/Catch.h>

#include "testhelpers.h"
#include <goofit/PDFs/basic/GaussianPdf.h>
#include <goofit/PDFs/combine/EventWeightedAddPdf.h>
#include <goofit/UnbinnedDataSet.h>

#include <goofit/Variable.h>

#include <random>

using namespace GooFit;

TEST_CASE("EventWeightedAdd", "[convert][fit]") {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Observable w1{"w1", -10, 10};
    Observable w2{"w2", -10, 10};

    std::vector<Observable> weights;
    weights.push_back(w1);
    weights.push_back(w2);

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
    Variable s1{"s1", 1, 0, 5};
    GaussianPdf gauss1{"gauss1", w1, a1, s1};

    Variable a2{"a2", 1, 0.1, -10, 10};
    Variable s2{"s2", 1, 0, 5};
    GaussianPdf gauss2{"gauss2", w2, a2, s2};

    std::vector<PdfBase *> pdfs;
    pdfs.push_back(&gauss1);
    pdfs.push_back(&gauss2);

    // GooPdf object
    EventWeightedAddPdf eventpdf{"eventweightedaddpdf", weights, pdfs};
    eventpdf.setData(&data);

    bool fitter = test_fitter(&eventpdf);

    CHECK(fitter);
    // CHECK(w1.getError() ==  .11);
    // CHECK(w1.getValue() == Approx(0.1).margin(w1.getError() * 3));
    // CHECK(w2.getError() ==  .11);
    // CHECK(w2.getValue() == Approx(0.1).margin(w2.getError() * 3));
}
