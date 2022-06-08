#include <goofit/Catch.h>

#include "testhelpers.h"
#include <goofit/BinnedDataSet.h>
#include <goofit/PDFs/basic/BinTransformPdf.h>

#include <goofit/Variable.h>

#include <random>

using namespace GooFit;

TEST_CASE("BinTransformTest", "[convert][fit]") {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Observable bt_1{"bt_1", 0.0, 10.0};
    Observable bt_2{"bt_2", 0.0, 10.0};

    // Data set
    BinnedDataSet data{{bt_1, bt_2}};

    // Generate toy events.
    for(int i = 0; i < 100000; ++i) {
        double val = d(gen);
        if(val < 10) {
            bt_1.setValue(val);
            bt_2.setValue(val);
            data.addEvent();
        }
    }

    // Fit parameter
    std::vector<double> list1;
    list1.push_back(0.343265);
    list1.push_back(0.497609);

    std::vector<double> list2;
    list2.push_back(0.596078);
    list2.push_back(0.123007);

    std::vector<int> listindex;
    listindex.push_back(100);
    listindex.push_back(100);

    std::vector<Observable> obslist;
    obslist.push_back(bt_1);
    obslist.push_back(bt_2);

    // GooPdf object
    BinTransformPdf bintransformpdf{"binlistindex", obslist, list1, list2, listindex};
    bintransformpdf.setData(&data);

    bool fitter = test_fitter_minuit1(&bintransformpdf);

    CHECK(fitter);
    // CHECK(bt_1.getError() ==  1.01);
    // CHECK(bt_1.getValue() == Approx(2.01).margin(bt_1.getError() * 3));
}
