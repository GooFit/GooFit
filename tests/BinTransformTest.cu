#include "gtest/gtest.h"

#include "goofit/fitting/FitManagerMinuit1.h"
#include "goofit/BinnedDataSet.h"
#include "goofit/PDFs/basic/BinTransformPdf.h"

#include "goofit/Variable.h"

#include <sys/time.h>
#include <sys/times.h>
#include <iostream>

#include <random>

using namespace std;
using namespace GooFit;

TEST(BinTransformTest, SimpleFit) {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Variable bt_1{"bt_1", 2.0, 0.0, 10.0};
    Variable bt_2{"bt_2", 3.0, 0.0, 10.0};

    std::vector<Variable *> varlist;
    varlist.push_back(&bt_1);
    varlist.push_back(&bt_2);

    // Data set
    BinnedDataSet data(varlist);

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

    // GooPdf object
    BinTransformPdf bintransformpdf{"binlistindex", varlist, list1, list2, listindex};
    bintransformpdf.setData(&data);

    GooFit::FitManagerMinuit1 fitter{&bintransformpdf};
    fitter.setVerbosity(2);
    fitter.fit();

    EXPECT_TRUE(fitter);
    EXPECT_LT(bt_1.getError(), 1.01);
    EXPECT_NEAR(2.01, bt_1.getValue(), bt_1.getError() * 3);
}
