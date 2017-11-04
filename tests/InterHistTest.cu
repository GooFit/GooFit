#include "gtest/gtest.h"

#include "goofit/fitting/FitManagerMinuit1.h"
#include "goofit/UnbinnedDataSet.h"
#include "goofit/PDFs/basic/InterHistPdf.h"

#include "goofit/Variable.h"

#include <sys/time.h>
#include <sys/times.h>
#include <iostream>

#include <random>

using namespace std;
using namespace GooFit;

TEST(InterHistPdf, SimpleFit) {
    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Variable xvar{"xvar", -10, 10};

    Variable a1{"a1", 0.1, -10, 10};
    Variable b1{"a1", 0.1, -10, 10};
    Variable c1{"a1", 0.1, -10, 10};

    std::vector<Variable *> var_list;
    var_list.push_back(&a1);
    var_list.push_back(&b1);
    var_list.push_back(&c1);

    std::vector<Variable *> obs_list;
    obs_list.push_back(&xvar);

    // Data set
    BinnedDataSet data(&xvar);

    // Generate toy events.
    for(int i = 0; i < 10000; ++i) {
        double val = d(gen);
        if(val < 10) {
            xvar.setValue(val);
            data.addEvent();
        }
    }

    // string, BinnedDataSet, vector<Variable*> params, vector<Variable*> observables
    InterHistPdf pdf{"InterHistPdf", &data, var_list, obs_list};
    pdf.setData(&data);

    GooFit::FitManagerMinuit1 fitter{&pdf};
    fitter.setVerbosity(2);
    fitter.fit();

    // EXPECT_TRUE(fitter);
    // EXPECT_LT(alpha.getError(), .1);
    // EXPECT_NEAR(0.665178392, alpha.getValue(), alpha.getError() * 3);
}
