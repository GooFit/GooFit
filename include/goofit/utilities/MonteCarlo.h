#pragma once

#include <goofit/PDFs/GooPdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>

#include <algorithm>
#include <numeric>
#include <random>

namespace GooFit {

inline void fillDataSetMC1D(GooPdf &pdf, Observable var, size_t nTotal, unsigned int seed = 0) {
    // Setup bins
    UnbinnedDataSet data{var};
    data.fillWithGrid();
    auto origdata = pdf.getData();

    if(origdata == nullptr)
        throw GeneralError("Can't run on a PDF with no DataSet to fill!");

    pdf.setData(&data);
    std::vector<double> pdfValues = pdf.getCompProbsAtDataPoints()[0];

    // Setup random numbers
    if(seed == 0) {
        std::random_device rd;
        seed = rd();
    }
    std::mt19937 gen(seed);

    // Poisson distribution
    std::poisson_distribution<> d(nTotal);
    size_t num_events = d(gen);

    // Uniform distribution
    std::uniform_real_distribution<> unihalf(-.5, .5);
    std::uniform_real_distribution<> uniwhole(0.0, 1.0);

    // CumSum in other languages
    std::vector<double> integral(pdfValues.size());
    std::partial_sum(pdfValues.begin(), pdfValues.end(), integral.begin());

    // Make this a 0-1 fraction by dividing by the end value
    std::for_each(integral.begin(), integral.end(), [&integral](double &val) { val /= integral.back(); });

    for(size_t i = 0; i < num_events; i++) {
        double r = uniwhole(gen);

        // Binary search for integral[cell-1] < r < integral[cell]
        size_t j = std::lower_bound(integral.begin(), integral.end(), r) - integral.begin();

        // Fill in the grid randomly
        double varValue = data.getValue(var, j) + var.getBinSize() * unihalf(gen);

        var.setValue(varValue);
        origdata->addEvent();
    }

    pdf.setData(origdata);
}

} // namespace GooFit
