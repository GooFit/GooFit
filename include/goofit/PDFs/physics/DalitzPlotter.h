#pragma once

#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/physics/Amp3Body.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/Version.h>

#include <algorithm>
#include <numeric>
#include <random>

#if GOOFIT_ROOT_FOUND
#include <TH2F.h>
#endif

namespace GooFit {

/// This class makes it easy to make plots over 3 body Dalitz PDFs. You can use ROOT style value access or bin numbers.
class DalitzPlotter {
    std::vector<size_t> xbins;
    std::vector<size_t> ybins;
    std::vector<std::vector<fptype>> pdfValues;
    Observable m13;
    Observable m23;
    EventNumber eventNumber;
    UnbinnedDataSet data;
    fptype mother;

  public:
    DalitzPlotter(GooPdf *overallSignal, Amp3Body *signalDalitz)
        : m13(signalDalitz->_m13)
        , m23(signalDalitz->_m23)
        , eventNumber(signalDalitz->_eventNumber)
        , data({m13, m23, eventNumber})
        , mother(signalDalitz->decayInfo.motherMass) {
        eventNumber.setValue(0);

        for(size_t i = 0; i < m13.getNumBins(); ++i) {
            m13.setValue(m13.getLowerLimit() + m13.getBinSize() * (i + 0.5));
            for(size_t j = 0; j < m23.getNumBins(); ++j) {
                m23.setValue(m23.getLowerLimit() + m23.getBinSize() * (j + 0.5));
                if(inDalitz2(m13.getValue(),
                            m23.getValue(),
                            signalDalitz->decayInfo.motherMass,
                            signalDalitz->decayInfo.daug1Mass,
                            signalDalitz->decayInfo.daug2Mass,
                            signalDalitz->decayInfo.daug3Mass)) {
                    xbins.push_back(i);
                    ybins.push_back(j);
                    data.addEvent();
                    eventNumber.setValue(eventNumber.getValue() + 1);
                }
            }
        }

        auto old = overallSignal->getData();
        overallSignal->setData(&data);
        signalDalitz->setDataSize(data.getNumEvents());
        pdfValues = overallSignal->getCompProbsAtDataPoints();
        overallSignal->setData(old);
    }

    /// Fill a dataset with MC events
    void fillDataSetMC(UnbinnedDataSet &dataset, size_t nTotal) {
        // Setup random numbers
        std::random_device rd;
        std::mt19937 gen(rd());

        // Uniform distribution
        std::uniform_real_distribution<> unihalf(-.5, .5);
        std::uniform_real_distribution<> uniwhole(0.0, 1.0);

        // CumSum in other languages
        std::vector<double> integral(pdfValues[0].size());
        std::partial_sum(pdfValues[0].begin(), pdfValues[0].end(), integral.begin());

        // Make this a 0-1 fraction by dividing by the end value
        std::for_each(integral.begin(), integral.end(), [&integral](double &val) { val /= integral.back(); });

        for(size_t i = 0; i < nTotal; i++) {
            double r = uniwhole(gen);

            // Binary search for integral[cell-1] < r < integral[cell]
            size_t j = std::lower_bound(integral.begin(), integral.end(), r) - integral.begin();

            // Fill in the grid square randomly
            double currm13 = data.getValue(m13, j) + m13.getBinSize() * unihalf(gen);
            double currm23 = data.getValue(m23, j) + m23.getBinSize() * unihalf(gen);

            m13.setValue(currm13);
            m23.setValue(currm23);
            eventNumber.setValue(i);
            dataset.addEvent();
        }
    }

    auto getNumEvents() const -> size_t { return data.getNumEvents(); }

    auto getX(size_t event) const -> size_t { return xbins.at(event); }

    auto getY(size_t event) const -> size_t { return ybins.at(event); }

    auto getXval(size_t event) const -> fptype { return data.getValue(m13, event); }

    auto getYval(size_t event) const -> fptype { return data.getValue(m23, event); }

    auto getZval(size_t event) const -> fptype { return POW2(mother) - POW2(getXval(event)) - POW2(getYval(event)); }

    auto getVal(size_t event, size_t num = 0) const -> fptype { return pdfValues.at(num).at(event); }

    auto getDataSet() -> UnbinnedDataSet * { return &data; }

    auto getM13() const -> const Observable & { return m13; }
    auto getM23() const -> const Observable & { return m23; }

#if GOOFIT_ROOT_FOUND
    /// Produce a TH2F over the contained evaluation
    auto make2D(std::string name = "dalitzplot", std::string title = "") -> TH2F * {
        auto *dalitzplot = new TH2F(name.c_str(),
                                    title.c_str(),
                                    m13.getNumBins(),
                                    m13.getLowerLimit(),
                                    m13.getUpperLimit(),
                                    m23.getNumBins(),
                                    m23.getLowerLimit(),
                                    m23.getUpperLimit());

        for(unsigned int j = 0; j < getNumEvents(); ++j) {
            size_t currm13 = getX(j);
            size_t currm23 = getY(j);
            double val     = getVal(j);

            dalitzplot->SetBinContent(1 + currm13, 1 + currm23, val);
        }

        return dalitzplot;
    }
#endif
};

} // namespace GooFit