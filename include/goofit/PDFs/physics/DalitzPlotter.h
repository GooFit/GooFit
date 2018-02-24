#pragma once

#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/DalitzPlotPdf.h>
#include <goofit/Version.h>

#include <algorithm>
#include <numeric>
#include <random>

#if GOOFIT_ROOT_FOUND
#include <TH2.h>
#endif

namespace GooFit {

/// This class makes it easy to make plots over 3 body Dalitz PDFs. You can use ROOT style value access or bin numbers.
class DalitzPlotter {
    std::vector<size_t> xbins;
    std::vector<size_t> ybins;
    std::vector<std::vector<fptype>> pdfValues;
    Observable m12;
    Observable m13;
    EventNumber eventNumber;
    UnbinnedDataSet data;
    fptype mother;

  public:
    DalitzPlotter(GooPdf *overallSignal, DalitzPlotPdf *signalDalitz)
        : m12(signalDalitz->_m12)
        , m13(signalDalitz->_m13)
        , eventNumber(signalDalitz->_eventNumber)
        , data({m12, m13, eventNumber})
        , mother(signalDalitz->decayInfo.motherMass) {
        eventNumber.setValue(0);

        for(size_t i = 0; i < m12.getNumBins(); ++i) {
            m12.setValue(m12.getLowerLimit()
                         + (m12.getUpperLimit() - m12.getLowerLimit()) * (i + 0.5) / m12.getNumBins());
            for(size_t j = 0; j < m13.getNumBins(); ++j) {
                m13.setValue(m13.getLowerLimit()
                             + (m13.getUpperLimit() - m13.getLowerLimit()) * (j + 0.5) / m13.getNumBins());
                if(inDalitz(m12.getValue(),
                            m13.getValue(),
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

        overallSignal->setData(&data);
        signalDalitz->setDataSize(data.getNumEvents());

        pdfValues = overallSignal->getCompProbsAtDataPoints();
    }

    /// Fill a dataset with MC events
    void fillDataSetMC(UnbinnedDataSet &dataset, size_t nTotal) {

        // Setup random numbers
        std::random_device rd;
        std::mt19937 gen(rd());

        // Poisson distribution
        std::poisson_distribution<> d(nTotal);
        size_t num_events = d(gen);

        // Uniform distribution
        std::uniform_real_distribution<> unihalf(-.5, .5);
        std::uniform_real_distribution<> uniwhole(0.0, 1.0);

        // CumSum in other languages
        std::vector<double> integral(pdfValues[0].size());
        std::partial_sum(pdfValues[0].begin(), pdfValues[0].end(), integral.begin());

        // Make this a 0-1 fraction by dividing by the end value
        std::for_each(integral.begin(), integral.end(), [&integral](double &val) { val /= integral.back(); });

        for(size_t i = 0; i < num_events; i++) {
            double r = uniwhole(gen);

            // Binary search for integral[cell-1] < r < integral[cell]
            size_t j = std::lower_bound(integral.begin(), integral.end(), r) - integral.begin();

            // Fill in the grid square randomly
            double currm12 = data.getValue(m12, j) + m12.getBinSize() * unihalf(gen);
            double currm13 = data.getValue(m13, j) + m13.getBinSize() * unihalf(gen);

            m12.setValue(currm12);
            m13.setValue(currm13);
            eventNumber.setValue(i);
            dataset.addEvent();
        }
    }

    size_t getNumEvents() const { return data.getNumEvents(); }

    size_t getX(size_t event) const { return xbins.at(event); }

    size_t getY(size_t event) const { return ybins.at(event); }

    fptype getXval(size_t event) const { return data.getValue(m12, event); }

    fptype getYval(size_t event) const { return data.getValue(m13, event); }

    fptype getZval(size_t event) const { return POW2(mother) - POW2(getXval(event)) - POW2(getYval(event)); }

    fptype getVal(size_t event, size_t num = 0) const { return pdfValues.at(num).at(event); }

    UnbinnedDataSet *getDataSet() { return &data; }

    const Observable &getM12() const { return m12; }
    const Observable &getM13() const { return m13; }

#if GOOFIT_ROOT_FOUND
    /// Produce a TH2F over the contained evaluation
    TH2F *make2D(std::string name = "dalitzplot", std::string title = "") {
        TH2F *dalitzplot = new TH2F(name.c_str(),
                                    title.c_str(),
                                    m12.getNumBins(),
                                    m12.getLowerLimit(),
                                    m12.getUpperLimit(),
                                    m13.getNumBins(),
                                    m13.getLowerLimit(),
                                    m13.getUpperLimit());

        for(unsigned int j = 0; j < getNumEvents(); ++j) {
            size_t currm12 = getX(j);
            size_t currm13 = getY(j);
            double val     = getVal(j);

            dalitzplot->SetBinContent(1 + currm12, 1 + currm13, val);
        }

        return dalitzplot;
    }
#endif
};

} // namespace GooFit
