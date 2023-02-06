#pragma once

#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/physics/Amp3BodySqDP.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/Version.h>

#include <algorithm>
#include <numeric>
#include <random>

#if GOOFIT_ROOT_FOUND
#include <TH2.h>
#endif

namespace GooFit {

/// This class makes it easy to make plots over 3 body Dalitz PDFs. You can use ROOT style value access or bin numbers.
class SqDalitzPlotter {
    std::vector<size_t> xbins;
    std::vector<size_t> ybins;
    std::vector<std::vector<fptype>> pdfValues;
    Observable mprime;
    Observable thetaprime;
    EventNumber eventNumber;
    UnbinnedDataSet data;
    fptype mother_mass;
    fptype d1_mass;
    fptype d2_mass;
    fptype d3_mass;

  public:
    SqDalitzPlotter(GooPdf *overallSignal, Amp3BodySqDP *signalDalitz)
        : mprime(signalDalitz->_mprime)
        , thetaprime(signalDalitz->_thetaprime)
        , eventNumber(signalDalitz->_eventNumber)
        , data({mprime, thetaprime, eventNumber})
        , mother_mass(signalDalitz->decayInfo.motherMass)
        , d1_mass(signalDalitz->decayInfo.daug1Mass)
        , d2_mass(signalDalitz->decayInfo.daug2Mass)
        , d3_mass(signalDalitz->decayInfo.daug3Mass) {
        eventNumber.setValue(0);

        std::default_random_engine generator;
        std::uniform_real_distribution<fptype> distribution(0.0,1.0);
        

        for(size_t i = 0; i < mprime.getNumBins(); ++i) {
            mprime.setValue(mprime.getLowerLimit() + mprime.getBinSize() * (i + 0.5));
            for(size_t j = 0; j < thetaprime.getNumBins(); ++j) {
                thetaprime.setValue(thetaprime.getLowerLimit() + thetaprime.getBinSize() * (j + 0.5));
                if(inSqDalitz(mprime.getValue(),
                            thetaprime.getValue())) {
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

        for(size_t i=0; i<data.getNumEvents();i++){
            data.loadEvent(i);
            fptype jacobian = calc_SqDp_Jacobian(mprime.getValue(), thetaprime.getValue(), mother_mass, d1_mass, d2_mass, d3_mass);
            //pdfValues.at(0).at(i) *= jacobian;
        }


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
        std::vector<fptype> integral(pdfValues[0].size());
        std::partial_sum(pdfValues[0].begin(), pdfValues[0].end(), integral.begin());

        // Make this a 0-1 fraction by dividing by the end value
        std::for_each(integral.begin(), integral.end(), [&integral](fptype &val) { val /= integral.back(); });

        for(size_t i = 0; i < nTotal; i++) {
            fptype r = uniwhole(gen);

            // Binary search for integral[cell-1] < r < integral[cell]
            size_t j = std::lower_bound(integral.begin(), integral.end(), r) - integral.begin();

            // Fill in the grid square randomly
            fptype currmprime = data.getValue(mprime, j) + mprime.getBinSize() ;//* unihalf(gen);
            fptype currthetaprime = data.getValue(thetaprime, j) + thetaprime.getBinSize();// * unihalf(gen);
            

            //printf("mprime=%.2f \t thetaprime=%.2f \t prob=%.4f \n", currmprime, currthetaprime, pdfValues[0][j]);

            mprime.setValue(currmprime);
            thetaprime.setValue(currthetaprime);
            eventNumber.setValue(i);
            dataset.addEvent();
        }
    }

    auto getNumEvents() const -> size_t { return data.getNumEvents(); }

    auto getX(size_t event) const -> size_t { return xbins.at(event); }

    auto getY(size_t event) const -> size_t { return ybins.at(event); }

    auto getXval(size_t event) const -> fptype { return data.getValue(mprime, event); }

    auto getYval(size_t event) const -> fptype { return data.getValue(thetaprime, event); }

    auto getZval(size_t event) const -> fptype { return POW2(mother_mass) - POW2(getXval(event)) - POW2(getYval(event)); }

    auto getVal(size_t event, size_t num = 0) const -> fptype { return pdfValues.at(num).at(event); }

    auto getDataSet() -> UnbinnedDataSet * { return &data; }

    auto getM12() const -> const Observable & { return mprime; }
    auto getM13() const -> const Observable & { return thetaprime; }

#if GOOFIT_ROOT_FOUND
    /// Produce a TH2F over the contained evaluation
    auto make2D(std::string name = "dalitzplot", std::string title = "") -> TH2F * {
        auto *dalitzplot = new TH2F(name.c_str(),
                                    title.c_str(),
                                    mprime.getNumBins(),
                                    mprime.getLowerLimit(),
                                    mprime.getUpperLimit(),
                                    thetaprime.getNumBins(),
                                    thetaprime.getLowerLimit(),
                                    thetaprime.getUpperLimit());

        

        for(unsigned int j = 0; j < getNumEvents(); ++j) {
           
            size_t currmprime = getX(j);
            size_t currthetaprime = getY(j);
            
            fptype val     = getVal(j);

            dalitzplot->SetBinContent(1 + currmprime, 1 + currthetaprime, val);
        }

        return dalitzplot;
    }
#endif
};

} // namespace GooFit
