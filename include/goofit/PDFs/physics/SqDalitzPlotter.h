#pragma once

#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/physics/Amp3BodySqDP.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/Version.h>

#include <algorithm>
#include <numeric>
#include <random>
#include <iostream>

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
    std::vector<fptype> m12_vec;
    std::vector<fptype> m13_vec;
    

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

        fptype m23_min = d2_mass+d3_mass;   
        fptype m13_min = d1_mass+d3_mass;  
        fptype m23_max = mother_mass-d1_mass;   
        fptype m13_max = mother_mass-d2_mass; 
        fptype m23 = 0.0;
        fptype m13 = 0.0;
        fptype m12 = 0.0;
        fptype NPTs = 10000000;

        std::random_device rd;  
        std::mt19937 gen(rd()); 
        std::uniform_real_distribution<> m23_gen(m23_min*m23_min, m23_max*m23_max);
        std::uniform_real_distribution<> m13_gen(m13_min*m13_min, m13_max*m13_max);

        while(data.getNumEvents()<NPTs){
            m23 = m23_gen(gen);
            m13 = m13_gen(gen);
            
            if(inDalitz2(m13,
                m23,
                mother_mass,
                d1_mass,
                d2_mass,
                d3_mass)) {
                    
                    fptype m12 = sqrt(mother_mass*mother_mass + d1_mass*d1_mass + d2_mass*d2_mass + d3_mass*d3_mass - m23- m13);
                    
                    fptype mp = calc_mprime(m12, mother_mass, d1_mass, d2_mass, d3_mass);
                    fptype th = calc_thetaprime(m12, sqrt(m13), mother_mass, d1_mass, d2_mass, d3_mass);

                    if(d2_mass==d3_mass){
                         if(th>0.5)
                            th = 1.0-th;
                    }
                  
                    
                    mprime.setValue(mp);
                    thetaprime.setValue(th);
                    data.addEvent();
                    eventNumber.setValue(eventNumber.getValue() + 1);
            }

        }

        // std::random_device rd;  
        // std::mt19937 gen(rd()); 
        // std::uniform_real_distribution<> random(0.0,1.0);
     

        // while(data.getNumEvents()<NPTs){
        //     fptype _mprime = random(gen);
        //     fptype _thetaprime = random(gen);
           
        //     if(_thetaprime>0.5)
        //         _thetaprime = 1.0-_thetaprime;
            
        //     mprime.setValue(_mprime);
        //     thetaprime.setValue(_thetaprime);
        //     data.addEvent();
        //     eventNumber.setValue(eventNumber.getValue() + 1);

        // }

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
        //std::uniform_real_distribution<> unihalf(-.5, +.5);
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
            
            fptype currmprime = data.getValue(mprime, j);
            fptype currthetaprime = data.getValue(thetaprime, j);
           
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
