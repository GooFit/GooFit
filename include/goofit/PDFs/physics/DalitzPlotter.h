#pragma once

#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/physics/Amp3Body.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/Version.h>

#include <algorithm>
#include <numeric>
#include <random>

#if GOOFIT_ROOT_FOUND
#include <TH2.h>
#include <TH1.h>
#include <TFile.h>
#include <TTree.h>
#include <TRatioPlot.h>
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
    fptype daug_1;
    fptype daug_2;
    fptype daug_3;

  public:
    DalitzPlotter(GooPdf *overallSignal, Amp3Body *signalDalitz)
        : m12(signalDalitz->_m12)
        , m13(signalDalitz->_m13)
        , eventNumber(signalDalitz->_eventNumber)
        , data({m12, m13, eventNumber})
        , mother(signalDalitz->decayInfo.motherMass)
        , daug_1(signalDalitz->decayInfo.daug1Mass)
        , daug_2(signalDalitz->decayInfo.daug2Mass)
        , daug_3(signalDalitz->decayInfo.daug3Mass) {
        eventNumber.setValue(0);

        for(size_t i = 0; i < m12.getNumBins(); ++i) {
            m12.setValue(m12.getLowerLimit() + m12.getBinSize() * (i + 0.5));
            for(size_t j = 0; j < m13.getNumBins(); ++j) {
                m13.setValue(m13.getLowerLimit() + m13.getBinSize() * (j + 0.5));
                  //if(m12.getValue()>m13.getValue()) continue;
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

        auto old_data = overallSignal->getData();
        overallSignal->setData(&data);
        signalDalitz->setDataSize(data.getNumEvents());
        pdfValues = overallSignal->getCompProbsAtDataPoints();
        overallSignal->setData(old_data);
    }

    auto getNumEvents() const -> size_t { return data.getNumEvents(); }

    auto getX(size_t event) const -> size_t { return xbins.at(event); }

    auto getY(size_t event) const -> size_t { return ybins.at(event); }

    auto getXval(size_t event) const -> fptype { return data.getValue(m12, event); }

    auto getYval(size_t event) const -> fptype { return data.getValue(m13, event); }

    auto getZval(size_t event) const -> fptype {
        return POW2(mother) + POW2(daug_1) + POW2(daug_2) + POW2(daug_3) - getXval(event) - getYval(event);
    }

    auto getVal(size_t event, size_t num = 0) const -> fptype { return pdfValues.at(num).at(event); }

    auto getDataSet() -> UnbinnedDataSet * { return &data; }

    auto getM12() const -> const Observable & { return m12; }
    auto getM13() const -> const Observable & { return m13; }

#if GOOFIT_ROOT_FOUND
    /// Produce a TH2F over the contained evaluation
    auto make2D(std::string name = "dalitzplot", std::string title = "", int nbins = 100) -> TH2F * {
        auto *dalitzplot = new TH2F(name.c_str(),
                                    title.c_str(),
                                    nbins,
                                    m12.getLowerLimit(),
                                    m12.getUpperLimit(),
                                    nbins,
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

    /// Fill a dataset with MC events
    void fillDataSetMC(UnbinnedDataSet &dataset,
                       size_t nTotal,
                       unsigned int seed     = 0,
                       bool poison           = false,
                       bool save_tree        = false,
                       std::string tree_name = "toyMC.root") {
        TFile *output = nullptr;
        TTree *tree   = nullptr;
        fptype _s12, _s13, _s23, _shigh, _slow;

        if(save_tree) {
            output = new TFile(tree_name.c_str(), "recreate");
            tree   = new TTree("toyTree", "");
            tree->Branch("s12", &_s12);
            tree->Branch("s13", &_s13);
            tree->Branch("s23", &_s23);
            tree->Branch("slow", &_slow);
            tree->Branch("shigh", &_shigh);
        }
        // Setup random numbers
        std::random_device rd;
        std::mt19937 gen;

        if(seed == 0)
            gen.seed(seed);
        else
            gen.seed(rd());

        // Uniform distribution
        std::uniform_real_distribution<> unihalf(-.5, .5);
        std::uniform_real_distribution<> uniwhole(0.0, 1.0);

        // CumSum in other languages
        std::vector<double> integral(pdfValues[0].size());
        std::partial_sum(pdfValues[0].begin(), pdfValues[0].end(), integral.begin());

        // Make this a 0-1 fraction by dividing by the end value
        std::for_each(integral.begin(), integral.end(), [&integral](double &val) { val /= integral.back(); });

        if(poison) {
            std::poisson_distribution<> d(nTotal);
            nTotal = d(gen);
        }

        for(size_t i = 0; i < nTotal; i++) {
            double r = uniwhole(gen);

            // Binary search for integral[cell-1] < r < integral[cell]
            size_t j = std::lower_bound(integral.begin(), integral.end(), r) - integral.begin();

            // Fill in the grid square randomly
            double currm12 = data.getValue(m12, j) + m12.getBinSize() * unihalf(gen);
            double currm13 = data.getValue(m13, j) + m13.getBinSize() * unihalf(gen);

            m12.setValue(currm12);
            m13.setValue(currm13);
            eventNumber.setValue(i);

            if(save_tree) {
                _s12 = currm12;
                _s13 = currm13;
                _s23 = pow(mother, 2) + pow(daug_1, 2) + pow(daug_2, 2) + pow(daug_3, 2) - currm12 - currm13;
                if(currm12 < currm13) {
                    _slow  = currm12;
                    _shigh = currm13;
                } else {
                    _slow  = currm13;
                    _shigh = currm12;
                }
                tree->Fill();
            }

            dataset.addEvent();
        }
        if(save_tree) {
            tree->Write(0, TObject::kOverwrite);
            output->Write(0, TObject::kOverwrite);
            output->Close();
        }
    }

    void Plot(std::string name, UnbinnedDataSet *data, int nbins = 100) {
        auto *Data_DP = new TH2F("Data_DP",
                                 "",
                                 nbins,
                                 m12.getLowerLimit(),
                                 m12.getUpperLimit(),
                                 nbins,
                                 m13.getLowerLimit(),
                                 m13.getUpperLimit());

        auto *Fitted_DP = new TH2F("Fitted_DP",
                                   "",
                                   nbins,
                                   m12.getLowerLimit(),
                                   m12.getUpperLimit(),
                                   nbins,
                                   m13.getLowerLimit(),
                                   m13.getUpperLimit());

        auto s23_min = pow(daug_2+daug_3,2);
        auto s23_max = pow(mother-daug_1,2);

        auto s23_data = new TH1F("s23_data", "", nbins, s23_min, s23_max);
        auto s23_pdf  = new TH1F("s23_pdf", "", nbins, s23_min, s23_max);

        for(int i = 0; i < data->getNumEvents(); i++) {
            data->loadEvent(i);
            Data_DP->Fill(m12.getValue(), m13.getValue());
            s23_data->Fill(pow(mother, 2) + pow(daug_1, 2) + pow(daug_2, 2) + pow(daug_3, 2) - m12.getValue()
                           - m13.getValue());
        }

        for(int i = 0; i < getNumEvents(); i++) {
            s23_pdf->Fill(getZval(i), getVal(i));
            Fitted_DP->Fill(getXval(i), getYval(i), getVal(i));
        }

        auto s12_data = (TH1D *)Data_DP->ProjectionX("s12_data");
        auto s13_data = (TH1D *)Data_DP->ProjectionY("s13_data");
        s12_data->Sumw2();
        s12_data->GetXaxis()->SetTitle("s12");
        s12_data->GetYaxis()->SetTitle("Candidates");
        s13_data->Sumw2();
        s13_data->GetXaxis()->SetTitle("s13");
        s13_data->GetYaxis()->SetTitle("Candidates");

        auto s12_pdf = (TH1D *)Fitted_DP->ProjectionX("s12_pdf");
        s12_pdf->Sumw2();
        auto s13_pdf = (TH1D *)Fitted_DP->ProjectionY("s13_pdf");
        s13_pdf->Sumw2();

        s12_pdf->SetLineColor(kRed);
        s12_pdf->SetLineWidth(2);
        s13_pdf->SetLineColor(kRed);
        s13_pdf->SetLineWidth(2);
        s23_pdf->SetLineColor(kRed);
        s23_pdf->SetLineWidth(2);

        TFile *output = new TFile(name.c_str(), "recreate");
        s12_data->Write(0, TObject::kOverwrite);
        s13_data->Write(0, TObject::kOverwrite);
        s23_data->Write(0, TObject::kOverwrite);
        s12_pdf->Write(0, TObject::kOverwrite);
        s13_pdf->Write(0, TObject::kOverwrite);
        s23_pdf->Write(0, TObject::kOverwrite);

        TCanvas ratio_s12("ratio_s12", "", 1000, 720);
        s12_pdf->Scale(s12_data->Integral() / s12_pdf->Integral());
        auto rt_s12 = new TRatioPlot(s12_data, s12_pdf);
        rt_s12->Draw("divsym");
        ratio_s12.Write(0, TObject::kOverwrite);
        TCanvas ratio_s13("ratio_s13", "", 1000, 720);
        s13_pdf->Scale(s13_data->Integral() / s13_pdf->Integral());
        auto rt_s13 = new TRatioPlot(s13_data, s13_pdf);
        rt_s13->Draw("divsym");
        ratio_s13.Write(0, TObject::kOverwrite);
        TCanvas ratio_s23("ratio_s23", "", 1000, 720);
        s23_pdf->Scale(s23_data->Integral() / s23_pdf->Integral());
        auto rt_s23 = new TRatioPlot(s23_data, s23_pdf);
        rt_s23->Draw("divsym");
        ratio_s23.Write(0, TObject::kOverwrite);
        output->Close();
    }

#else
    /// Fill a dataset with MC events
    void fillDataSetMC(UnbinnedDataSet &dataset, size_t nTotal, unsigned int seed = 0, bool poison = false) {
        // Setup random numbers
        std::random_device rd;
        std::mt19937 gen;

        if(seed == 0)
            gen.seed(seed);
        else
            gen.seed(rd());

        // Uniform distribution
        std::uniform_real_distribution<> unihalf(-.5, .5);
        std::uniform_real_distribution<> uniwhole(0.0, 1.0);

        // CumSum in other languages
        std::vector<double> integral(pdfValues[0].size());
        std::partial_sum(pdfValues[0].begin(), pdfValues[0].end(), integral.begin());

        // Make this a 0-1 fraction by dividing by the end value
        std::for_each(integral.begin(), integral.end(), [&integral](double &val) { val /= integral.back(); });

        if(poison) {
            std::poisson_distribution<> d(nTotal);
            nTotal = d(gen);
        }

        for(size_t i = 0; i < nTotal; i++) {
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

#endif
};

} // namespace GooFit
