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
#include <TH1.h>
#include <TFile.h>
#include <TTree.h>
#include <TRatioPlot.h>
#include <TCanvas.h>
#endif

namespace GooFit {

/// This class makes it easy to make plots over 3 body Dalitz PDFs. You can use ROOT style value access or bin numbers.
class SqDalitzPlotter {
    std::vector<size_t> xbins;
    std::vector<size_t> ybins;
    std::vector<std::vector<fptype>> pdfValues;
    std::vector<fptype> jacobianValues;
    Observable mprime;
    Observable thetaprime;
    EventNumber eventNumber;
    UnbinnedDataSet data;
    fptype mother_mass;
    fptype d1_mass;
    fptype d2_mass;
    fptype d3_mass;
    bool SymDP;
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
        , d3_mass(signalDalitz->decayInfo.daug3Mass)
        , SymDP(signalDalitz->decayInfo.SymDp) {
        eventNumber.setValue(0);

        fptype m23_min = d2_mass+d3_mass;   
        fptype m13_min = d1_mass+d3_mass;  
        fptype m23_max = mother_mass-d1_mass;   
        fptype m13_max = mother_mass-d2_mass; 
        fptype m23 = 0.0;
        fptype m13 = 0.0;
        fptype m12 = 0.0;
        fptype NPTs = 1000000;

        std::random_device rd;  
        std::mt19937 gen(rd()); 
        std::uniform_real_distribution<> random(0.0,1.0);
     

        while(data.getNumEvents()<NPTs){
            fptype _mprime = random(gen);
            fptype _thetaprime = random(gen);
           
            // if(_thetaprime>0.5)
            //     _thetaprime = 1.0-_thetaprime;

            if(!inSqDalitz(_mprime,_thetaprime,SymDP))
                continue;

            if(SymDP && _thetaprime>0.5) _thetaprime = 1.-_thetaprime;

            auto jac = calc_SqDp_Jacobian(_mprime, _thetaprime, mother_mass, d1_mass, d2_mass, d3_mass);
            jacobianValues.push_back(jac);
            
            mprime.setValue(_mprime);
            thetaprime.setValue(_thetaprime);
            data.addEvent();
            eventNumber.setValue(eventNumber.getValue() + 1);

        }

        auto old = overallSignal->getData();
        overallSignal->setData(&data);
        signalDalitz->setDataSize(data.getNumEvents());
        signalDalitz->normalise();
        pdfValues = overallSignal->getCompProbsAtDataPoints();
        std::transform(pdfValues[0].begin(), pdfValues[0].end(), jacobianValues.begin(), pdfValues[0].begin(), std::multiplies<fptype>());
      
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

    void Plot(std::string name, UnbinnedDataSet *data, int nbins = 100) {
        auto *Data_SqDP = new TH2D("Data_SqDP",
                                 "",
                                 nbins,
                                 mprime.getLowerLimit(),
                                 mprime.getUpperLimit(),
                                 nbins,
                                 thetaprime.getLowerLimit(),
                                 thetaprime.getUpperLimit());

        auto *Fitted_SqDP = new TH2D("Fitted_SqDP",
                                   "",
                                   nbins,
                                   mprime.getLowerLimit(),
                                   mprime.getUpperLimit(),
                                   nbins,
                                   thetaprime.getLowerLimit(),
                                   thetaprime.getUpperLimit());
        
        auto *Chi2_SqDP = new TH2D("Chi2_SqDP",
                                   "",
                                   nbins,
                                   mprime.getLowerLimit(),
                                   mprime.getUpperLimit(),
                                   nbins,
                                   thetaprime.getLowerLimit(),
                                   thetaprime.getUpperLimit());

        auto s12_min = pow(d1_mass+d2_mass,2);
        auto s12_max = pow(mother_mass-d3_mass,2);
        auto s12_data = new TH1D("s12_data", "", nbins, 0., s12_max);
        auto s12_pdf  = new TH1D("s12_pdf", "", nbins, 0., s12_max);

        auto s13_min = 0.5;//pow(d1_mass+d3_mass,2);
        auto s13_max = 14.;//pow(mother_mass-d2_mass,2);
        auto s13_data = new TH1D("s13_data", "", nbins, s13_min, s13_max);
        auto s13_pdf  = new TH1D("s13_pdf", "", nbins, s13_min, s13_max);

        auto s23_min = pow(d3_mass+d2_mass,2);
        auto s23_max = pow(mother_mass-d1_mass,2);
        auto s23_data = new TH1D("s23_data", "", nbins, s23_min, s23_max);
        auto s23_pdf  = new TH1D("s23_pdf", "", nbins, s23_min, s23_max);

        auto *Data_DP = new TH2D("Data_DP",
                                 "",
                                 nbins,
                                 s13_min,
                                 s13_max,
                                 nbins,
                                 s23_min,
                                 s23_max);

        auto *Fitted_DP = new TH2D("Fitted_DP",
                                 "",
                                 nbins,
                                 s13_min,
                                 s13_max,
                                 nbins,
                                 s23_min,
                                 s23_max);
        

        for(int i = 0; i < data->getNumEvents(); i++) {
            data->loadEvent(i);
            auto jacobian = 1.;//calc_SqDp_Jacobian(mprime.getValue(),  thetaprime.getValue(),mother_mass,d1_mass,d2_mass,d3_mass);
            Data_SqDP->Fill(mprime.getValue(), thetaprime.getValue());
            auto _m12 = calc_m12(mprime.getValue(),mother_mass,d1_mass,d2_mass,d3_mass);
            auto _m13 = calc_m13(_m12,cos(thetaprime.getValue()*M_PI),mother_mass,d1_mass,d2_mass,d3_mass);
            auto _s12 = _m12*_m12; s12_data->Fill(_s12);
            auto _s13 = _m13*_m13; s13_data->Fill(_s13);
            auto _s23 = mother_mass*mother_mass + d1_mass*d1_mass + d2_mass*d2_mass + d3_mass*d3_mass - _s12 - _s13; s23_data->Fill(_s23);
            Data_DP->Fill(_s13,_s23);
        }

        for(int i = 0; i < getNumEvents(); i++) {
            auto jacobian = 1.;//calc_SqDp_Jacobian(getXval(i), getYval(i),mother_mass,d1_mass,d2_mass,d3_mass);
            Fitted_SqDP->Fill(getXval(i), getYval(i), getVal(i)*jacobian);
            auto _m12 = calc_m12(getXval(i),mother_mass,d1_mass,d2_mass,d3_mass);
            auto _m13 = calc_m13(_m12,cos(getYval(i)*M_PI),mother_mass,d1_mass,d2_mass,d3_mass);
            auto _s12 = _m12*_m12; s12_pdf->Fill(_s12, getVal(i)*jacobian);
            auto _s13 = _m13*_m13; s13_pdf->Fill(_s13, getVal(i)*jacobian);
            auto _s23 = mother_mass*mother_mass + d1_mass*d1_mass + d2_mass*d2_mass + d3_mass*d3_mass - _s12 - _s13; s23_pdf->Fill(_s23, getVal(i)*jacobian);
            Fitted_DP->Fill(_s13,_s23, getVal(i)*jacobian);
        }

        double scale = Data_SqDP->GetEntries()/Fitted_SqDP->GetEntries();
        auto DP_toy = (TH2D*)Fitted_SqDP->Clone();
        Fitted_SqDP->Scale(scale);
        for(int i=0; i<nbins; i++){
                double dt = Data_SqDP->GetBinContent(i);
                double toy = DP_toy->GetBinContent(i);
                double toyerr = DP_toy->GetBinError(i);
                double pdf = Fitted_SqDP->GetBinContent(i);
                double errsq = pow((pdf/toy)*toyerr,2) + dt;
                double res = (pdf-dt)/sqrt(errsq);
                if(std::isnan(res))
                    continue;
                else
                    Chi2_SqDP->SetBinContent(i,res);
        }

        auto mprime_data = (TH1D *)Data_SqDP->ProjectionX("mprime_data");
        auto thetaprime_data = (TH1D *)Data_SqDP->ProjectionY("thetaprime_data");
        mprime_data->Sumw2();
        mprime_data->GetXaxis()->SetTitle("mprime_data");
        mprime_data->GetYaxis()->SetTitle("Candidates");
        thetaprime_data->GetXaxis()->SetTitle("thetaprime_data");
        thetaprime_data->GetYaxis()->SetTitle("Candidates");

        auto mprime_pdf = (TH1D *)Fitted_SqDP->ProjectionX("mprime_pdf");
        auto thetaprime_pdf = (TH1D *)Fitted_SqDP->ProjectionY("thetaprime_pdf");

        mprime_pdf->SetLineColor(kRed);
        mprime_pdf->SetLineWidth(2);
        mprime_pdf->SetLineColor(kRed);
        mprime_pdf->SetLineWidth(2);
        thetaprime_pdf->SetLineColor(kRed);
        thetaprime_pdf->SetLineWidth(2);
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
        mprime_data->Write(0, TObject::kOverwrite);
        thetaprime_data->Write(0, TObject::kOverwrite);
      
        s12_pdf->Write(0, TObject::kOverwrite);
        s13_pdf->Write(0, TObject::kOverwrite);
        s23_pdf->Write(0, TObject::kOverwrite);
        mprime_pdf->Write(0, TObject::kOverwrite);
        thetaprime_pdf->Write(0, TObject::kOverwrite);

        Data_DP->Write(0, TObject::kOverwrite);
        Fitted_DP->Write(0, TObject::kOverwrite);
        Data_SqDP->Write(0, TObject::kOverwrite);
        Fitted_SqDP->Write(0, TObject::kOverwrite);
        Chi2_SqDP->Write(0, TObject::kOverwrite);


        TCanvas ratio_s12("ratio_s12", "", 1000, 720);
        s12_pdf->Scale(s12_data->Integral() / s12_pdf->Integral());
        auto rt_s12 = new TRatioPlot(s12_data, s12_pdf);
        rt_s12->Draw("divsym");
        rt_s12->GetLowYaxis()->SetNdivisions(505);
        ratio_s12.Write(0, TObject::kOverwrite);

        TCanvas ratio_s13("ratio_s13", "", 1000, 720);
        s13_pdf->Scale(s13_data->Integral() / s13_pdf->Integral());
        auto rt_s13 = new TRatioPlot(s13_data, s13_pdf);
        rt_s13->Draw("divsym");
        rt_s13->GetLowYaxis()->SetNdivisions(505);
        ratio_s13.Write(0, TObject::kOverwrite);

        TCanvas ratio_s23("ratio_s23", "", 1000, 720);
        s23_pdf->Scale(s23_data->Integral() / s23_pdf->Integral());
        auto rt_s23 = new TRatioPlot(s23_data, s23_pdf);
        rt_s23->Draw("divsym");
        rt_s23->GetLowYaxis()->SetNdivisions(505);
        ratio_s23.Write(0, TObject::kOverwrite);

        TCanvas ratio_mprime("ratio_mprime", "", 1000, 720);
        mprime_pdf->Scale(mprime_data->Integral() / mprime_pdf->Integral());
        auto rt_mprime = new TRatioPlot(mprime_data, mprime_pdf);
        rt_mprime->Draw("divsym");
         rt_mprime->GetLowYaxis()->SetNdivisions(505);
        ratio_mprime.Write(0, TObject::kOverwrite);

        TCanvas ratio_thetaprime("ratio_thetaprime", "", 1000, 720);
        thetaprime_pdf->Scale(thetaprime_data->Integral() / thetaprime_pdf->Integral());
        auto rt_thetaprime = new TRatioPlot(thetaprime_data, thetaprime_pdf);
        rt_thetaprime->Draw("divsym");
        rt_thetaprime->GetLowYaxis()->SetNdivisions(505);
        ratio_thetaprime.Write(0, TObject::kOverwrite);

        output->Close();
    }
#endif
};

} // namespace GooFit
