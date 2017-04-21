// ROOT stuff
#include "TRandom.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TStyle.h"

// System stuff
#include <fstream>
#include <exception>
#include <sys/time.h>
#include <sys/times.h>

// RooFit stuff
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooMinuit.h"
#include "RooNLLVar.h"

// GooFit stuff
#include "goofit/Application.h"
#include "goofit/Variable.h"
#include "goofit/UnbinnedDataSet.h"
#include "goofit/BinnedDataSet.h"
#include "goofit/FitManager.h"
#include "goofit/FitControl.h"

#include "goofit/PDFs/KinLimitBWPdf.h"
#include "goofit/PDFs/ConvolutionPdf.h"
#include "goofit/PDFs/GaussianPdf.h"
#include "goofit/PDFs/ScaledGaussianPdf.h"
#include "goofit/PDFs/ArgusPdf.h"
#include "goofit/PDFs/AddPdf.h"


TCanvas* foo;
timeval startTime, stopTime, totalTime;
clock_t startCPU, stopCPU;
tms startProc, stopProc;
BinnedDataSet* binnedData = 0;
UnbinnedDataSet* data = 0;
int length = 0;

TH1F* data_hist = 0;
Variable* dm;

double pdf_int;

char histName[1000];
int numHists = 0;

TH1F* plotComponent(GooPdf* toPlot, double normFactor) {
//  static char name[1000];
//  static int numHists = 0;
    sprintf(histName, "%s_hist_%i", toPlot->getName().c_str(), numHists++);
    TH1F* ret = new TH1F(histName, "", dm->numbins, dm->lowerlimit, dm->upperlimit);
    std::vector<fptype> binValues;
    toPlot->evaluateAtPoints(dm, binValues);

    pdf_int = 0;
    double step = dm->upperlimit - dm->lowerlimit;
    step /= dm->numbins;

    for(int i = 1; i <= dm->numbins; ++i) {
        //std::cout << name << " " << i << " : " << binValues[i-1] << " " << (dm->lowerlimit + (i-1)*step) << std::endl;
        pdf_int += binValues[i-1];
    }

    for(int i = 1; i <= dm->numbins; ++i)
        ret->SetBinContent(i, binValues[i-1] * normFactor / pdf_int);

    return ret;
}

void getMCData() {
    data = new UnbinnedDataSet(dm);
    std::ifstream mcreader;
    mcreader.open("dataFiles/dstwidth_kpi_resMC.dat");
    TH1F* mchist = new TH1F("mchist", "", 300, 0.1365, 0.1665);

    double currDM = 0;

    while(true) {
        mcreader >> currDM;

        if(mcreader.eof())
            break;

        if(currDM < 0.13957)
            std::cout << "Bad DM\n";

        data->addEvent(currDM);
        mchist->Fill(currDM);
    }

    mchist->SetStats(false);
    mchist->SetMarkerStyle(8);
    mchist->SetMarkerSize(0.6);
    mchist->Draw("p");

    foo->SetLogy(true);
    foo->SaveAs("zach_mchist.png");

    std::cout << "MC: Got " << data->getNumEvents() << " events.\n";
}

void getData() {
    std::ifstream datareader;
    std::string filename = "dataFiles/zach/dstwidth_kpi_data.dat";

    if(!((bool) std::ifstream(filename))) {
        std::cerr << "The file " << filename << " does not exist! Exiting" << std::endl;
        throw std::runtime_error("Missing data file");
    }

    datareader.open(filename);

    binnedData = new BinnedDataSet(dm);
    delete data;
    data = new UnbinnedDataSet(dm);
    double currDM = 0;

    while(true) {
        datareader >> currDM;

        if(datareader.eof())
            break;

        if(currDM > dm->upperlimit)
            continue;

        if(currDM < dm->lowerlimit)
            continue;

        data->addEvent(currDM);
        data_hist->Fill(currDM);

        binnedData->addEvent(currDM);
    }

    std::cout << "Data events: " << data->getNumEvents() << std::endl;
    datareader.close();
}

void CudaMinimise(int fitType) {
    dm = new Variable("dm", 0.1395, 0.1665);
    dm->numbins = 2700;
    //dm->numbins = 540;

    getMCData();
    std::cout << "Done getting MC\n";

    Variable mean1("kpi_mc_mean1", 0.145402, 0.00001, 0.143, 0.148);
    Variable mean2("kpi_mc_mean2", 0.145465, 0.00001, 0.145, 0.1465);
    Variable mean3("kpi_mc_mean3", 0.145404, 0.00001, 0.144, 0.147);

    Variable sigma1("kpi_mc_sigma1", 0.00010, 0.00001, 0.000001, 0.002);
    Variable sigma2("kpi_mc_sigma2", 0.00075, 0.00001, 0.000001, 0.005);
    Variable sigma3("kpi_mc_sigma3", 0.00020, 0.00001, 0.000005, 0.001);

    Variable pimass("kpi_mc_pimass", 0.13957);
    Variable aslope("kpi_mc_aslope", -20.0, 1, -100.0, 10.0);
    Variable apower("kpi_mc_apower", 1.3, 0.1, 0.1, 10.0);
    Variable gfrac1("kpi_mc_gfrac1", 0.65, 0.01, 0.0, 0.9);
    Variable gfrac2("kpi_mc_gfrac2", 0.02, 0.001, 0.0, 0.12);
    Variable afrac("kpi_mc_afrac", 0.005, 0.003, 0.0, 0.1);

    //mean1.fixed = true;
    //mean2.fixed = true;
    //mean3.fixed = true;
    //sigma1.fixed = true;
    //sigma2.fixed = true;
    //sigma3.fixed = true;

    //aslope.fixed = true;
    //apower.fixed = true;
    //gfrac1.fixed = true;
    //gfrac2.fixed = true;
    //afrac.fixed = true;

    GaussianPdf gauss1("gauss1", dm, &mean1, &sigma1);
    GaussianPdf gauss2("gauss2", dm, &mean2, &sigma2);
    GaussianPdf gauss3("gauss3", dm, &mean3, &sigma3);
    ArgusPdf argus("argus", dm, &pimass, &aslope, false, &apower);

    std::vector<Variable*> weights;
    weights.push_back(&gfrac1);
    weights.push_back(&gfrac2);
    weights.push_back(&afrac);

    std::vector<PdfBase*> comps;
    comps.push_back(&gauss1);
    comps.push_back(&gauss2);
    comps.push_back(&argus);
    comps.push_back(&gauss3);

    AddPdf resolution("resolution", weights, comps);
    resolution.setData(data);
    FitManager mcpdf(&resolution);

    std::cout << "Done with data, starting minimisation" << std::endl;
    // Minimize
    //ROOT::Minuit2::FunctionMinimum* min = mcpdf.fit();
    mcpdf.fit();

    mean1.fixed = true;
    mean2.fixed = true;
    mean3.fixed = true;
    sigma1.fixed = true;
    sigma2.fixed = true;
    sigma3.fixed = true;
    pimass.fixed = true;
    aslope.fixed = true;
    gfrac1.fixed = true;
    gfrac2.fixed = true;
    afrac.fixed = true;
    apower.fixed = true;

    Variable dummyzero("kpi_rd_dummyzero", 0);
    Variable delta("kpi_rd_delta", 0.000002, -0.00005, 0.00005);
    Variable epsilon("kpi_rd_epsilon", 0.05, -0.1, 0.2);

    ScaledGaussianPdf resolution1("resolution1", dm, &dummyzero, &sigma1, &delta, &epsilon);
    ScaledGaussianPdf resolution2("resolution2", dm, &dummyzero, &sigma2, &delta, &epsilon);
    ScaledGaussianPdf resolution3("resolution3", dm, &dummyzero, &sigma3, &delta, &epsilon);

    Variable width_bw("kpi_rd_width_bw", 0.0001, 0.00001, 0.0005);
    KinLimitBWPdf rbw1("rbw1", dm, &mean1, &width_bw);
    KinLimitBWPdf rbw2("rbw2", dm, &mean2, &width_bw);
    KinLimitBWPdf rbw3("rbw3", dm, &mean3, &width_bw);

    //#define OTHERS 1
#ifdef OTHERS
    ConvolutionPdf signal1("signal1", dm, &rbw1, &resolution1, 2);
    ConvolutionPdf signal2("signal2", dm, &rbw2, &resolution2, 2);
    ConvolutionPdf signal3("signal3", dm, &rbw3, &resolution3, 2);
    std::vector<ConvolutionPdf*> convList;
    convList.push_back(&signal1);
    convList.push_back(&signal2);
    convList.push_back(&signal3);
    signal1.registerOthers(convList);
    signal2.registerOthers(convList);
    signal3.registerOthers(convList);
#else
    ConvolutionPdf signal1("signal1", dm, &rbw1, &resolution1);
    ConvolutionPdf signal2("signal2", dm, &rbw2, &resolution2);
    ConvolutionPdf signal3("signal3", dm, &rbw3, &resolution3);
#endif

    signal1.setIntegrationConstants(0.1395, 0.1665, 0.0000027);
    signal2.setIntegrationConstants(0.1395, 0.1665, 0.0000027);
    signal3.setIntegrationConstants(0.1395, 0.1665, 0.0000027);



    weights.clear();
    weights.push_back(&gfrac1);
    weights.push_back(&gfrac2);
    weights.push_back(&afrac);

    comps.clear();
    comps.push_back(&signal1);
    comps.push_back(&signal2);
    comps.push_back(&argus);
    comps.push_back(&signal3);
    AddPdf signal("signal", weights, comps);

    Variable slope("kpi_rd_slope", -1.0, 0.1, -35.0, 25.0);
    Variable* bpower = NULL;
    ArgusPdf bkg("bkg", dm, &pimass, &slope, false, bpower);

    weights.clear();
    comps.clear();

    Variable bkg_frac("kpi_rd_bkg_frac", 0.03, 0.0, 0.3);
    weights.push_back(&bkg_frac);
    comps.push_back(&bkg);
    comps.push_back(&signal);

    getData();

    AddPdf total("total", weights, comps);

    if(0 == fitType)
        total.setData(data);
    else {
        total.setData(binnedData);

        if(2 == fitType)
            total.setFitControl(new BinnedChisqFit());
    }

    FitManager datapdf(&total);

    std::cout << "Starting fit\n";
    gettimeofday(&startTime, NULL);
    startCPU = times(&startProc);
    //ROOT::Minuit2::FunctionMinimum* min2 = datapdf.fit();
    datapdf.fit();
    stopCPU = times(&stopProc);
    gettimeofday(&stopTime, NULL);

    //std::cout << "Minimum: " << *min2 << std::endl;
    /*
      double dat_int = 0;
      for (int i = 1; i <= 300; ++i) {
        dat_int += data_hist->GetBinContent(i);
      }

      signal1.setIntegrationConstants(0.1365, 0.1665, 0.00003);
      signal2.setIntegrationConstants(0.1365, 0.1665, 0.00003);
      signal3.setIntegrationConstants(0.1365, 0.1665, 0.00003);
      dm->numbins = 300;
      dm->lowerlimit = 0.1365;
      dm->upperlimit = 0.1665;
      std::cout << bkg_frac.value << std::endl;

      // plotComponent seems broken?
      TH1F* dpdf_hist = plotComponent(&total, dat_int);
      double totalIntegral = pdf_int;
      TH1F* barg_hist = plotComponent(&bkg,   dat_int*bkg_frac.value);
      double bkgIntegral = pdf_int;
      TH1F* sign_hist = plotComponent(&signal, dat_int*(1 - bkg_frac.value));
      double sigIntegral = pdf_int;

      double sig_int = 0;
      double bkg_int = 0;
      double tot_int = 0;
      for (int i = 1; i <= 300; ++i) {
        sig_int += sign_hist->GetBinContent(i);
        bkg_int += barg_hist->GetBinContent(i);
        tot_int += dpdf_hist->GetBinContent(i);

        dpdf_hist->SetBinContent(i, barg_hist->GetBinContent(i) + sign_hist->GetBinContent(i));
      }
    */

    dm->value = 0.1568;
    /*
      std::cout << "PDF: "
    	    << (dat_int/totalIntegral) * total.getValue() << " "
    	    << (1-bkg_frac.value)*(dat_int/sigIntegral)*signal.getValue() << " "
    	    << bkg_frac.value*(dat_int/bkgIntegral)*bkg.getValue() << " | "
    	    << dat_int << " " << sigIntegral << " " << bkgIntegral << " " << totalIntegral << " | "
    	    << sig_int << " " << bkg_int << " " << tot_int  << " | "
    	    << dpdf_hist->GetBinContent(204) << " " << sign_hist->GetBinContent(204) << " " << barg_hist->GetBinContent(204) << " "
    	    << std::endl;
    */

    /*
    data_hist->SetStats(false);
    data_hist->SetMarkerStyle(8);
    data_hist->SetMarkerSize(0.6);
    data_hist->Draw("p");

    dpdf_hist->SetLineColor(kViolet);
    dpdf_hist->SetLineWidth(3);
    dpdf_hist->Draw("lsame");
    //dpdf_hist->Draw("l");

    barg_hist->SetLineColor(kRed);
    barg_hist->SetLineWidth(3);
    barg_hist->SetLineStyle(kDashed);
    barg_hist->Draw("lsame");

    sign_hist->SetLineColor(kBlue);
    sign_hist->SetLineWidth(3);
    sign_hist->SetLineStyle(kDashed);
    sign_hist->Draw("lsame");

    foo->SetLogy(false);
    foo->SaveAs("zach_linear_CUDA_fit.png");
    foo->SetLogy(true);
    foo->SaveAs("zach_CUDA_fit.png");
    */
}

int main(int argc, char** argv) {
    GooFit::Application app("Zach-Fit example", argc, argv);
    
    int mode;
    app.add_set("-m,--mode,mode", mode, {0,1,2},
            "Program mode: 0-unbinned, 1-binned, 2-binned ChiSq")->required();


    try {
        app.run();
    } catch (const GooFit::ParseError &e) {
        return app.exit(e);
    }

    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasColor(10);
    gStyle->SetFrameFillColor(10);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetPadColor(0);
    gStyle->SetTitleColor(1);
    gStyle->SetStatColor(0);
    gStyle->SetFillColor(0);
    gStyle->SetFuncWidth(1);
    gStyle->SetLineWidth(1);
    gStyle->SetLineColor(1);
    gStyle->SetPalette(1, 0);
    foo = new TCanvas();

    data_hist = new TH1F("data_hist", "", 300, 0.1365, 0.1665);

    try {
        CudaMinimise(mode);
    } catch(const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        return 6;
    }

    //RooFitMinimise(atoi(argv[3]));

    data_hist->SetStats(false);
    data_hist->SetMarkerStyle(8);
    data_hist->SetMarkerSize(0.6);
    data_hist->Draw("p");
    foo->Draw();
    foo->SaveAs("zachDraw.png");


    // Print total minimization time
    double myCPU = stopCPU - startCPU;
    double totalCPU = myCPU;

    timersub(&stopTime, &startTime, &totalTime);
    std::cout << "Wallclock time  : " << totalTime.tv_sec + totalTime.tv_usec/1000000.0 << " seconds." << std::endl;
    std::cout << "CPU time: " << (myCPU / CLOCKS_PER_SEC) << std::endl;
    std::cout << "Total CPU time: " << (totalCPU / CLOCKS_PER_SEC) << std::endl;
    myCPU = stopProc.tms_utime - startProc.tms_utime;
    std::cout << "Processor time: " << (myCPU / CLOCKS_PER_SEC) << std::endl;

    delete binnedData;

    return 0;
}
