// ROOT stuff
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TRandom.h>
#include <TStyle.h>

// System stuff
#include <fstream>

// GooFit stuff
#include <goofit/Application.h>
#include <goofit/BinnedDataSet.h>
#include <goofit/FitControl.h>
#include <goofit/FitManager.h>
#include <goofit/Log.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>

#include <goofit/PDFs/basic/ArgusPdf.h>
#include <goofit/PDFs/basic/GaussianPdf.h>
#include <goofit/PDFs/basic/KinLimitBWPdf.h>
#include <goofit/PDFs/basic/ScaledGaussianPdf.h>
#include <goofit/PDFs/combine/AddPdf.h>
#include <goofit/PDFs/combine/ConvolutionPdf.h>

#include <fmt/format.h>

using namespace fmt::literals;
using namespace GooFit;

TH1D *getMCData(DataSet *data, Observable var, std::string filename, size_t reduce = 1) {
    TH1D *mchist = new TH1D{"mc_hist", "", 300, 0.1365, 0.1665};
    std::ifstream mcreader{filename};

    size_t val_read_in = 0;
    while(mcreader >> var) {
        if(!var)
            continue;
        if(val_read_in % reduce == 0) {
            data->addEvent();
            mchist->Fill(var.getValue());
        }
        val_read_in++;
    }

    mchist->SetStats(false);
    mchist->SetMarkerStyle(8);
    mchist->SetMarkerSize(0.6);

    GOOFIT_INFO("MC events: {}", data->getNumEvents());
    if(data->getNumEvents() == 0)
        throw GooFit::GeneralError("No MC events read in!");
    return mchist;
}

TH1D *getData(DataSet *data, Observable var, std::string filename, size_t reduce = 1) {
    TH1D *data_hist = new TH1D("data_hist", "", 300, 0.1365, 0.1665);
    std::ifstream datareader{filename};

    size_t val_read_in = 0;
    while(datareader >> var) {
        if(!var)
            continue;
        if(val_read_in % reduce == 0) {
            data->addEvent();
            data_hist->Fill(var.getValue());
        }
        val_read_in++;
    }

    data_hist->SetStats(false);
    data_hist->SetMarkerStyle(8);
    data_hist->SetMarkerSize(0.6);

    GOOFIT_INFO("Data events: {}", data->getNumEvents());
    if(data->getNumEvents() == 0)
        throw GooFit::GeneralError("No Data events read in!");

    return data_hist;
}

int main(int argc, char **argv) {
    GooFit::Application app{
        R"raw(This example performs a staged fit measuring the mass difference between the `D*(2010)+`
and `D0` using D*+ -> D0 pi+ events recorded by the BaBar detector (approximately 477
inverse femtobarn).
)raw",
        argc,
        argv};

    app.set_footer(R"raw(Dataset descriptions:
0-simple   Early testing sample for GooFit before nominal dataset was released.
           MC resolution sample and data for channel D*+ -> D0 pi+; D0 -> K- pi+
           Samples are composed of events that pass the majority of selection criteria, but
           fail at least one of the stricter tracking cuts. The resulting resolution is worse
           than in the events of the nominal samples used in the official analysis/publication
           marked below as data set options "1" and "2".
1-kpi      Nominal MC resolution sample and data for channel D*+ -> D0 pi+; D0 -> K- pi+
2-k3pi     Nominal MC resolution sample and data for channel D*+ -> D0 pi+; D0 -> K- pi+ pi- pi+)raw");

    int mode = 0, data = 0;
    bool plot;
    size_t reduce = 1;
    app.add_set("-m,--mode,mode", mode, {0, 1, 2}, "Program mode: 0-unbinned, 1-binned, 2-binned chisq", true);
    app.add_set("-d,--data,data", data, {0, 1, 2}, "Dataset: 0-simple, 1-kpi, 2-k3pi", true);
    app.add_option("--reduce", reduce, "Load every X line of data", true);
    app.add_flag("-p,--plot", plot, "Make and save plots of results");

    GOOFIT_PARSE(app);

    // Style
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

    TCanvas foo;
    foo.SetLogy(true);

    // Get the name of the files to use
    std::string mcfile, datafile;
    if(data == 0) {
        mcfile   = app.get_filename("dataFiles/dstwidth_kpi_resMC.dat", "examples/zachFit");
        datafile = app.get_filename("dataFiles/dstwidth_kpi_data.dat", "examples/zachFit");
    } else if(data == 1) {
        mcfile   = app.get_filename("dataFiles/DstarWidth_D0ToKpi_deltaM_MC.dat", "examples/zachFit");
        datafile = app.get_filename("dataFiles/DstarWidth_D0ToKpi_deltaM_Data.dat", "examples/zachFit");
    } else {
        mcfile   = app.get_filename("dataFiles/DstarWidth_D0ToK3pi_deltaM_MC.dat", "examples/zachFit");
        datafile = app.get_filename("dataFiles/DstarWidth_D0ToK3pi_deltaM_Data.dat", "examples/zachFit");
    }

    Observable dm{"dm", 0.1395, 0.1665};
    dm.setNumBins(2700);

    // This would be clearer with std::optional from C++17
    std::unique_ptr<DataSet> mc_dataset, data_dataset;

    if(mode == 0) {
        mc_dataset.reset(new UnbinnedDataSet{dm});
        data_dataset.reset(new UnbinnedDataSet{dm});
    } else {
        mc_dataset.reset(new BinnedDataSet{dm});
        data_dataset.reset(new BinnedDataSet{dm});
    }

    TH1D *mc_hist = getMCData(mc_dataset.get(), dm, mcfile, reduce);

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

    GaussianPdf gauss1("gauss1", dm, mean1, sigma1);
    GaussianPdf gauss2("gauss2", dm, mean2, sigma2);
    GaussianPdf gauss3("gauss3", dm, mean3, sigma3);
    ArgusPdf argus("argus", dm, pimass, aslope, false, apower);

    AddPdf resolution{"resolution", {gfrac1, gfrac2, afrac}, {&gauss1, &gauss2, &argus, &gauss3}};

    resolution.setData(mc_dataset.get());

    FitManager mcpdf{&resolution};

    GOOFIT_INFO("Done with collecting MC, starting minimisation");
    mcpdf.fit();

    if(plot) {
        GOOFIT_INFO("Plotting MC");
        mc_hist->SetLineColor(kBlack);
        mc_hist->Draw("e");

        double step   = mc_hist->GetXaxis()->GetBinWidth(2);
        auto tot_hist = resolution.plotToROOT(dm, mc_dataset->getNumEvents() * step);
        tot_hist->SetLineColor(kGreen);

        tot_hist->Draw("SAME");

        foo.SaveAs("MC_plot.png");
    }

    // Locking the MC variables
    mean1.setFixed(true);
    mean2.setFixed(true);
    mean3.setFixed(true);
    sigma1.setFixed(true);
    sigma2.setFixed(true);
    sigma3.setFixed(true);
    pimass.setFixed(true);
    aslope.setFixed(true);
    gfrac1.setFixed(true);
    gfrac2.setFixed(true);
    afrac.setFixed(true);
    apower.setFixed(true);

    Variable dummyzero("kpi_rd_dummyzero", 0);
    Variable delta("kpi_rd_delta", 0.000002, -0.00005, 0.00005);
    Variable epsilon("kpi_rd_epsilon", 0.05, -0.1, 0.2);

    ScaledGaussianPdf resolution1("resolution1", dm, dummyzero, sigma1, delta, epsilon);
    ScaledGaussianPdf resolution2("resolution2", dm, dummyzero, sigma2, delta, epsilon);
    ScaledGaussianPdf resolution3("resolution3", dm, dummyzero, sigma3, delta, epsilon);

    Variable width_bw("kpi_rd_width_bw", 0.0001, 0.00001, 0.0005);
    KinLimitBWPdf rbw1("rbw1", dm, mean1, width_bw);
    KinLimitBWPdf rbw2("rbw2", dm, mean2, width_bw);
    KinLimitBWPdf rbw3("rbw3", dm, mean3, width_bw);

    ConvolutionPdf signal1{"signal1", dm, &rbw1, &resolution1};
    ConvolutionPdf signal2{"signal2", dm, &rbw2, &resolution2};
    ConvolutionPdf signal3{"signal3", dm, &rbw3, &resolution3};

    signal1.setIntegrationConstants(0.1395, 0.1665, 0.0000027);
    signal2.setIntegrationConstants(0.1395, 0.1665, 0.0000027);
    signal3.setIntegrationConstants(0.1395, 0.1665, 0.0000027);

    AddPdf signal{"signal", {gfrac1, gfrac2, afrac}, {&signal1, &signal2, &argus, &signal3}};

    Variable slope("kpi_rd_slope", -1.0, 0.1, -35.0, 25.0);
    ArgusPdf bkg("bkg", dm, pimass, slope, false);

    Variable bkg_frac("kpi_rd_bkg_frac", 0.03, 0.0, 0.3);

    TH1D *data_hist = getData(data_dataset.get(), dm, datafile, reduce);

    AddPdf total("total", {bkg_frac}, {&bkg, &signal});

    total.setData(data_dataset.get());

    std::shared_ptr<BinnedChisqFit> chi_control;
    if(2 == mode) {
        chi_control.reset(new BinnedChisqFit);
        total.setFitControl(chi_control);
    }

    FitManager datapdf{&total};

    GOOFIT_INFO("Starting fit");

    datapdf.fit();

    if(plot) {
        GOOFIT_INFO("Plotting results");

        data_hist->SetLineColor(kBlack);
        data_hist->Draw("e");

        double scale = data_hist->GetXaxis()->GetBinWidth(2) * data_dataset->getNumEvents();

        auto sig_hist = signal.plotToROOT(dm, (1 - bkg_frac.getValue()) * scale);
        sig_hist->SetLineColor(kBlue);
        auto back_hist = bkg.plotToROOT(dm, bkg_frac.getValue() * scale);
        back_hist->SetLineColor(kRed);
        auto tot_hist = total.plotToROOT(dm, scale);
        tot_hist->SetLineColor(kGreen);

        tot_hist->Draw("SAME");
        sig_hist->Draw("SAME");
        back_hist->Draw("SAME");

        foo.SaveAs("ResultFit.png");
    }

    return datapdf;
}
