// ROOT stuff
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TChain.h>
#include <TLegend.h>
#include <TLine.h>
#include <TRandom.h>
#include <TRandom3.h>

// System stuff
#include <fstream>
#include <sys/time.h>
#include <sys/times.h>

// GooFit stuff
#include <goofit/Application.h>
#include <goofit/FitManager.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/basic/PolynomialPdf.h>
#include <goofit/PDFs/combine/AddPdf.h>
#include <goofit/PDFs/combine/ProdPdf.h>
#include <goofit/PDFs/physics/Amp3Body.h>
#include <goofit/PDFs/physics/DalitzPlotter.h>
#include <goofit/PDFs/physics/DalitzVetoPdf.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>
#include <goofit/utilities/Style.h>

using namespace std;
using namespace GooFit;

Variable fixedRhoMass("rho_mass", 0.7758, 0.01, 0.7, 0.8);
Variable fixedRhoWidth("rho_width", 0.1503, 0.01, 0.1, 0.2);
char strbuffer[1000]; 

const fptype _mDs       = 1.96828;
const fptype _mDs2      = _mDs * _mDs;
const fptype _mDs2inv   = 1. / _mDs2;
const fptype piPlusMass = 0.13957018;
const fptype piZeroMass = 0.1349766;

// Constants used in more than one PDF component.
Variable motherM("motherM", _mDs);
Variable chargeM("chargeM", piPlusMass);
Variable neutrlM("neutrlM", piZeroMass);
Variable massSum("massSum", _mDs *_mDs + 2 * piPlusMass * piPlusMass + piZeroMass * piZeroMass); // = 3.53481
Variable constantOne("constantOne", 1);
Variable constantZero("constantZero", 0);

vector<fptype> HH_bin_limits;
vector<Variable> pwa_coefs_reals;
vector<Variable> pwa_coefs_imags;

fptype cpuGetM23(fptype massPZ, fptype massPM) {
    return (_mDs2 + piZeroMass * piZeroMass + piPlusMass * piPlusMass + piPlusMass * piPlusMass - massPZ - massPM);
}

ResonancePdf* loadPWAResonance(const string fname = "PWA_COEFFS_30_SAda_reim.txt", bool fixAmp = false){
  std::ifstream reader;
  reader.open(fname.c_str()); 
  assert(reader.good());
  HH_bin_limits.clear();
  pwa_coefs_reals.clear();
  pwa_coefs_imags.clear();

  TRandom3 *gRandom = new TRandom3(time(NULL)); 
  double Rmag = 1., Rphs = 1.;
/*  Rmag = gRandom->Uniform(0.5,1.5);
  Rphs = gRandom->Uniform(0.8,1.2);*/

  double e1,e2,e3,e4;
  int i = 0;
//  while (reader >> e1 >> e2 >> e3 >> e4) {
  while (reader >> e1 >> e2 >> e3 ) {
//      if (i%3 ==0||i==29){
//      HH_bin_limits.push_back(e1*e1);
      HH_bin_limits.push_back(e1);
/*      sprintf(strbuffer, "pwa_coef_%d_real", i);
      Variable* vr = new Variable(strbuffer, e2, .1, -10,10);//e2>0?0:10*e2, e2>0?10*e2:0);
      sprintf(strbuffer, "pwa_coef_%d_imag", i);
      Variable* vi = new Variable(strbuffer, e3, .1, -10,10);//e3>0?0:10*e3, e3>0?10*e3:0);*/
      sprintf(strbuffer, "pwa_coef_%d_mag", i);
      Variable* vr = new Variable(strbuffer, Rmag*e2, .01, 0,1e2);
      sprintf(strbuffer, "pwa_coef_%d_phase", i);
//      Variable* vi = new Variable(strbuffer, Rphs*e3, .01, -3*TMath::Pi(), 3*TMath::Pi());
      Variable* vi = new Variable(strbuffer, Rphs*e3, .01, 0,0);
//      if (i==29) 
//      if (i==29) 
/*          vr->fixed = vi->fixed = true;
      if (i==-9) */
      vr->setFixed(true); vi->setFixed(true);

//      if (!(i>=1&&i<=10)) vr->fixed = vi->fixed = true;
//      if (i==10) vr->fixed = vi->fixed = true;
//      if (!(i>=1&&i<=2)) vr->fixed = vi->fixed = true;
//      if (i==29) vr->fixed = vi->fixed = true;
//      if (!(i>=26&&i<=29)) vr->fixed = vi->fixed = true;
      pwa_coefs_reals.push_back(*vr);
      pwa_coefs_imags.push_back(*vi);
//      }
      i++;
  }
/*  if (SEED>=0){  
  const double angspread = 4*TMath::Pi();
  double angintv = angspread/(i-1);
  for (int ii=0;ii<i;ii++){
      pwa_coefs_reals[ii]->value = 0.1*SEED;  
      pwa_coefs_imags[ii]->value = -0.5*angspread+angintv*ii;
//      pwa_coefs_imags[ii]->lowerlimit = pwa_coefs_imags[ii]->value -1.2*TMath::Pi();
//      pwa_coefs_imags[ii]->upperlimit = pwa_coefs_imags[ii]->value + 1.2*TMath::Pi();
  }}*/
  const fptype scale = 1;
  Variable* swave_amp_real = new Variable("swave_amp_real", -1.86396*scale,   0.001, 0, 0);
  Variable* swave_amp_imag = new Variable("swave_amp_imag", .775892256*scale,   0.001, 0, 0);
  swave_amp_real->setValue(3.4504137);//3.44482e+00; 
  swave_amp_real->setFixed(true);
  swave_amp_imag->setValue(0.056949170);//1.96424e-01; 
  swave_amp_imag->setFixed(true);
  if (fixAmp) { swave_amp_real->setValue( 1.); swave_amp_imag->setValue( 0.); swave_amp_real->setFixed(true); swave_amp_imag->setFixed(true); }
  cout<<"Numbers loaded: "<<HH_bin_limits.size()<<" / "<<i<<endl;
  ResonancePdf * swave = new Resonances::Spline("swave", *swave_amp_real,*swave_amp_imag, HH_bin_limits, pwa_coefs_reals, pwa_coefs_imags,
          PAIR_12, true);
  return swave;
}

void getToyData(std::string toyFileName, GooFit::Application &app, DataSet &data) {
    toyFileName = app.get_filename(toyFileName, "examples/dalitz");

    auto obs               = data.getObservables();
    Observable m12         = obs.at(0);
    Observable m13         = obs.at(1);
    Observable eventNumber = obs.at(2);

    TH2F dalitzplot("dalitzplot",
                    "Original Data",
                    m12.getNumBins(),
                    m12.getLowerLimit(),
                    m12.getUpperLimit(),
                    m13.getNumBins(),
                    m13.getLowerLimit(),
                    m13.getUpperLimit());
    std::vector<Observable> vars;
    vars.push_back(m12);
    vars.push_back(m13);
    vars.push_back(eventNumber);

    const string suffix = ".root";
    if (toyFileName.rfind(suffix)+suffix.length() == toyFileName.length()){
      std::cout<<"Reading file "<<toyFileName<<std::endl;
      TFile*f = TFile::Open(toyFileName.c_str());
      TTree*t = (TTree*)f->Get("DecayTree");
      assert(t);
      std::cout<<"Entries: "<<t->GetEntries()<<std::endl;
      double m2_12, m2_13, nsigwt;
      t->SetBranchAddress("s12", &m2_12);
      t->SetBranchAddress("s13", &m2_13);
//      t->SetBranchAddress("nsig_sw", &nsigwt);
/*      const double vmin = (veto_min->value);
      const double vmax = (veto_max->value);*/
      for (int i=0;i<t->GetEntries()/*&&i<MAXEVT*/;i++){
          t->GetEntry(i);
//          if (i%10) continue;// Only accept 10% of data
/*          if (m2_12>vmin&&m2_12<vmax) continue;
          if (m2_13>vmin&&m2_13<vmax) continue;*/
          m12.setValue( m2_12);
          m13.setValue( m2_13);
          eventNumber.setValue(data.getNumEvents()); 
//          wSig0->value = nsigwt;//1-toy_bkg_fraction;
          data.addEvent(); 
          dalitzplot.Fill(m12.getValue(), m13.getValue()); 
      }
      f->Close();
    }
    else{
    std::ifstream reader(toyFileName);
    std::string buffer;

    while(reader >> buffer) {
        if(buffer == "====")
            break;
        std::cout << buffer;
    }

    double dummy = 0;

    while(reader >> dummy
          // m23, m(pi+ pi-), called m12 in processToyRoot convention.
          // Already swapped according to D* charge. m12 = m(pi+pi0)
          >> dummy >> m12 >> m13

          // Errors on Dalitz variables
          >> dummy >> dummy >> dummy

          // Decay time, sigma_t
          >> dummy >> dummy

          // Md0, deltaM, Prob, Sig, Dst charge, Run, Event, Signal and four bkg fractions
          >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy) {
        // EXERCISE 1 (preliminary): Impose an artificial reconstruction efficiency
        // by throwing out events with a probability linear in m12.
        // NB! This exercise continues below.

        // EXERCISE 2: Instead of the above efficiency, impose a
        // K0 veto, by throwing out events with 0.475 < m23 < 0.505.

        // EXERCISE 3: Use both the above.

        eventNumber.setValue(data.getNumEvents());
        data.addEvent();

        dalitzplot.Fill(m12.getValue(), m13.getValue());
    }
    }
    GOOFIT_INFO("Read in {} events", data.getNumEvents());

    TCanvas foo;
    dalitzplot.SetStats(false);
    dalitzplot.Draw("colz");
    foo.SaveAs("dalitzplot.png");
}

void makeToyData(DalitzPlotter &dplotter, UnbinnedDataSet &data) {}

Amp3Body *makeSignalPdf(Observable m12, Observable m13, EventNumber eventNumber, GooPdf *eff = 0) {
    DecayInfo3 dtop0pp;
    dtop0pp.motherMass   = _mDs;
    dtop0pp.daug1Mass    = piPlusMass;
    dtop0pp.daug2Mass    = piPlusMass;
    dtop0pp.daug3Mass    = piPlusMass;
    dtop0pp.meson_radius = 1.5;

    ResonancePdf *f2_1270 = new Resonances::RBW(
        "f2_1270", Variable("f2_1270_amp_real", 1), Variable("f2_1270_amp_imag", 0), 
                              Variable("f2_1270_mass", 1.2755, 0.01, 1.0, 1.5),
                              Variable("f2_1270_width", 0.1867, 0.01, 0.01, 0.4),
        2, PAIR_12,1);


    bool fixAmps = false; // Takes ~400x longer

    ResonancePdf *rho0 = new Resonances::RBW(
        "rho0",
        fixAmps ? Variable("rho0_amp_real", 0.565) : Variable("rho0_amp_real", 0.565, 0.001, 0, 0),
        fixAmps ? Variable("rho0_amp_imag", 0.164) : Variable("rho0_amp_imag", 0.164, 0.1, 0, 0),
        fixedRhoMass,
        fixedRhoWidth,
        1,
        PAIR_12, 0);

/*    ResonancePdf *rho0_1450 = new Resonances::RBW(
        "rho0_1450",
        fixAmps ? Variable("rho0_1450_amp_real", 0.325) : Variable("rho0_1450_amp_real", 0.325, 0.001, 0, 0),
        fixAmps ? Variable("rho0_1450_amp_imag", 0.057) : Variable("rho0_1450_amp_imag", 0.057, 0.1, 0, 0),
        sharedMass,
        shareWidth,
        1,
        PAIR_12, 1);*/


    ResonancePdf *nonr = new Resonances::NonRes(
        "nonr",
        fixAmps ? Variable("nonr_amp_real", 0.5595 * (-1)) : Variable("nonr_amp_real", 0.5595 * (-1), 0.001, 0, 0),
        fixAmps ? Variable("nonr_amp_imag", -0.108761 * (-1)) : Variable("nonr_amp_imag", -0.108761 * (-1), 0.1, 0, 0));

    ResonancePdf* swave = loadPWAResonance("PiPiPi_PWA_COEFFS_30_SAda_ampphs.txt", fixAmps);

//    dtop0pp.resonances.push_back(nonr);
    dtop0pp.resonances.push_back(swave);
//    dtop0pp.resonances.push_back(rho0_1450);
//    dtop0pp.resonances.push_back(f2_1270);

    bool fitMasses = true;

    if(!fitMasses) {
        for(vector<ResonancePdf *>::iterator res = dtop0pp.resonances.begin(); res != dtop0pp.resonances.end(); ++res) {
            (*res)->setParameterConstantness(true);
        }
    }

    if(!eff) {
        // By default create a constant efficiency.
        vector<Variable> offsets;
        vector<Observable> observables;
        vector<Variable> coefficients;

        observables.push_back(m12);
        observables.push_back(m13);
        offsets.push_back(constantZero);
        offsets.push_back(constantZero);
        coefficients.push_back(constantOne);
        eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0);
    }

    return new Amp3Body("signalPDF", m12, m13, eventNumber, dtop0pp, eff);
}

int runToyFit(Amp3Body *signal, UnbinnedDataSet *data) {
    // EXERCISE 1 (real part): Create a PolynomialPdf which models
    // the efficiency you imposed in the preliminary, and use it in constructing
    // the signal PDF.

    // EXERCISE 2: Create a K0 veto function and use it as the efficiency.

    // EXERCISE 3: Make the efficiency a product of the two functions
    // from the previous exercises.

    signal->setData(data);
    signal->setDataSize(data->getNumEvents());
    FitManager datapdf(signal);

//    datapdf.fit();

    ProdPdf prodpdf{"prodpdf", {signal}};

    DalitzPlotter plotter(&prodpdf, signal);

    TCanvas foo("foo", "", 1200,600);
    foo.Divide(2);
    TH2F *dalitzplot = plotter.make2D();
    foo.cd(1);
    dalitzplot->Draw("colz");
    foo.cd(2);
    dalitzplot->ProjectionX()->Draw("colz");

    foo.SaveAs("dalitzpdf.png");

    return datapdf;
}

int main(int argc, char **argv) {
    GooFit::Application app("Dalitz example", argc, argv);

    std::string filename = "dalitz_toyMC_000.txt";
    app.add_option("-f,--filename,filename", filename, "File to read in", true)->check(GooFit::ExistingFile);

    bool make_toy;
    app.add_flag("-m,--make-toy", make_toy, "Make a toy instead of reading a file in");

    GOOFIT_PARSE(app);

    GooFit::setROOTStyle();

    // Observables setup
    Observable m12("m12", 0, 3.5);
    Observable m13("m13", 0, 3.5);
    EventNumber eventNumber("eventNumber");
    m12.setNumBins(350);
    m13.setNumBins(350);

    // Prepare the data
    UnbinnedDataSet data({m12, m13, eventNumber});

    // Set up the model
    Amp3Body *signal = makeSignalPdf(m12, m13, eventNumber);

    // A wrapper for plotting without complex number segfault
    ProdPdf prodpdf{"prodpdf", {signal}};

    // Add nice tool for making data or plotting
    DalitzPlotter dplotter{&prodpdf, signal};

    // Read in data
    if(make_toy) {
        dplotter.fillDataSetMC(data, 1000000);
    } else {
        getToyData(filename, app, data);
    }

    try {
        return runToyFit(signal, &data);
    } catch(const std::runtime_error &e) {
        std::cerr << e.what() << std::endl;
        return 7;
    }
}
