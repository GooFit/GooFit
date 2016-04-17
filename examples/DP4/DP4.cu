//ROOT
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

// GooFit stuff
#include "Variable.hh" 
#include "PolynomialPdf.hh" 
#include "AddPdf.hh"
#include "UnbinnedDataSet.hh"
#include "DP4Pdf.hh"
using namespace std;

const fptype _mD0 = 1.8645; 
const fptype piPlusMass = 0.13957018;
const fptype KmMass = .493677;
// Constants used in more than one PDF component. 

int main (int argc, char** argv) {

  unsigned int maxevents = 5;
  fptype* hostnorm = new fptype[5*maxevents];
  TFile* f1 = TFile::Open("phspMC.root");
  TTreeReader reader1("t1", f1);
  TTreeReaderValue<double> tm12(reader1, "m12");
  TTreeReaderValue<double> tm34(reader1, "m34");
  TTreeReaderValue<double> tcos12(reader1, "cos12");
  TTreeReaderValue<double> tcos34(reader1, "cos34");
  TTreeReaderValue<double> tphi(reader1, "phi");
  unsigned int MCevents = 0;

  while(MCevents<maxevents && reader1.Next()){
    hostnorm[MCevents*5]    = *tm12;
    hostnorm[1+MCevents*5]  = *tm34;
    hostnorm[2+MCevents*5]  = *tcos12;
    hostnorm[3+MCevents*5]  = *tcos34;
    hostnorm[4+MCevents*5]  = *tphi;
    // printf("%.10g, %.10g, %.10g, %.10g, %.10g \n",*tm12, *tm34, *tcos12, *tcos34, *tphi );
    MCevents++;
  }

  printf("read in %i events\n", MCevents );

  DecayInfo_DP* DK3P_DI = new DecayInfo_DP();
  DK3P_DI->meson_radius =1.5;
  DK3P_DI->particle_masses.push_back(_mD0);
  DK3P_DI->particle_masses.push_back(piPlusMass);
  DK3P_DI->particle_masses.push_back(piPlusMass);
  DK3P_DI->particle_masses.push_back(KmMass);
  DK3P_DI->particle_masses.push_back(piPlusMass);

  Variable* RhoMass  = new Variable("rho_mass", 0.77526, 0.01, 0.7, 0.8);
  Variable* RhoWidth = new Variable("rho_width", 0.1478, 0.01, 0.1, 0.2); 
  Variable* KstarM   = new Variable("KstarM", 0.89581, 0.01, 0.9, 0.1);
  Variable* KstarW   = new Variable("KstarW", 0.0474, 0.01, 0.1, 0.2); 

  bool useMINT = true;

  //Map of spinfactors, here we have two because of the bose symmetrization of the two pi+
  std::map<std::string, SpinFactor*> SFKRS;
  SFKRS["SVVS"]   = new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S, 0, 1, 2, 3);
  SFKRS["BSSVVS"] = new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S, 3, 1, 2, 0);

  std::map<std::string, SpinFactor*> SFKRP;
  SFKRP["SVVP"]   = new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_P, 0, 1, 2, 3);
  SFKRP["BSSVVP"] = new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_P, 3, 1, 2, 0);

  std::map<std::string, SpinFactor*> SFKRD;
  SFKRD["SVVD"]   = new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_D, 0, 1, 2, 3);
  SFKRD["BSSVVD"] = new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_D, 3, 1, 2, 0);

  //Map of lineshapes, also for both pi+ configurations
  std::map<std::string, Lineshape*> LSKRS;
  LSKRS["rho(770)"]    = new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_12, useMINT);
  LSKRS["K*(892)bar"]  = new Lineshape("K*(892)bar", KstarM, KstarW, 1, M_34, useMINT);
  LSKRS["BSrho(770)"]  = new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_24, useMINT);
  LSKRS["BSK*(892)bar"]= new Lineshape("K*(892)bar", KstarM, KstarW, 1, M_13, useMINT);

  std::map<std::string, Lineshape*> LSKRP;
  LSKRP["rho(770)"]    = new Lineshape("rho(770)P", RhoMass, RhoWidth, 1, M_12, useMINT);
  LSKRP["K*(892)bar"]  = new Lineshape("K*(892)barP", KstarM, KstarW, 1, M_34, useMINT);
  LSKRP["BSrho(770)"]  = new Lineshape("rho(770)P", RhoMass, RhoWidth, 1, M_24, useMINT);
  LSKRP["BSK*(892)bar"]= new Lineshape("K*(892)barP", KstarM, KstarW, 1, M_13, useMINT);

  std::map<std::string, Lineshape*> LSKRD;
  LSKRD["rho(770)"]    = new Lineshape("rho(770)D", RhoMass, RhoWidth, 1, M_12, useMINT);
  LSKRD["K*(892)bar"]  = new Lineshape("K*(892)barD", KstarM, KstarW, 1, M_34, useMINT);
  LSKRD["BSrho(770)"]  = new Lineshape("rho(770)D", RhoMass, RhoWidth, 1, M_24, useMINT);
  LSKRD["BSK*(892)bar"]= new Lineshape("K*(892)barD", KstarM, KstarW, 1, M_13, useMINT);

  // the very last parameter means that we have two permutations. so the first half of the Lineshapes 
  // and the first half of the spinfactors are amplitude 1, rest is amplitude 2
  // This means that it is important for symmetrized amplitueds that the spinfactors and lineshapes are in the "right" order
  Amplitude* Bose_symmetrized_AMP_S = new Amplitude( "K*(892)rho(770)_S", new Variable("amp_real1", -0.115177), new Variable("amp_imag1", 0.153976), LSKRS, SFKRS, 2);
  Amplitude* Bose_symmetrized_AMP_P = new Amplitude( "K*(892)rho(770)_P", new Variable("amp_real2", -0.0298697), new Variable("amp_imag2", -0.0722874), LSKRP, SFKRP, 2);
  Amplitude* Bose_symmetrized_AMP_D = new Amplitude( "K*(892)rho(770)_D", new Variable("amp_real3", -0.452212), new Variable("amp_imag3", 0.426521), LSKRD, SFKRD, 2);


  for (std::map<std::string, Lineshape*>::iterator res = LSKRS.begin(); res != LSKRS.end(); ++res) {
    (*res).second->setParameterConstantness(true); 
  }

  for (std::map<std::string, Lineshape*>::iterator res = LSKRP.begin(); res != LSKRP.end(); ++res) {
    (*res).second->setParameterConstantness(true); 
  }

  for (std::map<std::string, Lineshape*>::iterator res = LSKRD.begin(); res != LSKRD.end(); ++res) {
    (*res).second->setParameterConstantness(true); 
  }


  DK3P_DI->amplitudes.push_back(Bose_symmetrized_AMP_S);
  DK3P_DI->amplitudes.push_back(Bose_symmetrized_AMP_P);
  DK3P_DI->amplitudes.push_back(Bose_symmetrized_AMP_D);

  Variable* m12 = new Variable("m12", 0, 3);
  Variable* m34 = new Variable("m34", 0, 3); 
  Variable* cos12 = new Variable("cos12", -1, 1);
  Variable* cos34 = new Variable("m12", -1, 1);
  Variable* phi = new Variable("phi", -3.5, 3.5);
  m12->numbins = 250;
  m34->numbins = 250;
  cos12->numbins = 250;
  cos34->numbins = 250;
  phi->numbins = 250;

  Variable* constantOne = new Variable("constantOne", 1); 
  Variable* constantZero = new Variable("constantZero", 0);
  Variable* eventNumber = new Variable("eventNumber", 0, 99e9);

  vector<Variable*> observables;
  vector<Variable*> coefficients; 
  vector<Variable*> offsets;

  observables.push_back(m12);
  observables.push_back(m34);
  observables.push_back(cos12);
  observables.push_back(cos34);
  observables.push_back(phi);
  observables.push_back(eventNumber);
  offsets.push_back(constantZero);
  offsets.push_back(constantZero);
  coefficients.push_back(constantOne); 

  PolynomialPdf* eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0);
  DPPdf* dp = new DPPdf("test", observables, DK3P_DI, eff,25e5);

  std::vector<Variable*> vars;
  Variable* constant = new Variable("constant", 0.1); 
  Variable* constant2 = new Variable("constant", 1.0); 
  vars.push_back(constant);
  PolynomialPdf backgr("backgr", m12, vars); 
  AddPdf* signal = new AddPdf("signal",constant2,dp, &backgr);

  vars.clear();
  vars.push_back(m12);
  vars.push_back(m34);
  vars.push_back(cos12);
  vars.push_back(cos34);
  vars.push_back(phi);
  vars.push_back(eventNumber); 
  UnbinnedDataSet currData(vars); 
  int evtCounter = 0; 

//use 5e5 phasespace from the normalisation set as events
for (int i = 0; i < 5; ++i)
{ 
  m12->value = hostnorm[i*5]  ;
  m34->value = hostnorm[1+i*5];
  cos12->value = hostnorm[2+i*5];
  cos34->value = hostnorm[3+i*5];
  phi->value = hostnorm[4+i*5];
  eventNumber->value = evtCounter++; 
  currData.addEvent();
  // printf("%.5g %.5g %.5g %.5g %.5g\n",hostnorm[i*5], hostnorm[1+i*5], hostnorm[2+i*5], hostnorm[3+i*5], hostnorm[4+i*5] );
}

  //set phasespace Events for integration
  // dp->setphsp(hostnorm, MCevents);

  signal->setData(&currData);
  dp->setDataSize(currData.getNumEvents(), 6); 

  // FitManager datapdf(signal);
  // datapdf.fit();
  
  std::vector<std::vector<double> > pdfValues;
  signal->getCompProbsAtDataPoints(pdfValues);
  // for (int i = 0; i < 5e5; ++i)
  // {
  //   printf("%.10g\n", pdfValues[0][i]);
  // }


  return 0; 
}
