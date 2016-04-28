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

#include <thrust/count.h>

using namespace std;

// Constants used in more than one PDF component. 
const fptype _mD0 = 1.8645; 
const fptype piPlusMass = 0.13957018;
const fptype KmMass = .493677;

int main (int argc, char** argv) {

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
  LSKRS["rho(770)"]    = new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_12, LS::BW);
  LSKRS["K*(892)bar"]  = new Lineshape("K*(892)bar", KstarM, KstarW, 1, M_34, LS::BW);
  LSKRS["BSrho(770)"]  = new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_24, LS::BW);
  LSKRS["BSK*(892)bar"]= new Lineshape("K*(892)bar", KstarM, KstarW, 1, M_13, LS::BW);

  std::map<std::string, Lineshape*> LSKRP;
  LSKRP["rho(770)"]    = new Lineshape("rho(770)P", RhoMass, RhoWidth, 1, M_12, LS::BW);
  LSKRP["K*(892)bar"]  = new Lineshape("K*(892)barP", KstarM, KstarW, 1, M_34, LS::BW);
  LSKRP["BSrho(770)"]  = new Lineshape("rho(770)P", RhoMass, RhoWidth, 1, M_24, LS::BW);
  LSKRP["BSK*(892)bar"]= new Lineshape("K*(892)barP", KstarM, KstarW, 1, M_13, LS::BW);

  std::map<std::string, Lineshape*> LSKRD;
  LSKRD["rho(770)"]    = new Lineshape("rho(770)D", RhoMass, RhoWidth, 1, M_12, LS::BW);
  LSKRD["K*(892)bar"]  = new Lineshape("K*(892)barD", KstarM, KstarW, 1, M_34, LS::BW);
  LSKRD["BSrho(770)"]  = new Lineshape("rho(770)D", RhoMass, RhoWidth, 1, M_24, LS::BW);
  LSKRD["BSK*(892)bar"]= new Lineshape("K*(892)barD", KstarM, KstarW, 1, M_13, LS::BW);

  // the very last parameter means that we have two permutations. so the first half of the Lineshapes 
  // and the first half of the spinfactors are amplitude 1, rest is amplitude 2
  // This means that it is important for symmetrized amplitueds that the spinfactors and lineshapes are in the "right" order
  Amplitude* Bose_symmetrized_AMP_S = new Amplitude( "K*(892)rho(770)_S", new Variable("amp_real1", -0.115177), new Variable("amp_imag1", 0.153976), LSKRS, SFKRS, 2);
  Amplitude* Bose_symmetrized_AMP_P = new Amplitude( "K*(892)rho(770)_P", new Variable("amp_real2", -0.0298697), new Variable("amp_imag2", -0.0722874), LSKRP, SFKRP, 2);
  Amplitude* Bose_symmetrized_AMP_D = new Amplitude( "K*(892)rho(770)_D", new Variable("amp_real3", -0.452212), new Variable("amp_imag3", 0.426521), LSKRD, SFKRD, 2);

  DK3P_DI->amplitudes.push_back(Bose_symmetrized_AMP_S);
  DK3P_DI->amplitudes.push_back(Bose_symmetrized_AMP_P);
  DK3P_DI->amplitudes.push_back(Bose_symmetrized_AMP_D);


  Variable* m12 = new Variable("m12", 0, 3);
  Variable* m34 = new Variable("m34", 0, 3); 
  Variable* cos12 = new Variable("cos12", -1, 1);
  Variable* cos34 = new Variable("m12", -1, 1);
  Variable* phi = new Variable("phi", -3.5, 3.5);
  Variable* eventNumber = new Variable("eventNumber", 0, INT_MAX);
  Variable* constantOne = new Variable("constantOne", 1); 
  Variable* constantZero = new Variable("constantZero", 0);
 
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
  DPPdf* dp = new DPPdf("test", observables, DK3P_DI, eff,2e6);


  int numEvents = 1e6;
  auto tuple = dp->GenerateSig(numEvents);
  
  auto variables = std::get<1>(tuple);
  auto weights = std::get<2>(tuple);
  auto flags = std::get<3>(tuple);
  int accepted = thrust::count_if(flags.begin(), flags.end(), thrust::identity<bool>());
  fprintf(stderr,"Using accept-reject method would leave you with %i out of %i events\n", accepted, numEvents);

  for (int i = 0; i < weights.size(); ++i)
  {
    printf("%.5g %.5g %.5g %.5g %.5g %.5g %.5g\n", (*(variables[0]))[i], (*(variables[1]))[i], (*(variables[2]))[i], (*(variables[3]))[i], (*(variables[4]))[i], weights[i], flags[i]);
  }

  return 0; 
}
