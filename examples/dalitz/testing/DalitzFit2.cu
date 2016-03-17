// System stuff
#include <fstream> 
#include <sys/time.h>
#include <sys/times.h>

// GooFit stuff
#include "Variable.hh" 
#include "PolynomialPdf.hh" 
#include "DalitzPlotPdf.hh" 
#include "DalitzVetoPdf.hh" 
#include "ResonancePdf.hh" 
#include "AddPdf.hh"
#include "ProdPdf.hh"
#include "GooPdf.hh" 
#include "FitManager.hh" 
#include "UnbinnedDataSet.hh"

#include "DPPdf.hh"

using namespace std;

timeval startTime, stopTime, totalTime;
clock_t startCPU, stopCPU; 
tms startProc, stopProc; 
UnbinnedDataSet* data = 0; 

Variable* m12 = 0;
Variable* m13 = 0;
Variable* eventNumber = 0; 
bool fitMasses = false; 
Variable* fixedRhoMass  = new Variable("rho_mass", 0.7758, 0.01, 0.7, 0.8);
Variable* fixedRhoWidth = new Variable("rho_width", 0.1503, 0.01, 0.1, 0.2); 

const fptype _mD0 = 1.86484; 
const fptype _mD02 = _mD0 *_mD0;
const fptype _mD02inv = 1./_mD02; 
const fptype piPlusMass = 0.13957018;
const fptype piZeroMass = 0.1349766;

// Constants used in more than one PDF component. 
Variable* motherM = new Variable("motherM", _mD0);
Variable* chargeM = new Variable("chargeM", piPlusMass);
Variable* neutrlM = new Variable("neutrlM", piZeroMass);
Variable* massSum = new Variable("massSum", _mD0*_mD0 + 2*piPlusMass*piPlusMass + piZeroMass*piZeroMass); // = 3.53481 
Variable* constantOne = new Variable("constantOne", 1); 
Variable* constantZero = new Variable("constantZero", 0); 

GooPdf* kzero_veto = 0; 

fptype cpuGetM23 (fptype massPZ, fptype massPM) {
  return (_mD02 + piZeroMass*piZeroMass + piPlusMass*piPlusMass + piPlusMass*piPlusMass - massPZ - massPM); 
}

bool cpuDalitz (fptype m12, fptype m13, fptype bigM, fptype dm1, fptype dm2, fptype dm3) {
  if (m12 < POW(dm1 + dm2, 2)) return false; // This m12 cannot exist, it's less than the square of the (1,2) particle mass.
  if (m12 > POW(bigM - dm3, 2)) return false;   // This doesn't work either, there's no room for an at-rest 3 daughter. 
  
  // Calculate energies of 1 and 3 particles in m12 rest frame. 
  fptype e1star = 0.5 * (m12 - dm2*dm2 + dm1*dm1) / SQRT(m12); 
  fptype e3star = 0.5 * (bigM*bigM - m12 - dm3*dm3) / SQRT(m12); 

  // Bounds for m13 at this value of m12.
  fptype minimum = POW(e1star + e3star, 2) - POW(SQRT(e1star*e1star - dm1*dm1) + SQRT(e3star*e3star - dm3*dm3), 2);
  if (m13 < minimum) return false;
  fptype maximum = POW(e1star + e3star, 2) - POW(SQRT(e1star*e1star - dm1*dm1) - SQRT(e3star*e3star - dm3*dm3), 2);
  if (m13 > maximum) return false;

  return true; 
}

void getToyData (std::string toyFileName) {
  std::vector<Variable*> vars;
  vars.push_back(m12);
  vars.push_back(m13);
  vars.push_back(eventNumber); 
  data = new UnbinnedDataSet(vars); 

  std::ifstream reader;
  reader.open(toyFileName.c_str()); 
  std::string buffer;
  while (!reader.eof()) {
    reader >> buffer;
    if (buffer == "====") break; 
    std::cout << buffer; 
  }

  double dummy = 0; 
  double m23 = 0;
  while (!reader.eof()) {
    reader >> dummy;
    reader >> m23;      // m23, m(pi+ pi-), called m12 in processToyRoot convention. 
    reader >> m12->value; // Already swapped according to D* charge. m12 = m(pi+pi0)
    reader >> m13->value;

    // Errors on Dalitz variables
    reader >> dummy; 
    reader >> dummy; 
    reader >> dummy; 

    reader >> dummy; // Decay time
    reader >> dummy; // sigma_t

    reader >> dummy; // Md0
    reader >> dummy; // deltaM
    reader >> dummy; // ProbSig
    reader >> dummy; // Dst charge
    reader >> dummy; // Run
    reader >> dummy; // Event
    reader >> dummy; // Signal and four bkg fractions. 
    reader >> dummy; 
    reader >> dummy; 
    reader >> dummy; 
    reader >> dummy; 

    // EXERCISE 1 (preliminary): Impose an artificial reconstruction efficiency
    // by throwing out events with a probability linear in m12. 
    // NB! This exercise continues below. 

    // EXERCISE 2: Instead of the above efficiency, impose a 
    // K0 veto, by throwing out events with 0.475 < m23 < 0.505. 

    // EXERCISE 3: Use both the above. 
    if(!(m23 > 0.475 and m23 < 0.505 ) ){
	    eventNumber->value = data->getNumEvents(); 
	    data->addEvent(); 
  	}
  }
}

GooPdf* makeKzeroVeto () {
  if (kzero_veto) return kzero_veto; 

  VetoInfo* kVetoInfo = new VetoInfo();
  kVetoInfo->cyclic_index = PAIR_23; 
  kVetoInfo->minimum = new Variable("veto_min", 0.475*0.475);
  kVetoInfo->maximum = new Variable("veto_max", 0.505*0.505);
  vector<VetoInfo*> vetos; vetos.push_back(kVetoInfo); 
  kzero_veto = new DalitzVetoPdf("kzero_veto", m12, m13, motherM, neutrlM, chargeM, chargeM, vetos); 
  return kzero_veto;
}

DalitzPlotPdf* makeSignalPdf (GooPdf* eff = 0) {
  DecayInfo* dtop0pp = new DecayInfo();
  dtop0pp->motherMass  = _mD0; 
  dtop0pp->daug1Mass  = piZeroMass;
  dtop0pp->daug2Mass  = piPlusMass;
  dtop0pp->daug3Mass  = piPlusMass;
  dtop0pp->meson_radius  = 1.5; 
 
  ResonancePdf* rhop  = new ResonancePdf("rhop",
							     new Variable("rhop_amp_real", 1),
							     new Variable("rhop_amp_imag", 0),
							     fixedRhoMass,
							     fixedRhoWidth,
							     1,
							     PAIR_12);


  bool fixAmps = false;

  ResonancePdf* rhom  = new ResonancePdf("rhom", 
							     fixAmps ? new Variable("rhom_amp_real", 0.714) : 
							     new Variable("rhom_amp_real",  0.714, 0.001, 0, 0),
							     fixAmps ? new Variable("rhom_amp_imag", -0.025) :
							     new Variable("rhom_amp_imag", -0.025, 0.1, 0, 0),
							     fixedRhoMass,
							     fixedRhoWidth,
							     1,
							     PAIR_13);

  ResonancePdf* rho0  = new ResonancePdf("rho0", 
							     fixAmps ? new Variable("rho0_amp_real", 0.565) : 
							     new Variable("rho0_amp_real", 0.565, 0.001, 0, 0),
							     fixAmps ? new Variable("rho0_amp_imag", 0.164) :
							     new Variable("rho0_amp_imag", 0.164, 0.1, 0, 0),
							     fixedRhoMass,
							     fixedRhoWidth,
							     1,
							     PAIR_23);

  Variable* sharedMass = new Variable("rhop_1450_mass", 1.465, 0.01, 1.0, 2.0);
  Variable* shareWidth = new Variable("rhop_1450_width", 0.400, 0.01, 0.01, 5.0); 

  ResonancePdf* rhop_1450  = new ResonancePdf("rhop_1450", 
								  fixAmps ? new Variable("rhop_1450_amp_real", -0.174) : 
								  new Variable("rhop_1450_amp_real", -0.174, 0.001, 0, 0),
								  fixAmps ? new Variable("rhop_1450_amp_imag", -0.117) :
								  new Variable("rhop_1450_amp_imag", -0.117, 0.1, 0, 0),
								  sharedMass,
								  shareWidth,
								  1,
								  PAIR_12);

  ResonancePdf* rho0_1450  = new ResonancePdf("rho0_1450", 
								  fixAmps ? new Variable("rho0_1450_amp_real", 0.325) : 
								  new Variable("rho0_1450_amp_real", 0.325, 0.001, 0, 0),
								  fixAmps ? new Variable("rho0_1450_amp_imag", 0.057) : 
								  new Variable("rho0_1450_amp_imag", 0.057, 0.1, 0, 0),  
								  sharedMass,
								  shareWidth,
								  1,
								  PAIR_23);

  ResonancePdf* rhom_1450  = new ResonancePdf("rhom_1450", 
								  fixAmps ? new Variable("rhom_1450_amp_real", 0.788) : 
								  new Variable("rhom_1450_amp_real", 0.788, 0.001, 0, 0),
								  fixAmps ? new Variable("rhom_1450_amp_imag", 0.226) : 
								  new Variable("rhom_1450_amp_imag", 0.226, 0.1, 0, 0),  
								  sharedMass,
								  shareWidth,
								  1,
								  PAIR_13);

  sharedMass = new Variable("rhop_1700_mass",  1.720, 0.01, 1.6, 1.9);
  shareWidth = new Variable("rhop_1700_width", 0.250, 0.01, 0.1, 1.0); 

  
  ResonancePdf* rhop_1700  = new ResonancePdf("rhop_1700", 
								  fixAmps ? new Variable("rhop_1700_amp_real", 2.151) : 
								  new Variable("rhop_1700_amp_real",  2.151, 0.001, 0, 0),
								  fixAmps ? new Variable("rhop_1700_amp_imag", -0.658) : 
								  new Variable("rhop_1700_amp_imag", -0.658, 0.1, 0, 0),  
								  sharedMass,
								  shareWidth,
								  1,
								  PAIR_12);
  
  ResonancePdf* rho0_1700  = new ResonancePdf("rho0_1700", 
								  fixAmps ? new Variable("rho0_1700_amp_real",  2.400) : 
								  new Variable("rho0_1700_amp_real",  2.400, 0.001, 0, 0),
								  fixAmps ? new Variable("rho0_1700_amp_imag", -0.734) : 
								  new Variable("rho0_1700_amp_imag", -0.734, 0.1, 0, 0),  
								  sharedMass,
								  shareWidth,
								  1,
								  PAIR_23);
  
  ResonancePdf* rhom_1700  = new ResonancePdf("rhom_1700", 
								  fixAmps ? new Variable("rhom_1700_amp_real",  1.286) : 
								  new Variable("rhom_1700_amp_real",  1.286, 0.001, 0, 0),
								  fixAmps ? new Variable("rhom_1700_amp_imag", -1.532) : 
								  new Variable("rhom_1700_amp_imag", -1.532, 0.1, 0, 0),  
								  sharedMass,
								  shareWidth,
								  1,
								  PAIR_13);
  
  ResonancePdf* f0_980  = new ResonancePdf("f0_980", 
							       fixAmps ? new Variable("f0_980_amp_real",  0.008 * (-_mD02)) : 
							       new Variable("f0_980_amp_real",  0.008 * (-_mD02), 0.001, 0, 0),
							       fixAmps ? new Variable("f0_980_amp_imag", -0.013 * (-_mD02)) : 
							       new Variable("f0_980_amp_imag", -0.013 * (-_mD02), 0.1, 0, 0),  
							       new Variable("f0_980_mass",     0.980, 0.01, 0.8, 1.2),
							       new Variable("f0_980_width",    0.044, 0.001, 0.001, 0.08),
							       0,
							       PAIR_23);
  
  ResonancePdf* f0_1370  = new ResonancePdf("f0_1370", 
								fixAmps ? new Variable("f0_1370_amp_real", -0.058 * (-_mD02)) : 
								new Variable("f0_1370_amp_real", -0.058 * (-_mD02), 0.001, 0, 0),
								fixAmps ? new Variable("f0_1370_amp_imag",  0.026 * (-_mD02)) : 
								new Variable("f0_1370_amp_imag",  0.026 * (-_mD02), 0.1, 0, 0),  
								new Variable("f0_1370_mass",     1.434, 0.01, 1.2, 1.6),
								new Variable("f0_1370_width",    0.173, 0.01, 0.01, 0.4),
								0,
								PAIR_23);
  
  ResonancePdf* f0_1500  = new ResonancePdf("f0_1500", 
								fixAmps ? new Variable("f0_1500_amp_real", 0.057 * (-_mD02)) : 
								new Variable("f0_1500_amp_real", 0.057 * (-_mD02), 0.001, 0, 0),
								fixAmps ? new Variable("f0_1500_amp_imag", 0.012 * (-_mD02)) : 
								new Variable("f0_1500_amp_imag", 0.012 * (-_mD02), 0.1, 0, 0),  
								new Variable("f0_1500_mass",     1.507, 0.01, 1.3, 1.7),
								new Variable("f0_1500_width",    0.109, 0.01, 0.01, 0.3),
								0,
								PAIR_23);
  
  ResonancePdf* f0_1710  = new ResonancePdf("f0_1710", 
								fixAmps ? new Variable("f0_1710_amp_real", 0.070 * (-_mD02)) : 
								new Variable("f0_1710_amp_real", 0.070 * (-_mD02), 0.001, 0, 0),
								fixAmps ? new Variable("f0_1710_amp_imag", 0.087 * (-_mD02)) : 
								new Variable("f0_1710_amp_imag", 0.087 * (-_mD02), 0.1, 0, 0),  
								new Variable("f0_1710_mass",     1.714, 0.01, 1.5, 2.9), 
								new Variable("f0_1710_width",    0.140, 0.01, 0.01, 0.5),
								0,
								PAIR_23);
  
  ResonancePdf* f2_1270  = new ResonancePdf("f2_1270", 
								fixAmps ? new Variable("f2_1270_amp_real", -1.027 * (-_mD02inv)) : 
								new Variable("f2_1270_amp_real", -1.027 * (-_mD02inv), 0.001, 0, 0),
								fixAmps ? new Variable("f2_1270_amp_imag", -0.162 * (-_mD02inv)) : 
								new Variable("f2_1270_amp_imag", -0.162 * (-_mD02inv), 0.1, 0, 0),  
								new Variable("f2_1270_mass",     1.2754, 0.01, 1.0, 1.5),
								new Variable("f2_1270_width",    0.1851, 0.01, 0.01, 0.4),
								2,
								PAIR_23);
  
  ResonancePdf* f0_600  = new ResonancePdf("f0_600", 
							       fixAmps ? new Variable("f0_600_amp_real", 0.068 * (-_mD02)) : 
							       new Variable("f0_600_amp_real", 0.068 * (-_mD02), 0.001, 0, 0),
							       fixAmps ? new Variable("f0_600_amp_imag", 0.010 * (-_mD02)) : 
							       new Variable("f0_600_amp_imag", 0.010 * (-_mD02), 0.1, 0, 0),  
							       new Variable("f0_600_mass",     0.500, 0.01, 0.3, 0.7),
							       new Variable("f0_600_width",    0.400, 0.01, 0.2, 0.6), 
							       0,
							       PAIR_23);
  
  ResonancePdf* nonr  = new ResonancePdf("nonr",
							     fixAmps ? new Variable("nonr_amp_real", 0.5595 * (-1)) : 
							     new Variable("nonr_amp_real", 0.5595 * (-1),   0.001, 0, 0),
							     fixAmps ? new Variable("nonr_amp_imag", -0.108761 * (-1)) : 
							     new Variable("nonr_amp_imag", -0.108761* (-1), 0.1, 0, 0)); 

//  dtop0pp->resonances.push_back(nonr); 
  dtop0pp->resonances.push_back(rhop);
  dtop0pp->resonances.push_back(rhop);
//  dtop0pp->resonances.push_back(rho0); 
//  dtop0pp->resonances.push_back(rhom); 
//  dtop0pp->resonances.push_back(rhop_1450); 
//  dtop0pp->resonances.push_back(rho0_1450); 
//  dtop0pp->resonances.push_back(rhom_1450); 
//  dtop0pp->resonances.push_back(rhop_1700); 
//  dtop0pp->resonances.push_back(rho0_1700); 
//  dtop0pp->resonances.push_back(rhom_1700); 
//  dtop0pp->resonances.push_back(f0_980); 
//  dtop0pp->resonances.push_back(f0_1370); 
//  dtop0pp->resonances.push_back(f0_1500); 
//  dtop0pp->resonances.push_back(f0_1710); 
//  dtop0pp->resonances.push_back(f2_1270); 
//  dtop0pp->resonances.push_back(f0_600); 

  if (!fitMasses) {
    for (vector<ResonancePdf*>::iterator res = dtop0pp->resonances.begin(); res != dtop0pp->resonances.end(); ++res) {
      (*res)->setParameterConstantness(true); 
    }
  }

  if (!eff) {
    // By default create a constant efficiency. 
    vector<Variable*> offsets;
    vector<Variable*> observables;
    vector<Variable*> coefficients; 

    observables.push_back(m12);
    observables.push_back(m13);
    offsets.push_back(constantZero);
    offsets.push_back(constantZero);
    coefficients.push_back(constantOne); 
    eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0);
  }
  return new DalitzPlotPdf("signalPDF", m12, m13, eventNumber, dtop0pp, eff);
  
}

void runToyFit (std::string toyFileName) {
  m12 = new Variable("m12", 0, 3);
  m13 = new Variable("m13", 0, 3); 
  m12->numbins = 240;
  m13->numbins = 240;
  eventNumber = new Variable("eventNumber", 0, INT_MAX);
  getToyData(toyFileName);

  // EXERCISE 1 (real part): Create a PolynomialPdf which models
  // the efficiency you imposed in the preliminary, and use it in constructing
  // the signal PDF. 

  // EXERCISE 2: Create a K0 veto function and use it as the efficiency. 
  // EXERCISE 3: Make the efficiency a product of the two functions
  // from the previous exercises.


  DalitzPlotPdf* dalitzpdf = makeSignalPdf(); 
  dalitzpdf->setDataSize(data->getNumEvents()); 
  std::vector<Variable*> vars;
  Variable* constant = new Variable("constant", 0.1); 
  Variable* constant2 = new Variable("constant", 1.0); 
  vars.push_back(constant);
  PolynomialPdf backgr("backgr", m12, vars); 
  AddPdf* signal = new AddPdf("signal",constant2,dalitzpdf, &backgr);

  signal->setData(data); 
  FitManager datapdf(signal); 
  
  gettimeofday(&startTime, NULL);
  startCPU = times(&startProc);
  // datapdf.fit(); 
  stopCPU = times(&stopProc);
  gettimeofday(&stopTime, NULL);

  vars.clear();
  vars.push_back(m12);
  vars.push_back(m13);
  vars.push_back(eventNumber); 
  UnbinnedDataSet currData(vars); 
  int evtCounter = 0; 


  for (int i = 0; i < m12->numbins; ++i) {
    m12->value = m12->lowerlimit + (m12->upperlimit - m12->lowerlimit)*(i + 0.5) / m12->numbins; 
    for (int j = 0; j < m13->numbins; ++j) {
      m13->value = m13->lowerlimit + (m13->upperlimit - m13->lowerlimit)*(j + 0.5) / m13->numbins; 
      if (!cpuDalitz (m12->value, m13->value, _mD0, piZeroMass, piPlusMass, piPlusMass)) continue;
    eventNumber->value = evtCounter; 
    evtCounter++;
    currData.addEvent(); 
    }
  }
  signal->setData(&currData);
  dalitzpdf->setDataSize(currData.getNumEvents()); 
  std::vector<std::vector<double> > pdfValues;
  signal->getCompProbsAtDataPoints(pdfValues);
  for (unsigned int j = 0; j < pdfValues[0].size(); ++j) {
    double currm12 = currData.getValue(m12, j);
    double currm13 = currData.getValue(m13, j);
  }

  std::cout << pdfValues[0].size() << std::endl;

}

int main (int argc, char** argv) {
  cudaSetDevice(0);
  // runToyFit(argv[1]);

  DecayInfo_DP* dtop0pp = new DecayInfo_DP();
  dtop0pp->meson_radius =1.5;
  dtop0pp->particle_masses[211] = _mD0;
  dtop0pp->particle_masses[212] = piZeroMass;
  dtop0pp->particle_masses[213] = piPlusMass;
  dtop0pp->particle_masses[214] = piPlusMass;
  dtop0pp->particle_masses[215] = 139.1;

  unsigned int a = 333;
  std::map<std::string, Lineshape*> LS;
  LS["12311"]= new Lineshape("rhop",
                   211,
                   fixedRhoMass,
                   fixedRhoWidth,
                   1,
                   PAIR_12);
  // LS["12345"]= new Lineshape("rho12",
  //                  211,
  //                  fixedRhoMass,
  //                  fixedRhoWidth,
  //                  1,
  //                  PAIR_12);
  // LS["1211"]= new Lineshape("rhop",
  //                  211,
  //                  fixedRhoMass,
  //                  fixedRhoWidth,
  //                  1,
  //                  PAIR_12);
  std::map<std::string, SpinFactor*> SC;
  unsigned int svsp = 211;
  SC["0101"] = new SpinFactor("spin0", 0);
  SC["0000"] = new SpinFactor("spin0", 0);

  Amplitude* test = new Amplitude( "decazstring", new Variable("amp_real", 1),new Variable("amp_imag", 0), LS, SC);
  Amplitude* test2 = new Amplitude( "decazstring2", new Variable("amp_real", 1),new Variable("amp_imag", 0), LS, SC);

  dtop0pp->amplitudes.push_back(test);
  dtop0pp->amplitudes.push_back(test2);


  m12 = new Variable("m12", 0, 3);
  m13 = new Variable("m13", 0, 3); 
  m12->numbins = 10;
  m13->numbins = 10;
  vector<Variable*> observables;
  vector<Variable*> coefficients; 
  vector<Variable*> offsets;

  observables.push_back(m12);
  observables.push_back(m13);
  offsets.push_back(constantZero);
  offsets.push_back(constantZero);
  coefficients.push_back(constantOne); 
  PolynomialPdf* eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0);

  eventNumber = new Variable("eventNumber", 0, INT_MAX);


  DalitzPlotPdf* dalitzpdf = makeSignalPdf(); 
  DPPdf* dp = new DPPdf("test", observables, eventNumber, dtop0pp, eff);

  std::vector<Variable*> vars;
  Variable* constant = new Variable("constant", 0.1); 
  Variable* constant2 = new Variable("constant", 1.0); 
  vars.push_back(constant);
  PolynomialPdf backgr("backgr", m12, vars); 
  // AddPdf* signal = new AddPdf("signal",constant2,dp, &backgr);
  AddPdf* signal = new AddPdf("signal",constant2,dalitzpdf, &backgr);

  

  vars.clear();
  vars.push_back(m12);
  vars.push_back(m13);
  vars.push_back(eventNumber); 
  UnbinnedDataSet currData(vars); 
  int evtCounter = 0; 


  for (int i = 0; i < m12->numbins; ++i) {
    m12->value = m12->lowerlimit + (m12->upperlimit - m12->lowerlimit)*(i + 0.5) / m12->numbins; 
    for (int j = 0; j < m13->numbins; ++j) {
      m13->value = m13->lowerlimit + (m13->upperlimit - m13->lowerlimit)*(j + 0.5) / m13->numbins; 
      if (!cpuDalitz (m12->value, m13->value, _mD0, piZeroMass, piPlusMass, piPlusMass)) continue;
    eventNumber->value = evtCounter; 
    evtCounter++;
    currData.addEvent(); 
    }
  }


//  std::vector<double> v;
//  dp->evaluateAtPoints(m12,v);
//  for (int i = 0; i < v.size(); ++i)
//  {
//    std::cout<<v[i]<<" "<<std::endl;
//  }

  signal->setData(&currData);
  dp->setDataSize(currData.getNumEvents()); 
  dalitzpdf->setDataSize(currData.getNumEvents());
  std::vector<std::vector<double> > pdfValues;
  signal->getCompProbsAtDataPoints(pdfValues);

  std::cout << pdfValues[0].size() << "  " << pdfValues[0][0] <<std::endl;
  // Print total minimization time
  double myCPU = stopCPU - startCPU;
  double totalCPU = myCPU; 

  timersub(&stopTime, &startTime, &totalTime);
  std::cout << "Wallclock time  : " << totalTime.tv_sec + totalTime.tv_usec/1000000.0 << " seconds." << std::endl;
  std::cout << "CPU time: " << (myCPU / CLOCKS_PER_SEC) << std::endl; 
  std::cout << "Total CPU time: " << (totalCPU / CLOCKS_PER_SEC) << std::endl; 
  myCPU = stopProc.tms_utime - startProc.tms_utime;
  std::cout << "Processor time: " << (myCPU / CLOCKS_PER_SEC) << std::endl;

  return 0; 
}
