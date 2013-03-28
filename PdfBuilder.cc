#include "FunctorBase.hh"
#include "PdfBuilder.hh"
#include "ThrustPdfFunctor.hh" 
#include <cstdio> 
#include <cassert> 
#include <limits> 
#include <typeinfo> 
#include <set>
#include "Variable.hh" 

using namespace std; 

char stringbuffer[100];
int numDummies = 0; 
#ifdef OMP_ON
#pragma omp threadprivate (stringbuffer, numDummies)
#pragma omp threadprivate (host_params)
#endif

#if MINUIT_VERSION == 3
FunctorBase* pdfPointer; 
PdfFunctor* currGlue = 0; 
int numPars = 0; 
#ifdef OMP_ON
#pragma omp threadprivate (pdfPointer)
#pragma omp threadprivate (currGlue)
#pragma omp threadprivate (numPars)
#endif
#endif
#if MINUIT_VERSION == 1
FunctorBase* pdfPointer; 
PdfFunctor* currGlue = 0; 
int numPars = 0; 
#ifdef OMP_ON
#pragma omp threadprivate (pdfPointer)
#pragma omp threadprivate (currGlue)
#pragma omp threadprivate (numPars)
#endif
#endif

void specialTddpPrint (double fun); 

#if MINUIT_VERSION == 2
ROOT::Minuit2::FunctionMinimum* PdfFunctor::fit () {
  host_callnumber = 0; 
  params = new ROOT::Minuit2::MnUserParameters();
  vars.clear();
  pdfPointer->getParameters(vars); 

  numPars = vars.size();
  int maxIndex = 0;
  for (set<Variable*>::iterator i = vars.begin(); i != vars.end(); ++i) {
    if ((*i)->lowerlimit == (*i)->upperlimit) params->Add((*i)->name, (*i)->value, (*i)->error); 
    else params->Add((*i)->name, (*i)->value, (*i)->error, (*i)->lowerlimit, (*i)->upperlimit); 
    if ((*i)->fixed) params->Fix(params->Index((*i)->name)); 

    if (maxIndex < (*i)->getIndex()) maxIndex = (*i)->getIndex();
  }

  numPars = maxIndex+1; 

  migrad = new ROOT::Minuit2::MnMigrad(*this, *params); 
  ROOT::Minuit2::FunctionMinimum* ret = new ROOT::Minuit2::FunctionMinimum((*migrad)()); 

  return ret; 
}

double PdfFunctor::operator () (const vector<double>& pars) const {
  vector<double> gooPars; // Translates from Minuit indexing to GooFit indexing
  gooPars.resize(numPars); 
  int counter = 0; 
  for (set<Variable*>::iterator i = vars.begin(); i != vars.end(); ++i) {
    gooPars[(*i)->index] = pars[counter++]; 
  }

  pdfPointer->copyParams(gooPars); 
  double nll = pdfPointer->calculateNLL();
  host_callnumber++; 

#ifdef PRINTCALLS
  double edm = migrad->State().Edm(); 
  cout.precision(8); 
  cout << "State at call " 
	    << host_callnumber << " : "
	    << nll << " "
	    << edm << " Pars: ";
  set<Variable*> vars;
  pdfPointer->getParameters(vars); 
  for (set<Variable*>::iterator i = vars.begin(); i != vars.end(); ++i) {
    if (0 > (*i)->getIndex()) continue;
    if ((*i)->fixed) continue;
    cout << "(" << (*i)->name << " " << pars[(*i)->getIndex()] << ") "; // migrad->Value((*i)->getIndex()) << ") ";
  }

  //for (vector<double>::const_iterator i = pars.begin(); i != pars.end(); ++i) {
  //cout << (*i) << " ";
  //}
  cout << endl; 
#endif 

  return nll; 
}
#elif MINUIT_VERSION == 3
PdfFunctor::PdfFunctor (FunctorBase* dat) {
  pdfPointer = dat;
  currGlue = this; 
} 

#ifdef OMP_ON
set<Variable*> vars[MAX_THREADS]; 
#else
set<Variable*> vars; 
#endif 

#include "TMinuit.hh" 
void PdfFunctor::fit () {
  host_callnumber = 0; 
  pdfPointer->getParameters(vars); 

  numPars = vars.size();
  fitter = TVirtualFitter::Fitter(0, numPars);

  int maxIndex = 0; 
  int counter = 0; 
  for (set<Variable*>::iterator i = vars.begin(); i != vars.end(); ++i) {
    fitter->SetParameter(counter, (*i)->name.c_str(), (*i)->value, (*i)->error, (*i)->lowerlimit, (*i)->upperlimit);
    if ((*i)->fixed) fitter->FixParameter(counter); 
    counter++; 
    if (maxIndex < (*i)->getIndex()) maxIndex = (*i)->getIndex();
  }

  numPars = maxIndex+1; 
  pdfPointer->copyParams();   

  // Hah! gMinuit is global, we can avoid the annoying wrapper. 
  //gMinuit->fIdbg[2] = 1; // Debug for mnderi 

  fitter->SetFCN(FitFun); 
  fitter->ExecuteCommand("MIGRAD", 0, 0); 
}

void PdfFunctor::getMinuitValues () const {
  int counter = 0; 
#ifdef OMP_ON
  int tid = omp_get_thread_num();
  for (set<Variable*>::iterator i = vars[tid].begin(); i != vars[tid].end(); ++i) {
    (*i)->value = fitter->GetParameter(counter);
    (*i)->error = fitter->GetParError(counter);
    counter++;
  }
#else
  for (set<Variable*>::iterator i = vars.begin(); i != vars.end(); ++i) {
    (*i)->value = fitter->GetParameter(counter);
    (*i)->error = fitter->GetParError(counter);
    counter++;
  }
#endif
}

void FitFun (int &npar, double *gin, double &fun, double *fp, int iflag) { // MINUIT 3 version 
  vector<double> pars; // Translates from Minuit to GooFit indices
  pars.resize(numPars); 
  int counter = 0; 
  for (set<Variable*>::iterator i = vars.begin(); i != vars.end(); ++i) {
    pars[(*i)->getIndex()] = fp[counter++]; 
  }

  pdfPointer->copyParams(pars); 
  fun = pdfPointer->calculateNLL();
  host_callnumber++; 

#ifdef PRINTCALLS
  specialTddpPrint(fun); 
#endif 
}
#else // MINUIT_VERSION is not 2 or 3

#ifdef OMP_ON
set<Variable*> vars[MAX_THREADS]; 
#else
set<Variable*> vars; 
#endif 

PdfFunctor::PdfFunctor (FunctorBase* dat) 
  : minuit(0)
  , overrideCallLimit(-1)
{
  pdfPointer = dat;
  currGlue = this; 
} 

PdfFunctor::~PdfFunctor () {
  if (minuit) delete minuit; 
}

void PdfFunctor::setupMinuit () {
#ifdef OMP_ON
  int tid = omp_get_thread_num(); 
  vars[tid].clear(); 
  pdfPointer->getParameters(vars[tid]); 

  numPars = vars[tid].size();
  if (minuit) delete minuit;
  minuit = new TMinuit(numPars); 
  int maxIndex = 0; 
  int counter = 0; 
  for (set<Variable*>::iterator i = vars[tid].begin(); i != vars[tid].end(); ++i) {
    minuit->DefineParameter(counter, (*i)->name.c_str(), (*i)->value, (*i)->error, (*i)->lowerlimit, (*i)->upperlimit); 
    if ((*i)->fixed) minuit->FixParameter(counter);
    counter++; 
    if (maxIndex < (*i)->getIndex()) maxIndex = (*i)->getIndex();
  }
#else
  vars.clear(); 
  pdfPointer->getParameters(vars); 

  numPars = vars.size();
  if (minuit) delete minuit;
  minuit = new TMinuit(numPars); 
  int maxIndex = 0; 
  int counter = 0; 
  for (set<Variable*>::iterator i = vars.begin(); i != vars.end(); ++i) {
    minuit->DefineParameter(counter, (*i)->name.c_str(), (*i)->value, (*i)->error, (*i)->lowerlimit, (*i)->upperlimit); 
    if ((*i)->fixed) minuit->FixParameter(counter);
    counter++; 
    if (maxIndex < (*i)->getIndex()) maxIndex = (*i)->getIndex();
  }
#endif
  numPars = maxIndex+1; 
  pdfPointer->copyParams();   
  minuit->SetFCN(FitFun); 
}

void PdfFunctor::fit () {
  setupMinuit();
  runMigrad();
}

void PdfFunctor::runMigrad () { 
  assert(minuit);
  host_callnumber = 0;
  if (0 < overrideCallLimit) {
    std::cout << "Calling MIGRAD with call limit " << overrideCallLimit << std::endl; 
    double plist[1];
    plist[0] = overrideCallLimit;
    int err = 0; 
    minuit->mnexcm("MIGRAD", plist, 1, err);
  }
  else minuit->Migrad(); 
}

void PdfFunctor::getMinuitValues () const {
  int counter = 0; 
#ifdef OMP_ON
  int tid = omp_get_thread_num();
  for (set<Variable*>::iterator i = vars[tid].begin(); i != vars[tid].end(); ++i) {
    //(*i)->value = host_params[(*i)->getIndex()];
    minuit->GetParameter(counter++, (*i)->value, (*i)->error);
  }
#else
  for (set<Variable*>::iterator i = vars.begin(); i != vars.end(); ++i) {
    //(*i)->value = host_params[(*i)->getIndex()];
    minuit->GetParameter(counter++, (*i)->value, (*i)->error);
  }
#endif
}

void FitFun(int &npar, double *gin, double &fun, double *fp, int iflag) {
  //cout << "FitFun call " << host_callnumber << endl; 

  vector<double> pars;
  // Notice that npar is number of variable parameters, not total. 
  pars.resize(numPars); 
  int counter = 0; 
#ifdef OMP_ON
  int tid = omp_get_thread_num();
  for (set<Variable*>::iterator i = vars[tid].begin(); i != vars[tid].end(); ++i) {
    pars[(*i)->getIndex()] = fp[counter++]; 
  }
#else
  for (set<Variable*>::iterator i = vars.begin(); i != vars.end(); ++i) {
    if (isnan(fp[counter])) cout << "Variable " << (*i)->name << " " << (*i)->index << " is NaN\n"; 
    pars[(*i)->getIndex()] = fp[counter++]; 
  }
#endif // OMP_ON
  
  pdfPointer->copyParams(pars); 
  fun = pdfPointer->calculateNLL(); 
  host_callnumber++; 

#ifdef PRINTCALLS
  specialTddpPrint(fun); 
#endif 
}
#endif 


void specialTddpPrint (double fun) {
  // Stupid amplitude-fit debugging method. 
  cout << "Function call " << host_callnumber << ": " << fun << "\n";
  currGlue->getMinuitValues();
  int varCount = 1; 
  for (set<Variable*>::iterator v = vars.begin(); v != vars.end(); ++v) {
    if (!(*v)) cout << "Null!" << endl; 
    if ((*v)->fixed) continue; 

    const fptype _mD0 = 1.86484; 
    const fptype _mD02 = _mD0 *_mD0;
    const fptype _mD02inv = 1./_mD02; 
    double stupidSpecialModifier = 1; // Mikhail interprets some of the weights differently. 
    if (((*v)->name == "f0_980_amp_real") || 
	((*v)->name == "f0_980_amp_imag") ||
	((*v)->name == "f0_1370_amp_real") || 
	((*v)->name == "f0_1370_amp_imag") ||
	((*v)->name == "f0_1500_amp_real") || 
	((*v)->name == "f0_1500_amp_imag") ||
	((*v)->name == "f0_1710_amp_real") || 
	((*v)->name == "f0_1710_amp_imag") ||
	((*v)->name == "f0_600_amp_real") || 
	((*v)->name == "f0_600_amp_imag")) stupidSpecialModifier = -_mD02; 
    else if (((*v)->name == "f2_1270_amp_real") ||
	     ((*v)->name == "f2_1270_amp_imag")) stupidSpecialModifier = -_mD02inv; 
    else if (((*v)->name == "nonr_amp_real") ||
	     ((*v)->name == "nonr_amp_imag")) stupidSpecialModifier = -1; 

    cout.width(20); 
    cout << (*v)->name;
    cout.setf(ios_base::right,ios_base::adjustfield);
    cout.width(3);
    cout << varCount++;
    cout.setf(ios_base::right,ios_base::adjustfield); cout.precision(8);
    cout << "  ";         cout.width(12);
    cout << (*v)->value / stupidSpecialModifier;
    cout.setf(ios_base::right,ios_base::adjustfield); cout.precision(8);
    cout << "  ";         cout.width(12);
    cout << (*v)->error;
    cout << endl; 
  }

  cout << endl; 
}
