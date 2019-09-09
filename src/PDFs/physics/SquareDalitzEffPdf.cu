#include <goofit/PDFs/physics/SquareDalitzEffPdf.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <vector>


namespace GooFit {

__device__ fptype inPS(fptype m12, fptype m13, fptype mD, fptype mKS0, fptype mh1, fptype mh2) {

  if (m12 < pow(mKS0 + mh1, 2)) return 0;
  if (m12 > pow(mD - mh2, 2)) return 0;

  // Calculate energies of 1 and 3 particles in m12 rest frame. 
  fptype e1star = 0.5 * (m12 - mh1*mh1 + mKS0*mKS0) / sqrt(m12);
  fptype e3star = 0.5 * (mD*mD - m12 - mh2*mh2) / sqrt(m12);

  fptype minimum = pow(e1star + e3star, 2) - pow(sqrt(e1star*e1star - mKS0*mKS0) + sqrt(e3star*e3star - mh2*mh2), 2);
  if (m13 < minimum) return 0;
  fptype maximum = pow(e1star + e3star, 2) - pow(sqrt(e1star*e1star - mKS0*mKS0) - sqrt(e3star*e3star - mh2*mh2), 2);
  if (m13 > maximum) return 0;

  return 1;
}

__device__ fptype mprime (fptype m12, fptype m13, fptype mD, fptype mKS0, fptype mh1, fptype mh2) {
  // Helper function to calculate m'^2
  fptype m23 = mD*mD + mKS0*mKS0 + mh1*mh1 + mh2*mh2 - m12 - m13; 
  //fptype rootPi = -2.*ATAN2(-1.0,0.0); // Pi

  if (m23 < 0) return -99;
  fptype tmp = ((2.0*(sqrt(m23) - (mh1 + mh2))/(mD - mKS0 - (mh1 + mh2))) - 1.0);
  if (isnan(tmp)) tmp = -99;
  return tmp;
}

__device__ fptype thetaprime (fptype m12, fptype m13, fptype mD, fptype mKS0, fptype mh1, fptype mh2) {
  // Helper function to calculate theta'
  fptype m23 = mD*mD + mKS0*mKS0 + mh1*mh1 + mh2*mh2 - m12 - m13; 
  if (m23 < 0) return -99;

  fptype num = m23*( m12 - m13) + (mh2*mh2 - mh1*mh1)*(mD*mD - mKS0*mKS0);
  fptype denum = sqrt(((m23 - mh1*mh1 + mh2*mh2)*(m23 - mh1*mh1 + mh2*mh2) - 4*m23*mh2*mh2))*sqrt(((mD*mD - mKS0*mKS0 - m23)*(mD*mD - mKS0*mKS0 -m23) - 4*m23*mKS0*mKS0));
  fptype theta = -99 ;
  if (isnan(denum)) return -99;

  if (denum != 0.){
    theta = num/denum;
  }

  return theta;
}

__device__ fptype device_SquareDalitzEff (fptype* evt, ParameterContainer &pc) {

  // Define observables 
  int idx = pc.getObservable(0);
  int idy = pc.getObservable(1);
  int id_num = pc.getObservable(2);

  fptype x = RO_CACHE(evt[idx]);
  fptype y = RO_CACHE(evt[idy]);
  fptype evtIndex = RO_CACHE(evt[id_num]);

  // Define coefficients
  fptype c0 = pc.getParameter(0);
  fptype c1 = pc.getParameter(1);
  fptype c2 = pc.getParameter(2);
  fptype c3 = pc.getParameter(3);
  fptype c4 = pc.getParameter(4);
  fptype c5 = pc.getParameter(5);
  fptype c6 = pc.getParameter(6);

  // Define constvals
  fptype mD = pc.getParameter(7);
  fptype mKS0 = pc.getParameter(8);
  fptype mh1 = pc.getParameter(9);
  fptype mh2 = pc.getParameter(10);

  pc.incrementIndex(1, pc.getNumParameters(), pc.getNumConstants(), pc.getNumObservables(), 1);

  // Check phase space
  if (inPS == 0) return 0;
  
  // Call helper functions
  fptype thetap = thetaprime(x,y,mD,mKS0,mh1,mh2); 
  if (thetap > 1. || thetap < -1.) return 0; 

  fptype m23 = mD*mD + mKS0*mKS0 + mh1*mh1 + mh2*mh2 - x - y; 
  if (m23 < 0) return 0;

  fptype ret = c0*m23*m23 + c1*m23 + c2*m23*thetap*thetap + c3*thetap*thetap + c4*thetap + c5 + c6*m23*thetap; 
  
  return ret; 
}

__device__ device_function_ptr ptr_to_SquareDalitzEff = device_SquareDalitzEff; 

__host__ __device__ SquareDalitzEffPdf::SquareDalitzEffPdf (std::string n, 
				        std::vector<Observable> obses, 
					std::vector<Variable> coeffs, 
					std::vector<Variable> constvals) 
  : GooPdf("SquareDalitzEffPdf", n, obses, coeffs, constvals) {

  // Register observables - here m12, m13 and dtime
  for (unsigned int i = 0; i < obses.size(); ++i) {
    registerObservable(obses[i]);
  }

  // Register constvals
  for (std::vector<Variable>::iterator v = constvals.begin(); v != constvals.end(); ++v) {
    registerParameter(v);
  }

  // Register coefficients
  for (std::vector<Variable>::iterator c = coeffs.begin(); c != coeffs.end(); ++c) {
    registerParameter(c);
  }

  registerFunction("ptr_to_SquareDalitzEff", ptr_to_SquareDalitzEff);

  initialize();

}

} // namespace GooFit