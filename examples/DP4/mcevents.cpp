#include <iostream>
#include <cmath>
#include <TGenPhaseSpace.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>

int main (int argc, char** argv) {

  unsigned int MCevents = 1e7;
  double* data = new double[MCevents*6];
  TLorentzVector D(0.0, 0.0, 0.0, 1.8645);
  TLorentzVector d1,d2,d3,d4;
  Double_t m[4]={.13957018, .13957018, .493677, .13957018};
  TGenPhaseSpace event;
  event.SetDecay(D, 4, m);

  double weight;
  double m12, m34, cos12, cos34, phi;
  double m13, m14, m23, m24;
  double  wmax = 0;

  for (Int_t n=0;n<MCevents;n++) {
     weight = event.Generate();
     wmax = (wmax<weight)? weight : wmax;

     d1 = *(event.GetDecay(0));
     d2 = *(event.GetDecay(1));
     d3 = *(event.GetDecay(2));
     d4 = *(event.GetDecay(3));

     m12 = (d1+d2).M();
     m13 = (d1+d3).M();
     m14 = (d1+d4).M();
     m23 = (d2+d3).M();
     m24 = (d2+d4).M();
     m34 = (d3+d4).M();

     TVector3 d1n = d1.Vect().Unit();
     TVector3 d2n = d2.Vect().Unit();
     TVector3 d3n = d3.Vect().Unit();
     TVector3 d4n = d4.Vect().Unit();
     TVector3 d12n = (d1+d2).Vect().Unit();
     TVector3 d34n = (d3+d4).Vect().Unit();

     TVector3 n1 = d1n.Cross(d2n);
     TVector3 n2 = d3n.Cross(d4n);
     TVector3 n3 = n1.Unit().Cross(n2.Unit());

     // Calculation of the angle Phi between the planes.
     double cosp = n1.Unit().Dot(n2.Unit());
     double sinp = n3.Dot(d34n);
     phi = acos(cosp); 
     if (sinp <0) phi *= -1;

     //Vectors in rest fram of their resonance. 
     TLorentzVector d1r = d1;
     TLorentzVector d3r = d3;
     d1r.Boost(-(d1+d2).BoostVector());
     d3r.Boost(-(d3+d4).BoostVector());
     TVector3 d1rn = d1r.Vect().Unit();
     TVector3 d3rn = d3r.Vect().Unit();
     
     // helicity angle for d12 and d34 frame
     cos12 = d12n.Dot(d1rn);
     cos34 = d34n.Dot(d3rn);
     // printf("%.5g %.5g %.5g %.5g %.5g %.5g %.5g\n",m12, m13, m14, m23, m24, m34, weight);
     data[n*6] = m12;
     data[n*6+1] = m34;
     data[n*6+2] = cos12;
     data[n*6+3] = cos34;
     data[n*6+4] = phi;
     data[n*6+5] = weight;
  }

  double tm12, tm34, tcos12, tcos34, tphi, tweight;
  TFile f2("phspMC.root","recreate");
  TTree t1("t1","a simple Tree with simple variables");
  t1.Branch("m12",&tm12,"m12/D");
  t1.Branch("m34",&tm34,"m34/D");
  t1.Branch("cos12",&tcos12,"cos12/D");   
  t1.Branch("cos34",&tcos34,"cos34/D");
  t1.Branch("phi",&tphi,"phi/D");

  TRandom3 rnd(0);
  wmax*=1.3;
  for (int i = 0; i < MCevents ; ++i)
  { 
    if (rnd.Uniform(wmax)>=data[i*6+5]){
      continue;
    }
    tm12 = data[i*6];
    tm34 = data[i*6+1];
    tcos12 = data[i*6+2];
    tcos34 = data[i*6+3];
    tphi = data[i*6+4];
    t1.Fill();
  }
  t1.Write();

}