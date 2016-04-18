/*
 * Generate.cu
 *
 * Copyright 2016 Antonio Augusto Alves Junior
 *
 * Created on : 10/03/2016
 *      Author: augalves
 */

/*
    This file is part of MCBooster.

    MCBooster is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MCBooster is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MCBooster.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include <time.h>
#include <string>
#include <map>

#include <math.h>

#define CUDA_API_PER_THREAD_DEFAULT_STREAM

//this lib
#include <mcbooster/GTypes.h>
#include <mcbooster/Vector4R.h>
#include <mcbooster/Generate.h>
#include <mcbooster/GContainers.h>
#include <mcbooster/Evaluate.h>
#include <mcbooster/EvaluateArray.h>
#include <mcbooster/GFunctional.h>
//root
#include <TROOT.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TApplication.h>
#include <TCanvas.h>

using namespace std;

using namespace MCBooster;



struct Dataset: public IFunctionArray
{
	Dataset()
	{
		dim = 4;
	}

	__host__   __device__ GReal_t cosHELANG(const Vector4R p, const Vector4R q,
			const Vector4R d)
	{
		GReal_t pd = p * d;
		GReal_t pq = p * q;
		GReal_t qd = q * d;
		GReal_t mp2 = p.mass2();
		GReal_t mq2 = q.mass2();
		GReal_t md2 = d.mass2();

		return (pd * mq2 - pq * qd)
				/ sqrt((pq * pq - mq2 * mp2) * (qd * qd - mq2 * md2));

	}

	__host__   __device__ GReal_t deltaAngle(const Vector4R& p4_p,
			const Vector4R& p4_d1, const Vector4R& p4_d2, const Vector4R& p4_h1,
			const Vector4R& p4_h2)
	{

		Vector4R p4_d1p, p4_h1p, p4_h2p, p4_d2p;

		Vector4R d1_perp, d1_prime, h1_perp;
		Vector4R D;

		D = p4_d1 + p4_d2;

		d1_perp = p4_d1 - (D.dot(p4_d1) / D.dot(D)) * D;
		h1_perp = p4_h1 - (D.dot(p4_h1) / D.dot(D)) * D;

		// orthogonal to both D and d1_perp

		d1_prime = D.cross(d1_perp);

		d1_perp = d1_perp / d1_perp.d3mag();
		d1_prime = d1_prime / d1_prime.d3mag();

		GReal_t x, y;

		x = d1_perp.dot(h1_perp);
		y = d1_prime.dot(h1_perp);

		GReal_t chi = atan2(y, x);

		if (chi < 0.0)
			chi += 2.0*CUDART_PI_HI;

		return chi;

	}

	__host__ __device__
	void operator()(const GInt_t n, Vector4R** particles, GReal_t* variables)
	{
		Vector4R pJpsi = *particles[0];
		Vector4R pK    = *particles[1];
		Vector4R ppi   = *particles[2];
		Vector4R pMup  = *particles[3];
		Vector4R pMum  = *particles[4];

		//K* helicity angle
		Vector4R pB0 = pJpsi + pK + ppi;
		Vector4R pKpi = pK + ppi;
		Vector4R pJpsipi = pJpsi + ppi;


		variables[0] = pKpi.mass();
		variables[1] = pJpsipi.mass();
		variables[2] = cosHELANG(pB0, pKpi, pK);
		variables[3] = cosHELANG(pB0, pJpsi, pMup);
		variables[4] = deltaAngle(pB0, pK, ppi, pMup, pMum);

	}

};



GInt_t main(void)
{

	TApplication *myapp=new TApplication("myapp",0,0);

	GLong_t events = 10000000;
	size_t  ndaughters = 3;
	GReal_t mass0 = 5.2795;

	//Particles mothers(events, Vector4R(5.2795,0.0,0.0,0.0) );

	//{ "J/psi", "K", "pi", "pi" }

	vector<string> namesB0;
	namesB0.push_back("J/#psi");
	namesB0.push_back("K");
	namesB0.push_back("pi");

	//{ 3.096916, 0.493677, 0.13957018 }
	vector<GReal_t> massesB0;
	massesB0.push_back(3.096916);
	massesB0.push_back(0.493677);
	massesB0.push_back(0.13957018);

	//chronometer
	timespec time1, time2;

	cout << "=========================================================" << endl;
	cout << "===================   B0 -> J/psi K pi=  ================" << endl;
	cout << "=========================================================" << endl;
	cout << "Number of events: " << events << endl;

	clock_gettime(CLOCK_REALTIME, &time1);

	/// Create PhaseSpace object for B0-> K pi J/psi
	PhaseSpace phsp(5.2795, massesB0, events);

	clock_gettime(CLOCK_REALTIME, &time2);

	GReal_t phsp_time_used = ((GReal_t) (time_diff(time1, time2).tv_sec
			+ time_diff(time1, time2).tv_nsec * 1.0e-9));

	cout << "|\t PhaseSpace ctor time [B0]:\t " << phsp_time_used << " s"
			<< endl;

	clock_gettime(CLOCK_REALTIME, &time1);

	/// Generate events B0-> K pi J/psi
	phsp.Generate(Vector4R(5.2795, 0.0, 0.0, 0.0));

	clock_gettime(CLOCK_REALTIME, &time2);

	phsp.Unweight();

	GReal_t cpu_time_used;
	cpu_time_used = ((GReal_t) (time_diff(time1, time2).tv_sec
			+ time_diff(time1, time2).tv_nsec * 1.0e-9));

	cout << "|\t Generate time [B0]:\t " << cpu_time_used << " s" << endl;

	clock_gettime(CLOCK_REALTIME, &time1);

	/// Create Events container
	Events *MyEvents = new Events(ndaughters, events);

	clock_gettime(CLOCK_REALTIME, &time2);

	GReal_t evt_time_used = ((GReal_t) (time_diff(time1, time2).tv_sec
			+ time_diff(time1, time2).tv_nsec * 1.0e-9));

	cout << "|\t Event ctor time [B0]:\t " << evt_time_used << " s" << endl;

	clock_gettime(CLOCK_REALTIME, &time1);

	/// Export events
	phsp.Export(MyEvents);

	clock_gettime(CLOCK_REALTIME, &time2);

	GReal_t exp_time_used = ((GReal_t) (time_diff(time1, time2).tv_sec
			+ time_diff(time1, time2).tv_nsec * 1.0e-9));

	cout << "|\t Export time [B0]:\t\t " << exp_time_used << " s" << endl;

	cout << "=========================================================" << endl;

	/// Print the first 10 events with corresponding weights
	for (GInt_t event = 0; event < 10; event++)
	{
		cout << "Event: " << event << endl << "\t| weight "
				<< MyEvents->fWeights[event]
				<< "\t| flag " << MyEvents->fAccRejFlags[event] <<endl;

		for (GInt_t daughter = 0; daughter < 3; daughter++)
		{
			cout << " \t| " << namesB0[daughter] << " : mass "
					<< MyEvents->fDaughters[daughter][event].mass()
					<< " 4-momentum ( "
					<< MyEvents->fDaughters[daughter][event].get(0) << ", "
					<< MyEvents->fDaughters[daughter][event].get(1) << ", "
					<< MyEvents->fDaughters[daughter][event].get(2) << ", "
					<< MyEvents->fDaughters[daughter][event].get(3) << ")  "
					<< endl;

		}
		cout << endl;

	}


	//{ "J/psi", "K", "pi", "pi" }

	vector<string> namesJpsi;
	namesJpsi.push_back("mu+");
	namesJpsi.push_back("mu-");

	//{ 3.096916, 0.493677, 0.13957018 }
	vector<GReal_t> massesJpsi;
	massesJpsi.push_back(0.100);
	massesJpsi.push_back(0.100);

	cout << "=========================================================" << endl;
	cout << "====================   J/psi -> mu mu  ==================" << endl;
	cout << "=========================================================" << endl;
    cout << "Number of events: " << events << endl;

	clock_gettime(CLOCK_REALTIME, &time1);
	Events *MyEventsJpsi = new Events(2, events);
	clock_gettime(CLOCK_REALTIME, &time2);
	GReal_t evt_time_usedJpsi = ((GReal_t) (time_diff(time1, time2).tv_sec
			+ time_diff(time1, time2).tv_nsec * 1.0e-9));

	cout << "|\t Event ctor time [J/psi]:\t " << evt_time_usedJpsi << " s"
			<< endl;

	//Decays trees(mothers, names, masses );

	clock_gettime(CLOCK_REALTIME, &time1);
	PhaseSpace phspJpsi(3.096916, massesJpsi, events);
	clock_gettime(CLOCK_REALTIME, &time2);
	GReal_t phsp_time_usedJpsi = ((GReal_t) (time_diff(time1, time2).tv_sec
			+ time_diff(time1, time2).tv_nsec * 1.0e-9));

	cout << "|\t PhaseSpace ctor time [J/psi]:\t " << phsp_time_usedJpsi << " s"
			<< endl;

	clock_gettime(CLOCK_REALTIME, &time1);
	phspJpsi.Generate(phsp.GetDaughters(0));
	clock_gettime(CLOCK_REALTIME, &time2);

	phspJpsi.Unweight();

	GReal_t cpu_time_usedJpsi;
	cpu_time_usedJpsi = ((GReal_t) (time_diff(time1, time2).tv_sec
			+ time_diff(time1, time2).tv_nsec * 1.0e-9));

	cout << "|\t Generate time [J/psi]:\t " << cpu_time_usedJpsi << " s"
			<< endl;

	clock_gettime(CLOCK_REALTIME, &time1);
	phspJpsi.Export(MyEventsJpsi);
	clock_gettime(CLOCK_REALTIME, &time2);

	GReal_t exp_time_usedJpsi = ((GReal_t) (time_diff(time1, time2).tv_sec
			+ time_diff(time1, time2).tv_nsec * 1.0e-9));

	cout << "|\t Export time [J/psi]:\t\t " << exp_time_usedJpsi << " s"
			<< endl;

	cout << "=========================================================" << endl;

	for (GInt_t event = 0; event < 10; event++)
	{
		cout << "Event: " << event << endl << "\t| weight "
				<< MyEventsJpsi->fWeights[event]
				<< "\t| flag " << MyEventsJpsi->fAccRejFlags[event] <<endl;

		for (GInt_t daughter = 0; daughter < 2; daughter++)
		{
			cout << " \t| " << namesJpsi[daughter] << " : mass "
					<< MyEventsJpsi->fDaughters[daughter][event].mass()
					<< " 4-momentum ( "
					<< MyEventsJpsi->fDaughters[daughter][event].get(0) << ", "
					<< MyEventsJpsi->fDaughters[daughter][event].get(1) << ", "
					<< MyEventsJpsi->fDaughters[daughter][event].get(2) << ", "
					<< MyEventsJpsi->fDaughters[daughter][event].get(3) << ")  "
					<< endl;

		}
		cout << endl;

	}

	VariableSet_h Var(5);
	RealVector_h result_MKpi(events);
	RealVector_h result_MJpsipi(events);
	RealVector_h result_CosThetaK(events);
	RealVector_h result_CosThetaMu(events);
	RealVector_h result_DeltaAngle(events);

	Var[0] = &result_MKpi;
	Var[1] = &result_MJpsipi;
	Var[2] = &result_CosThetaK;
	Var[3] = &result_CosThetaMu;
	Var[4] = &result_DeltaAngle;

	ParticlesSet_d JpsiKpiMuMu(5);
	JpsiKpiMuMu[0] = &phsp.GetDaughters(0);
	JpsiKpiMuMu[1] = &phsp.GetDaughters(1);
	JpsiKpiMuMu[2] = &phsp.GetDaughters(2);
	JpsiKpiMuMu[3] = &phspJpsi.GetDaughters(0);
	JpsiKpiMuMu[4] = &phspJpsi.GetDaughters(1);

	Dataset DataJpsiKpi = Dataset();

	clock_gettime(CLOCK_REALTIME, &time1);
	EvaluateArray<Dataset>(DataJpsiKpi, JpsiKpiMuMu, Var);
	clock_gettime(CLOCK_REALTIME, &time2);
	GReal_t Dataset_time = ((GReal_t) (time_diff(time1, time2).tv_sec
			+ time_diff(time1, time2).tv_nsec * 1.0e-9));
	cout << "=========================================================" << endl;
	cout << "=================   Evaluate Dataset   ===============" << endl;
	cout << "=========================================================" << endl;
	cout << "|\t Dataset time : " << Dataset_time << " s" << endl;
	cout << "=========================================================" << endl;


	 for(GInt_t event=0; event<10; event++ )
	 {

	 cout
	 << event
	 <<" " <<  result_MKpi[event]
	 <<" " <<  result_MJpsipi[event]
	 <<" " <<  result_CosThetaK[event]
	 <<" " <<  result_CosThetaMu[event]
	 <<" " <<  result_DeltaAngle[event]
	 << endl;
	 }


	 TH1D *cosThetaK  = new TH1D( "cosThetaK", ";Cos(#theta_{K});Events", 100, -1.0, 1.0);
	 TH1D *cosThetaMu = new TH1D("cosThetaMu", ";Cos(#theta_{#Mu});Events", 100, -1.0, 1.0);
	 TH1D *deltaAngle = new TH1D("deltaAngle", ";#Delta #phi;Events", 100, 0.0, 6.3);
	 TH1D *MKpi       = new TH1D("MKpi", ";M(K,#pi);Events", 100, massesB0[1]+massesB0[2], mass0 - massesB0[0]  );
	 TH1D *MJpsipi    = new TH1D("MJpsipi", ";M(Jpsi/#psi,#pi);Events", 100, massesB0[0]+massesB0[2], mass0 - massesB0[1]  );


	 for(GInt_t event=0; event<events; event++ )
	 {

	 cosThetaK->Fill( result_CosThetaK[event], MyEvents->fWeights[event]  );
	 cosThetaMu->Fill(result_CosThetaMu[event], MyEvents->fWeights[event]  );
	 deltaAngle->Fill(result_DeltaAngle[event], MyEvents->fWeights[event]  );
	 MKpi->Fill(result_MKpi[event], MyEvents->fWeights[event]  );
	 MJpsipi->Fill(result_MJpsipi[event], MyEvents->fWeights[event]  );

	 }

	 TCanvas *c_cosThetaK = new TCanvas( "cosThetaK", "", 600, 500 );
	 cosThetaK->Draw("HIST");
	 cosThetaK->SetLineColor(4);
	 cosThetaK->SetLineWidth(2);
	 cosThetaK->SetStats(0);
	 cosThetaK->SetMinimum(0);
	 c_cosThetaK->Print( "cosThetaK.pdf" );

	 TCanvas *c_cosThetaMu = new TCanvas( "cosThetaMu", "", 600, 500 );
	 cosThetaMu->Draw("HIST");
	 cosThetaMu->SetLineColor(4);
	 cosThetaMu->SetLineWidth(2);
	 cosThetaMu->SetStats(0);
	 cosThetaMu->SetMinimum(0);
	 c_cosThetaMu->Print( "cosThetaMu.pdf" );

	 TCanvas *c_deltaAngle = new TCanvas( "deltaAngle", "", 600, 500 );
	 deltaAngle->Draw("HIST");
	 deltaAngle->SetLineColor(4);
	 deltaAngle->SetLineWidth(2);
	 deltaAngle->SetStats(0);
	 deltaAngle->SetMinimum(0);
	 c_deltaAngle->Print( "deltaAngle.pdf" );

	 TCanvas *c_MKpi = new TCanvas( "MKpi", "", 600, 500 );
	 MKpi->Draw("HIST");
	 MKpi->SetLineColor(4);
	 MKpi->SetLineWidth(2);
	 MKpi->SetStats(0);
	 MKpi->SetMinimum(0);
	 c_MKpi->Print( "MKpi.pdf" );

	 TCanvas *c_MJpsipi = new TCanvas( "MJpsipi", "", 600, 500 );
	 MJpsipi->Draw("HIST");
	 MJpsipi->SetLineColor(4);
	 MJpsipi->SetLineWidth(2);
	 MJpsipi->SetStats(0);
	 MJpsipi->SetMinimum(0);
	 c_MJpsipi->Print( "MJpsipi.pdf" );




	 TH2D *dalitz = new TH2D("dalitz",
			 TString::Format("Weigted;M^{2}(%s,%s) [GeV^{2}/c^{4}]; M^{2}(%s,%s) [GeV^{2}/c^{4}]"
					 , namesB0[0].c_str(), namesB0[1].c_str()
					 , namesB0[1].c_str(), namesB0[2].c_str() ).Data(),
					 100, pow(massesB0[0]+massesB0[2],2), pow(mass0 - massesB0[1],2),
					 100, pow(massesB0[1]+massesB0[2],2), pow(mass0 - massesB0[0],2) );


	 for(GInt_t event=0; event<events; event++ )
	 {

		 dalitz->Fill( result_MJpsipi[event]*result_MJpsipi[event],
				 result_MKpi[event]*result_MKpi[event],
				 MyEvents->fWeights[event]  );

	 }

	 TCanvas *c2 = new TCanvas( "dalitz", "", 600, 500 );
	 dalitz->Draw("COLZ");
	 dalitz->SetStats(0);
	 c2->Print( "dalitzW.pdf" );

	 TH2D *dalitz2 = new TH2D("dalitz2",
			 TString::Format("Unweigted;M^{2}(%s,%s) [GeV^{2}/c^{4}]; M^{2}(%s,%s) [GeV^{2}/c^{4}]"
					 , namesB0[0].c_str(), namesB0[1].c_str()
					 , namesB0[1].c_str(), namesB0[2].c_str() ).Data(),
					 100, pow(massesB0[0]+massesB0[2],2), pow(mass0 - massesB0[1],2),
					 100, pow(massesB0[1]+massesB0[2],2), pow(mass0 - massesB0[0],2) );


	 for(GInt_t event=0; event<events; event++ )
	 {

		 dalitz2->Fill( result_MJpsipi[event]*result_MJpsipi[event],
				 result_MKpi[event]*result_MKpi[event],
				 MyEvents->fAccRejFlags[event]  );

	 }


	 TCanvas *c3 = new TCanvas( "dalitz2", "", 600, 500 );
		 dalitz2->Draw("COLZ");
		 dalitz2->SetStats(0);
		 c2->Print( "dalitzU.pdf" );



	 myapp->Run();
	 return 0;

}
