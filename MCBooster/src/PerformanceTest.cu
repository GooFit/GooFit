/*
 * PerformanceTest.cu
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

//command line
#include <tclap/CmdLine.h>
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
#include "TFile.h"
#include "TString.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TGraph.h"

using namespace std;

using namespace MCBooster;

void RunMCGen (GInt_t nfinal, GInt_t nevents,  Double_t *time )
{


	vector<string> names;
	vector<GReal_t> masses;

	for(GInt_t i=0; i<nfinal; i++ )
	{
		names.push_back(TString::Format("f%d",i).Data());
		masses.push_back(0.13957018);
	}

	PhaseSpace phsp(5.2795, masses, nevents);
	Events _Events(nfinal, nevents);
    timespec time1, time2;


	clock_gettime(CLOCK_REALTIME, &time1);

	phsp.Generate(Vector4R(5.2795, 0.0, 0.0, 0.0));

	clock_gettime(CLOCK_REALTIME, &time2);

	GReal_t cpu_time_used;
	cpu_time_used = ((GReal_t)( time_diff (time1,time2).tv_sec + time_diff(time1,time2).tv_nsec*1.0e-9));

	clock_gettime(CLOCK_REALTIME, &time1);

	/// Create Events container
	phsp.Export(&_Events);

	clock_gettime(CLOCK_REALTIME, &time2);

	GReal_t exp_time_used = ((GReal_t) (time_diff(time1, time2).tv_sec
				+ time_diff(time1, time2).tv_nsec * 1.0e-9));


	time[0] =  phsp.GetRndTime();
	time[1] =  phsp.GetEvtTime();
	time[2] =  exp_time_used;
	time[3] =  cpu_time_used;


}


GInt_t main(int argv, char** argc)
{

	GULong_t nevents_scan=0;
	string suffix="";
	string output_dir="";

		try {

			TCLAP::CmdLine cmd("Command line arguments for PerformanceTest", ' ', "0.9");

			TCLAP::ValueArg<GULong_t> eArg("e", "events", "Number of events in time profile",
					false, 1e6, "long integer");
			cmd.add(eArg);

			TCLAP::ValueArg<string> sArg("s", "suffix",
					"Suffix for plots name",
					false, "", "string");
			cmd.add(sArg);

			TCLAP::ValueArg<string> oArg("o", "output_dir", "output directory",
					false, ".", "string");
			cmd.add(oArg);

			// Parse the argv array.
			cmd.parse(argv, argc);

			// Get the value parsed by each arg.
			nevents_scan   = eArg.getValue();
			suffix         = sArg.getValue();
			output_dir     = oArg.getValue();

		} catch (TCLAP::ArgException &e)  // catch any exceptions
		{
			std::cerr << "error: " << e.error() << " for arg " << e.argId()
					<< std::endl;
		}

	vector<GDouble_t>     _NParticles;
	vector<GDouble_t>      _NEvents;

	vector<GDouble_t>    _Rnd_Time1;
	vector<GDouble_t>    _Evt_Time1;
	vector<GDouble_t>    _Out_Time1;
	vector<GDouble_t>    _Cpu_Time1;

	vector<GDouble_t>    _Rnd_Time2;
	vector<GDouble_t>    _Evt_Time2;
	vector<GDouble_t>    _Out_Time2;
	vector<GDouble_t>    _Cpu_Time2;




	GDouble_t _time[4];

	for(GInt_t i=1;i<100;i++)
	{
		GULong_t nevt = i*500000;

		RunMCGen(3, nevt , &_time[0] );

		 _NEvents.push_back(     nevt );
		_Rnd_Time1.push_back( _time[0] );
		_Evt_Time1.push_back( _time[1] );
		_Out_Time1.push_back( _time[2] );
		_Cpu_Time1.push_back( _time[3] );

	}

	for(GInt_t i=2;i<10;i++)
		{


			RunMCGen(i, nevents_scan , &_time[0] );
			_NParticles.push_back( i );
			_Rnd_Time2.push_back( _time[0] );
			_Evt_Time2.push_back( _time[1] );
			_Out_Time2.push_back( _time[2] );
			_Cpu_Time2.push_back( _time[3] );

		}




	TApplication *myapp=new TApplication("myapp",0,0);

	TGraph  *Rnd_Time1 = new TGraph( _NEvents.size(), _NEvents.data(), _Rnd_Time1.data());
	TGraph  *Evt_Time1 = new TGraph( _NEvents.size(), _NEvents.data(), _Evt_Time1.data());
	TGraph  *Out_Time1 = new TGraph( _NEvents.size(), _NEvents.data(), _Out_Time1.data());
	TGraph  *Cpu_Time1 = new TGraph( _NEvents.size(), _NEvents.data(), _Cpu_Time1.data());

	TCanvas *c1 = new TCanvas( "Rnd_Time1", "Rnd_Time1", 600, 500 );
	Rnd_Time1->SetLineColor(2);
	Rnd_Time1->SetLineWidth(4);
	Rnd_Time1->SetMarkerColor(4);
	Rnd_Time1->SetMarkerSize(1.0);
	Rnd_Time1->SetMarkerStyle(21);
	Rnd_Time1->SetTitle("Randon number generation time");
	Rnd_Time1->GetXaxis()->SetTitle("Number of events");
	Rnd_Time1->GetYaxis()->SetTitle("Time (s)");
	Rnd_Time1->GetYaxis()->SetTitleOffset(1.4);
	Rnd_Time1->Draw("ACP");

	TCanvas *c2 = new TCanvas( "Evt_Time1", "Evt_Time1", 600, 500 );
	Evt_Time1->SetLineColor(2);
	Evt_Time1->SetLineWidth(4);
	Evt_Time1->SetMarkerColor(4);
	Evt_Time1->SetMarkerSize(1.0);
	Evt_Time1->SetMarkerStyle(21);
	Evt_Time1->SetTitle("Events generation time");
	Evt_Time1->GetXaxis()->SetTitle("Number of events");
	Evt_Time1->GetYaxis()->SetTitle("Time (s)");
	Evt_Time1->GetYaxis()->SetTitleOffset(1.4);
	Evt_Time1->Draw("ACP");


	TCanvas *c3 = new TCanvas( "Out_Time1", "Out_Time1", 600, 500 );
	Out_Time1->SetLineColor(2);
	Out_Time1->SetLineWidth(4);
	Out_Time1->SetMarkerColor(4);
	Out_Time1->SetMarkerSize(1.0);
	Out_Time1->SetMarkerStyle(21);
	Out_Time1->SetTitle("Output time");
	Out_Time1->GetXaxis()->SetTitle("Number of events");
	Out_Time1->GetYaxis()->SetTitle("Time (s)");
	Out_Time1->GetYaxis()->SetTitleOffset(1.4);
	Out_Time1->Draw("ACP");

	TCanvas *c4 = new TCanvas( "Cpu_Time", "Cpu_Time", 600, 500 );
	Cpu_Time1->SetLineColor(2);
	Cpu_Time1->SetLineWidth(4);
	Cpu_Time1->SetMarkerColor(4);
	Cpu_Time1->SetMarkerSize(1.0);
	Cpu_Time1->SetMarkerStyle(21);
	Cpu_Time1->SetTitle("Generation time");
	Cpu_Time1->GetXaxis()->SetTitle("Number of events");
	Cpu_Time1->GetYaxis()->SetTitle("Time (s)");
	Cpu_Time1->GetYaxis()->SetTitleOffset(1.4);
	Cpu_Time1->Draw("ACP");

	c1->SaveAs(TString::Format("%s/Rnd_Time_%s.pdf", output_dir.c_str(), suffix.c_str()));
	c2->SaveAs(TString::Format("%s/Evt_Time_%s.pdf", output_dir.c_str(), suffix.c_str()));
	c3->SaveAs(TString::Format("%s/Out_Time_%s.pdf", output_dir.c_str(), suffix.c_str()));
	c4->SaveAs(TString::Format("%s/Cpu_Time_%s.pdf", output_dir.c_str(), suffix.c_str()));



	TGraph  *Rnd_Time2 = new TGraph(_NParticles.size(), _NParticles.data(), _Rnd_Time2.data());
	TGraph  *Evt_Time2 = new TGraph(_NParticles.size(), _NParticles.data(), _Evt_Time2.data());
	TGraph  *Out_Time2 = new TGraph(_NParticles.size(), _NParticles.data(), _Out_Time2.data());
	TGraph  *Cpu_Time2 = new TGraph(_NParticles.size(), _NParticles.data(), _Cpu_Time2.data());

	TCanvas *d1 = new TCanvas( "Rnd_Time2", "Rnd_Time2", 600, 500 );
	Rnd_Time2->SetLineColor(2);
	Rnd_Time2->SetLineWidth(4);
	Rnd_Time2->SetMarkerColor(4);
	Rnd_Time2->SetMarkerSize(1.0);
	Rnd_Time2->SetMarkerStyle(21);
	Rnd_Time2->SetTitle("Randon number generation time");
	Rnd_Time2->GetXaxis()->SetTitle("Number of particles");
	Rnd_Time2->GetYaxis()->SetTitle("Time (s)");
	Rnd_Time2->GetYaxis()->SetTitleOffset(1.4);
	Rnd_Time2->Draw("ACP");

	TCanvas *d2 = new TCanvas( "Evt_Time2", "Evt_Time2", 600, 500 );
	Evt_Time2->SetLineColor(2);
	Evt_Time2->SetLineWidth(4);
	Evt_Time2->SetMarkerColor(4);
	Evt_Time2->SetMarkerSize(1.0);
	Evt_Time2->SetMarkerStyle(21);
	Evt_Time2->SetTitle("Events generation time");
	Evt_Time2->GetXaxis()->SetTitle("Number of particles");
	Evt_Time2->GetYaxis()->SetTitle("Time (s)");
	Evt_Time2->GetYaxis()->SetTitleOffset(1.4);
	Evt_Time2->Draw("ACP");


	TCanvas *d3 = new TCanvas( "Out_Time2", "Out_Time2", 600, 500 );
	Out_Time2->SetLineColor(2);
	Out_Time2->SetLineWidth(4);
	Out_Time2->SetMarkerColor(4);
	Out_Time2->SetMarkerSize(1.0);
	Out_Time2->SetMarkerStyle(21);
	Out_Time2->SetTitle("Output time");
	Out_Time2->GetXaxis()->SetTitle("Number of particles");
	Out_Time2->GetYaxis()->SetTitle("Time (s)");
	Out_Time2->GetYaxis()->SetTitleOffset(1.4);
	Out_Time2->Draw("ACP");

	TCanvas *d4 = new TCanvas( "Cpu_Time2", "Cpu_Time", 600, 500 );
	Cpu_Time2->SetLineColor(2);
	Cpu_Time2->SetLineWidth(4);
	Cpu_Time2->SetMarkerColor(4);
	Cpu_Time2->SetMarkerSize(1.0);
	Cpu_Time2->SetMarkerStyle(21);
	Cpu_Time2->SetTitle("Generation time");
	Cpu_Time2->GetXaxis()->SetTitle("Number of particles");
	Cpu_Time2->GetYaxis()->SetTitle("Time (s)");
	Cpu_Time2->GetYaxis()->SetTitleOffset(1.4);
	Cpu_Time2->Draw("ACP");

	d1->SaveAs(TString::Format("%s/Rnd_Time2_%s.pdf", output_dir.c_str(), suffix.c_str()));
	d2->SaveAs(TString::Format("%s/Evt_Time2_%s.pdf", output_dir.c_str(), suffix.c_str()));
	d3->SaveAs(TString::Format("%s/Out_Time2_%s.pdf", output_dir.c_str(), suffix.c_str()));
	d4->SaveAs(TString::Format("%s/Cpu_Time2_%s.pdf", output_dir.c_str(), suffix.c_str()));


	myapp->Run();
		return 0;
}
