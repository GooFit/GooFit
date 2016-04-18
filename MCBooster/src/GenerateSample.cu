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
#include <sstream>
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

//ROOT
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

using namespace std;

using namespace MCBooster;


void splitString(const std::string &s, const char delim, std::vector<std::string> &elems)
{
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim))
	{
		elems.push_back(item);
	}
	//return elems;
	return;
}

void splitReal(const std::string &s, const char delim, std::vector<GReal_t> &elems)
{
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim))
	{
		elems.push_back(std::stod(item));
	}
	//return elems;
	return;
}



GInt_t main(int argv, char** argc)
{

	GULong_t nevents=0;
	string   output_dir="";
	GReal_t  mother;
	string _names;
	string _masses;

	try {

		TCLAP::CmdLine cmd("Command line arguments for GenerateSample", '=');

		TCLAP::ValueArg<GULong_t> eArg("n", "number-of-events",
				"Number of events",
				true, 1e6, "long");
		cmd.add(eArg);


		TCLAP::ValueArg<string> oArg("o", "output-file",
				"Output file",
				false, "./phsp.root", "string");
		cmd.add(oArg);


		TCLAP::ValueArg<string> pArg("p", "particles",
				"List of particles. First particle is the mother.Example: D0->Kpipi is 'D0;K;pi+;pi-",
				true, "", "string");
		cmd.add(pArg);


		TCLAP::ValueArg<string> mArg("m", "masses",
				"Particle mass.  First particle is the mother. Example: D0->Kpipi is '1.865;0.439;0.139;0.139",
				true, "" ,"string");
		cmd.add(mArg);

		// Parse the argv array.
		cmd.parse(argv, argc);

		// Get the value parsed by each arg.
		nevents        = eArg.getValue();
		output_dir     = oArg.getValue();
		_masses        = mArg.getValue();
		_names         = pArg.getValue();

	} catch (TCLAP::ArgException &e)  // catch any exceptions
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId()
							<< std::endl;
	}


	vector<string> names_temp;
	vector<GReal_t> masses_temp;

	splitReal(_masses, ';'  , masses_temp);
	splitString(_names, ';' , names_temp );

	if(masses_temp.size() < 3 || masses_temp.size()>9)
	{
		cout << "Exit. Number of particles is (< 2) or (> 9)." << endl;
		exit(0);
	}



	if( masses_temp.size() != names_temp.size() )
	{
		cout << "Exit. Number of particles is different of number of names." << endl;
		exit(0);
	}


	//dump configuration
	cout << "-----------------------------------------------------"<< endl;
	cout << "---------------------- MCBooster --------------------"<< endl;
	cout << "- Mother name: " << names_temp[0] << " mass: "<< masses_temp[0] << endl;
	for(GInt_t i=1;i<masses_temp.size();i++)
	{
		cout << "- Daughter name: " << names_temp[i] << " mass: "<< masses_temp[i] << endl;
	}
	cout << "- Number of events: "<<nevents<< endl;
	cout << "- Output file: "<< output_dir << endl;
	cout << "-----------------------------------------------------"<< endl;


	// Generation
	GReal_t mother_mass = masses_temp[0];
	vector<GReal_t> masses(masses_temp.size()-1);
	std::copy ( masses_temp.begin() +1, masses_temp.end(), masses.begin() );

	PhaseSpace phsp(mother_mass, masses, nevents);
	Events GenEvents(masses.size(), nevents);

	timespec time1, time2;
	phsp.Generate(Vector4R(mother_mass, 0.0, 0.0, 0.0));

	///-------------------------------------
	//unweight
	clock_gettime(CLOCK_REALTIME, &time1);
	phsp.Unweight();

	clock_gettime(CLOCK_REALTIME, &time2);

		GReal_t unweight_time_used = ((GReal_t) (time_diff(time1, time2).tv_sec
				+ time_diff(time1, time2).tv_nsec * 1.0e-9));

	//-------------------------------------
	clock_gettime(CLOCK_REALTIME, &time1);

	/// Create Events container
	phsp.Export(&GenEvents);

	clock_gettime(CLOCK_REALTIME, &time2);

	GReal_t exp_time_used = ((GReal_t) (time_diff(time1, time2).tv_sec
			+ time_diff(time1, time2).tv_nsec * 1.0e-9));

	phsp.Export(&GenEvents);

	cout << "-----------------------------------------------------"<< endl;
	cout << "----------------------- Timing ----------------------"<< endl;
	cout << "Event generation: " << phsp.GetEvtTime() << endl;
	cout << "Unweight generation: " << unweight_time_used << endl;
	cout << "Export events to host: " << exp_time_used << endl;
	cout << "-----------------------------------------------------"<< endl;

	TFile *file = new TFile( output_dir.c_str() , "RECREATE");
	TTree *tree = new TTree("events", "events");

	TLorentzVector* decayVectors = new TLorentzVector[masses.size()+1 ];

	for (GInt_t p = 0; p < masses.size() ; p++) {

		tree->Branch(names_temp[p+1].c_str(), names_temp[p+1].c_str(), &decayVectors[p]);
	}

	GReal_t wevt, wmax;
	GInt_t flag;

	tree->Branch("weightEvt",  &wevt, "weightEvt/D");
	tree->Branch("weightMax",  &wmax, "weightMax/D");
	tree->Branch("AccRej"   ,  &flag, "AccRej/I");

	cout << "\n Storing events in Tree..."<<endl;

	for (GInt_t evt = 0; evt < nevents; evt++) {
		for (GInt_t p = 0; p < masses.size(); p++) {

			decayVectors[p].SetPxPyPzE(
					GenEvents.fDaughters[p][evt].get(1),
					GenEvents.fDaughters[p][evt].get(2),
					GenEvents.fDaughters[p][evt].get(3),
					GenEvents.fDaughters[p][evt].get(0)
					);

		}

		wevt = GenEvents.fWeights[evt];
		wmax = GenEvents.fMaxWeight;
        flag = GenEvents.fAccRejFlags[evt];

		tree->Fill();

	}
	cout << "Done. \n"<<endl;
	tree->Write();
	file->Close();

	return 0;
}
