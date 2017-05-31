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
#include <time.h>
#include <string>
#include <sstream>
#include <map>
#include <vector>
// command line
#include <tclap/CmdLine.h>
// this lib
#include <mcbooster/GTypes.h>
#include <mcbooster/Vector4R.h>
#include <mcbooster/Generate.h>
#include <mcbooster/GContainers.h>
#include <mcbooster/Evaluate.h>
#include <mcbooster/EvaluateArray.h>
#include <mcbooster/GFunctional.h>

// ROOT
#include "TROOT.h"
#include "TString.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TString.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TGraph.h"
#include "TLegend.h"

using namespace std;

using namespace mcbooster;

GInt_t factorial(GInt_t n) { return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n; }

void splitString(const std::string &s, const char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    // return elems;
    return;
}

void splitReal(const std::string &s, const char delim, std::vector<GReal_t> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(std::stod(item));
    }
    // return elems;
    return;
}

GInt_t main(int argv, char **argc) {
    GULong_t nevents  = 0;
    string output_dir = "";
    GReal_t mother;
    string _names;
    string _masses;

    try {
        TCLAP::CmdLine cmd("Command line arguments for GenerateSample", '=');

        TCLAP::ValueArg<GULong_t> eArg("n", "number-of-events", "Number of events", true, 1e6, "long");
        cmd.add(eArg);

        TCLAP::ValueArg<std::string> oArg("o", "output-file", "Output file", false, "./phsp.root", "string");
        cmd.add(oArg);

        TCLAP::ValueArg<std::string> pArg(
            "p",
            "particles",
            "List of particles. First particle is the mother.Example: D0->Kpipi is 'D0;K;pi+;pi-",
            true,
            "",
            "string");
        cmd.add(pArg);

        TCLAP::ValueArg<std::string> mArg(
            "m",
            "masses",
            "Particle mass.  First particle is the mother. Example: D0->Kpipi is '1.865;0.439;0.139;0.139",
            true,
            "",
            "string");
        cmd.add(mArg);

        // Parse the argv array.
        cmd.parse(argv, argc);

        // Get the value parsed by each arg.
        nevents    = eArg.getValue();
        output_dir = oArg.getValue();
        _masses    = mArg.getValue();
        _names     = pArg.getValue();

    } catch(TCLAP::ArgException &e) // catch any exceptions
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }

    vector<std::string> names_temp;
    vector<GReal_t> masses_temp;

    splitReal(_masses, ';', masses_temp);
    splitString(_names, ';', names_temp);

    if(masses_temp.size() < 3 || masses_temp.size() > 9) {
        cout << "Exit. Number of particles is (< 2) or (> 9)." << std::endl;
        exit(0);
    }

    if(masses_temp.size() != names_temp.size()) {
        cout << "Exit. Number of particles is different of number of names." << std::endl;
        exit(0);
    }

    // dump configuration
    cout << "-----------------------------------------------------" << std::endl;
    cout << "---------------------- MCBooster --------------------" << std::endl;
    cout << "- Mother name: " << names_temp[0] << " mass: " << masses_temp[0] << std::endl;
    for(GInt_t i = 1; i < masses_temp.size(); i++) {
        cout << "- Daughter name: " << names_temp[i] << " mass: " << masses_temp[i] << std::endl;
    }
    cout << "- Number of events: " << nevents << std::endl;
    cout << "- Output file: " << output_dir << std::endl;
    cout << "-----------------------------------------------------" << std::endl;

    // Generation
    GReal_t mother_mass = masses_temp[0];
    vector<GReal_t> masses(masses_temp.size() - 1);
    std::copy(masses_temp.begin() + 1, masses_temp.end(), masses.begin());

    PhaseSpace phsp(mother_mass, masses, nevents);
    Events GenEvents(masses.size(), nevents);

    timespec time1, time2;
    phsp.Generate(Vector4R(mother_mass, 0.0, 0.0, 0.0));

    ///-------------------------------------
    // unweight
    clock_gettime(CLOCK_REALTIME, &time1);
    phsp.Unweight();

    clock_gettime(CLOCK_REALTIME, &time2);

    GReal_t unweight_time_used = ((GReal_t)(time_diff(time1, time2).tv_sec + time_diff(time1, time2).tv_nsec * 1.0e-9));

    //-------------------------------------
    clock_gettime(CLOCK_REALTIME, &time1);

    /// Create Events container
    phsp.Export(&GenEvents);

    clock_gettime(CLOCK_REALTIME, &time2);

    GReal_t exp_time_used = ((GReal_t)(time_diff(time1, time2).tv_sec + time_diff(time1, time2).tv_nsec * 1.0e-9));

    phsp.Export(&GenEvents);

    cout << "-----------------------------------------------------" << std::endl;
    cout << "----------------------- Timing ----------------------" << std::endl;
    cout << "Event generation: " << phsp.GetEvtTime() << std::endl;
    cout << "Unweight generation: " << unweight_time_used << std::endl;
    cout << "Export events to host: " << exp_time_used << std::endl;
    cout << "-----------------------------------------------------" << std::endl;

    TApplication *myapp = new TApplication("myapp", 0, 0);

    GInt_t number_of_combinations = factorial(masses.size()) / (2 * factorial(masses.size() - 2));

    map<GInt_t, TH1D *> H1D;
    map<GInt_t, TH1D *> H1D_ROOT;

    map<GInt_t, GReal_t> min_histo_limits;
    map<GInt_t, GReal_t> max_histo_limits;

    for(GInt_t i = 0; i < masses.size(); i++) {
        for(GInt_t j = 0; j < masses.size(); j++) {
            if(j >= i)
                continue;

            GInt_t index = i + j * masses.size();
            TString name = TString::Format("M(%s,%s)", names_temp[j + 1].c_str(), names_temp[i + 1].c_str());
            GReal_t min  = masses[j] + masses[i];
            GReal_t max  = mother_mass;
            for(GInt_t k = 0; k < masses.size(); k++) {
                if((k != j) && (k != i))
                    max -= masses[k];
            }
            max_histo_limits[index] = max;
            min_histo_limits[index] = min;
            cout << "Adding [" << index << "] ---> M( " << names_temp[j + 1] << ", " << names_temp[i + 1] << ")"
                 << std::endl;
            cout << "Histogram " << name.Data() << " min " << min << " max " << max << std::endl;
            TH1D *h
                = new TH1D(name.Data(), TString::Format(";%s [GeV/c^{2}];Yield", name.Data()).Data(), 100, min, max);
            TH1D *h2 = new TH1D(TString::Format("%s_ROOT", name.Data()).Data(),
                                TString::Format(";%s [GeV/c^{2}];Yield", name.Data()).Data(),
                                100,
                                min,
                                max);
            h->Sumw2();
            h2->Sumw2();
            H1D.insert(std::make_pair(index, h));
            H1D_ROOT.insert(std::make_pair(index, h2));
        }
    }

    // MCBooster
    for(GInt_t event = 0; event < nevents; event++) {
        for(GInt_t i = 0; i < masses.size(); i++) {
            for(GInt_t j = 0; j < masses.size(); j++) {
                if(j >= i)
                    continue;

                GInt_t index = i + j * masses.size();

                Vector4R p = GenEvents.fDaughters[j][event] + GenEvents.fDaughters[i][event];

                H1D[index]->Fill(p.mass(), GenEvents.fWeights[event]);
            }
        }
    }

    // ROOT
    TGenPhaseSpace phsp_root;
    TLorentzVector mother_root(0.0, 0.0, 0.0, mother_mass);
    phsp_root.SetDecay(mother_root, masses.size(), masses.data());
    std::vector<std::vector<TLorentzVector>> particles_root(masses.size(), std::vector<TLorentzVector>(nevents));
    std::vector<Double_t> weigts_root(nevents);
    clock_gettime(CLOCK_REALTIME, &time1);
    for(GInt_t event = 0; event < nevents; event++) {
        weigts_root[event] = phsp_root.Generate();
        for(GInt_t i = 0; i < masses.size(); i++) {
            TLorentzVector pr        = *phsp_root.GetDecay(i);
            particles_root[i][event] = pr;
        }
    }
    clock_gettime(CLOCK_REALTIME, &time2);
    GReal_t root_time_used = ((GReal_t)(time_diff(time1, time2).tv_sec + time_diff(time1, time2).tv_nsec * 1.0e-9));

    cout << "-----------------------------------------------------" << std::endl;
    cout << "----------------------- Timing ----------------------" << std::endl;
    cout << "Event generation in ROOT: " << root_time_used << std::endl;
    cout << "-----------------------------------------------------" << std::endl;
    cout << "----------------------- Timing ----------------------" << std::endl;

    for(GInt_t event = 0; event < nevents; event++) {
        for(GInt_t i = 0; i < masses.size(); i++) {
            for(GInt_t j = 0; j < masses.size(); j++) {
                if(j >= i)
                    continue;

                GInt_t index = i + j * masses.size();

                TLorentzVector pr = particles_root[i][event] + particles_root[j][event];
                H1D_ROOT[index]->Fill(pr.M(), weigts_root[event]);
            }
        }
    }

    TLegend *leg;

    for(GInt_t i = 0; i < masses.size(); i++) {
        for(GInt_t j = 0; j < masses.size(); j++) {
            if(j >= i)
                continue;
            GInt_t index = i + j * masses.size();
            TCanvas *c   = new TCanvas(H1D[index]->getName(), H1D[index]->getName(), 600, 500);

            H1D[index]->Draw("e0");
            H1D[index]->SetMarkerColor(kBlue);
            H1D[index]->SetMarkerSize(1.0);
            H1D[index]->SetMarkerStyle(8);
            H1D[index]->SetStats(0);

            H1D_ROOT[index]->Draw("HISTsame");
            H1D_ROOT[index]->SetLineColor(kRed);
            H1D_ROOT[index]->SetLineWidth(2);
            H1D_ROOT[index]->SetStats(0);

            leg = new TLegend(0.60, 0.8, 0.90, 0.9);
            leg->AddEntry(H1D[index], "mcbooster", "l");
            leg->AddEntry(H1D_ROOT[index], "TGenPhaseSpace", "l");
            leg->Draw();

            c->Print(TString::Format("./histo_%d%d.pdf", i, j));
        }
    }

    myapp->Run();
    return 0;
}
