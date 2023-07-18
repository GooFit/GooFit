#include <iostream>
#include <fstream>
#include <string>

#include <mcbooster/Vector4R.h>
#include <mcbooster/GContainers.h>

#include <goofit/TD4BodyOutput.h>

namespace GooFit {
TD4BodyOutput::TD4BodyOutput(const std::string &outFile, bool write4Vecs) {
    _write4Vecs = write4Vecs;
#ifdef ROOT_FOUND
    // build output root file/ntuple
    _rootFile = TFile(outFile, "RECREATE");

    _tree = TTree("events", "events");
    _tree->Branch("m12_MeV", &tm12_MeV_2, "m12_MeV/D");
    _tree->Branch("m34_MeV", &tm34_MeV_2, "m34_MeV/D");
    _tree->Branch("c12", &tc12_2, "c12/D");
    _tree->Branch("c34", &tc34_2, "c34/D");
    _tree->Branch("phi_negPi_to_pi", &tphi_negPi_to_pi_2, "phi_negPi_to_pi/D");
    _tree->Branch("dtime_ps", &tdtime_ps_2, "dtime_ps/D");
    if(_write4Vecs) {
        _tree->Branch("d0_Mass_MeV", &td0_Mass_MeV, "d0_Mass_MeV/D");
        _tree->Branch("kMinus_PE_MeV", &tkMinus_PE_MeV, "kMinus_PE_MeV/D");
        _tree->Branch("kMinus_Px_MeV", &tkMinus_Px_MeV, "kMinus_Px_MeV/D");
        _tree->Branch("kMinus_Py_MeV", &tkMinus_Py_MeV, "kMinus_Py_MeV/D");
        _tree->Branch("kMinus_Pz_MeV", &tkMinus_Pz_MeV, "kMinus_Pz_MeV/D");
        _tree->Branch("piPlus1_PE_MeV", &tpiPlus1_PE_MeV, "piPlus1_PE_MeV/D");
        _tree->Branch("piPlus1_Px_MeV", &tpiPlus1_Px_MeV, "piPlus1_Px_MeV/D");
        _tree->Branch("piPlus1_Py_MeV", &tpiPlus1_Py_MeV, "piPlus1_Py_MeV/D");
        _tree->Branch("piPlus1_Pz_MeV", &tpiPlus1_Pz_MeV, "piPlus1_Pz_MeV/D");
        _tree->Branch("piPlus2_PE_MeV", &tpiPlus2_PE_MeV, "piPlus2_PE_MeV/D");
        _tree->Branch("piPlus2_Px_MeV", &tpiPlus2_Px_MeV, "piPlus2_Px_MeV/D");
        _tree->Branch("piPlus2_Py_MeV", &tpiPlus2_Py_MeV, "piPlus2_Py_MeV/D");
        _tree->Branch("piPlus2_Pz_MeV", &tpiPlus2_Pz_MeV, "piPlus2_Pz_MeV/D");
        _tree->Branch("piMinus_PE_MeV", &tpiMinus_PE_MeV, "piMinus_PE_MeV/D");
        _tree->Branch("piMinus_Px_MeV", &tpiMinus_Px_MeV, "piMinus_Px_MeV/D");
        _tree->Branch("piMinus_Py_MeV", &tpiMinus_Py_MeV, "piMinus_Py_MeV/D");
        _tree->Branch("piMinus_Pz_MeV", &tpiMinus_Pz_MeV, "piMinus_Pz_MeV/D");
    } // end if writing 4 vecs
#else
    _csvFile.open(outFile, std::ios::out | std::ios::trunc);
    // TODO check units, particle order
    _csvFile << "m12_MeV" << _DELIM << "m34_MeV" << _DELIM << "c12" << _DELIM << "c34" << _DELIM << "phi_negPi_to_pi" << _DELIM << "t_ps";
    if(_write4Vecs) {
        _csvFile << _DELIM << "D_mass_MeV" << _DELIM;
        _csvFile << "K_PE_MeV" << _DELIM << "K_Px_MeV" << _DELIM << "K_Py_MeV" << _DELIM << "K_Pz_MeV" << _DELIM;
        _csvFile << "OS1_Pi_PE_MeV" << _DELIM << "OS1_Pi_Px_MeV" << _DELIM << "OS1_Pi_Py_MeV" << _DELIM << "OS1_Pi_Pz_MeV" << _DELIM;
        _csvFile << "OS2_Pi_PE_MeV" << _DELIM << "OS2_Pi_Px_MeV" << _DELIM << "OS2_Pi_Py_MeV" << _DELIM << "OS2_Pi_Pz_MeV" << _DELIM;
        _csvFile << "SS_Pi_PE_MeV" << _DELIM << "SS_Pi_Px_MeV" << _DELIM << "SS_Pi_Py_MeV" << _DELIM << "SS_Pi_Pz_MeV"; 
    }
    _csvFile << std::endl;
    _csvFile.flush();
#endif
}

TD4BodyOutput::~TD4BodyOutput() {
#ifdef ROOT_FOUND
    _tree->Write();
    _rootFile->Close();
#else
    _csvFile << _ss.rdbuf();
    _csvFile.close();
#endif
}

void TD4BodyOutput::flushBuffer()
{
#ifdef ROOT_FOUND
    // do nothing
#else
    _csvFile << _ss.rdbuf();
    _csvFile.flush();
    _ss=std::stringstream();
#endif    
}

bool TD4BodyOutput::flushNeeded()
{
    if (_eventCount != 0 && _eventCount % _MAX_NUM_EVTS_BUFFER == 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void TD4BodyOutput::writeEvent(double d0Mass,
                               const mcbooster::ParticlesSet_h &particles,
                               const mcbooster::VariableSet_h &variables,
                               unsigned int evtNum) 
{
    _eventCount++;

    // get the values for this event
    tm12_MeV_2   = (*(variables[0]))[evtNum];
    tm34_MeV_2   = (*(variables[1]))[evtNum];
    tc12_2   = (*(variables[2]))[evtNum];
    tc34_2   = (*(variables[3]))[evtNum];
    tphi_negPi_to_pi_2   = (*(variables[4]))[evtNum];
    tdtime_ps_2 = (*(variables[5]))[evtNum];
    if(_write4Vecs) {
        td0_Mass_MeV                = d0Mass;
        mcbooster::Vector4R piPlus1 = (*(particles[0]))[evtNum];
        tpiPlus1_PE_MeV             = piPlus1.get(0);
        tpiPlus1_Px_MeV             = piPlus1.get(1);
        tpiPlus1_Py_MeV             = piPlus1.get(2);
        tpiPlus1_Pz_MeV             = piPlus1.get(3);
        mcbooster::Vector4R piMinus = (*(particles[1]))[evtNum];
        tpiMinus_PE_MeV             = piMinus.get(0);
        tpiMinus_Px_MeV             = piMinus.get(1);
        tpiMinus_Py_MeV             = piMinus.get(2);
        tpiMinus_Pz_MeV             = piMinus.get(3);
        mcbooster::Vector4R kMinus  = (*(particles[2]))[evtNum];
        tkMinus_PE_MeV              = kMinus.get(0);
        tkMinus_Px_MeV              = kMinus.get(1);
        tkMinus_Py_MeV              = kMinus.get(2);
        tkMinus_Pz_MeV              = kMinus.get(3);
        mcbooster::Vector4R piPlus2 = (*(particles[3]))[evtNum];
        tpiPlus2_PE_MeV             = piPlus2.get(0);
        tpiPlus2_Px_MeV             = piPlus2.get(1);
        tpiPlus2_Py_MeV             = piPlus2.get(2);
        tpiPlus2_Pz_MeV             = piPlus2.get(3);
    } // end if writing 4 vecs

#ifdef ROOT_FOUND
    _tree->Fill();
#else
    _ss << tm12_MeV_2 << _DELIM << tm34_MeV_2 << _DELIM << tc12_2 << _DELIM << tc34_2 << _DELIM << tphi_negPi_to_pi_2 << _DELIM << tdtime_ps_2;
    if(_write4Vecs) 
    {
        _ss << _DELIM << td0_Mass_MeV 
            << _DELIM << tkMinus_PE_MeV << _DELIM << tkMinus_Px_MeV << _DELIM << tkMinus_Py_MeV << _DELIM << tkMinus_Pz_MeV 
            << _DELIM << tpiPlus1_PE_MeV << _DELIM << tpiPlus1_Px_MeV << _DELIM << tpiPlus1_Py_MeV << _DELIM << tpiPlus1_Pz_MeV 
            << _DELIM << tpiPlus2_PE_MeV << _DELIM << tpiPlus2_Px_MeV << _DELIM << tpiPlus2_Py_MeV << _DELIM << tpiPlus2_Pz_MeV 
            << _DELIM << tpiMinus_PE_MeV << _DELIM << tpiMinus_Px_MeV << _DELIM << tpiMinus_Py_MeV << _DELIM << tpiMinus_Pz_MeV;
    }
    _ss << std::endl;
#endif
    if (flushNeeded())
    {
        flushBuffer();
    }
}
} // end namespace GooFit