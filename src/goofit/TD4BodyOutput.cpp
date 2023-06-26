#include <iostream>
#include <fstream>
#include <string>

#include <mcbooster/Vector4R.h>
#include <mcbooster/GContainers.h>

#include <goofit/TD4BodyOutput.h>

namespace GooFit
{
    TD4BodyOutput::TD4BodyOutput(const std::string &outFile, bool write4Vecs)
    {
        _write4Vecs = write4Vecs;
#ifdef ROOT_FOUND
        // build output root file/ntuple
        _rootFile = TFile(outFile, "RECREATE");

        _tree = TTree("events", "events");
        _tree->Branch("m12", &tm12_2, "m12/D");
        _tree->Branch("m34", &tm34_2, "m34/D");
        _tree->Branch("c12", &tc12_2, "c12/D");
        _tree->Branch("c34", &tc34_2, "c34/D");
        _tree->Branch("phi", &tphi_2, "phi/D");
        _tree->Branch("dtime", &tdtime_2, "dtime/D");
        if (_write4Vecs)
        {
            _tree->Branch("d0_Mass_GeV", &td0_Mass_GeV, "d0_Mass_GeV/D");
            _tree->Branch("kMinus_PE_GeV", &tkMinus_PE_GeV, "kMinus_PE_GeV/D");
            _tree->Branch("kMinus_Px_GeV", &tkMinus_Px_GeV, "kMinus_Px_GeV/D");
            _tree->Branch("kMinus_Py_GeV", &tkMinus_Py_GeV, "kMinus_Py_GeV/D");
            _tree->Branch("kMinus_Pz_GeV", &tkMinus_Pz_GeV, "kMinus_Pz_GeV/D");
            _tree->Branch("piPlus1_PE_GeV", &tpiPlus1_PE_GeV, "piPlus1_PE_GeV/D");
            _tree->Branch("piPlus1_Px_GeV", &tpiPlus1_Px_GeV, "piPlus1_Px_GeV/D");
            _tree->Branch("piPlus1_Py_GeV", &tpiPlus1_Py_GeV, "piPlus1_Py_GeV/D");
            _tree->Branch("piPlus1_Pz_GeV", &tpiPlus1_Pz_GeV, "piPlus1_Pz_GeV/D");
            _tree->Branch("piPlus2_PE_GeV", &tpiPlus2_PE_GeV, "piPlus2_PE_GeV/D");
            _tree->Branch("piPlus2_Px_GeV", &tpiPlus2_Px_GeV, "piPlus2_Px_GeV/D");
            _tree->Branch("piPlus2_Py_GeV", &tpiPlus2_Py_GeV, "piPlus2_Py_GeV/D");
            _tree->Branch("piPlus2_Pz_GeV", &tpiPlus2_Pz_GeV, "piPlus2_Pz_GeV/D");
            _tree->Branch("piMinus_PE_GeV", &tpiMinus_PE_GeV, "piMinus_PE_GeV/D");
            _tree->Branch("piMinus_Px_GeV", &tpiMinus_Px_GeV, "piMinus_Px_GeV/D");
            _tree->Branch("piMinus_Py_GeV", &tpiMinus_Py_GeV, "piMinus_Py_GeV/D");
            _tree->Branch("piMinus_Pz_GeV", &tpiMinus_Pz_GeV, "piMinus_Pz_GeV/D");
        } // end if writing 4 vecs
#else
        _csvFile.open(outFile, std::ios::out | std::ios::trunc);
        // TODO check units, particle order
        _csvFile << "m12, m34, c12, c34, phi, t";
        if (_write4Vecs)
        {
            _csvFile << ", D mass, ";
            _csvFile << "K PE [GeV], K Px [GeV], K Py [GeV], K Pz [GeV], ";
            _csvFile << "OS1 Pi PE [GeV], OS1 Pi Px [GeV],  OS1 Pi Py [GeV],  OS1 Pi Pz [GeV], ";
            _csvFile << "OS2 Pi PE [GeV], OS2 Pi Px [GeV],  OS2 Pi Py [GeV],  OS2 Pi Pz [GeV], ";
            _csvFile << "SS Pi PE [GeV], SS Pi Px [GeV],  SS Pi Py [GeV],  SS Pi Pz [GeV]";
        }
        _csvFile << std::endl;
        _csvFile.flush();
#endif
    }

    TD4BodyOutput::~TD4BodyOutput()
    {
#ifdef ROOT_FOUND
        _tree->Write();
        _rootFile->Close();
#else
        _csvFile.close();
#endif
    }

    void TD4BodyOutput::writeEvent(
        double d0Mass,
        const mcbooster::ParticlesSet_h& particles,
        const mcbooster::VariableSet_h& variables,
        unsigned int evtNum)
    {
        // get the values for this event
        tm12_2 = (*(variables[0]))[evtNum];
        tm34_2 = (*(variables[1]))[evtNum];
        tc12_2 = (*(variables[2]))[evtNum];
        tc34_2 = (*(variables[3]))[evtNum];
        tphi_2 = (*(variables[4]))[evtNum];
        tdtime_2 = (*(variables[5]))[evtNum];
        if (_write4Vecs)
        {
            td0_Mass_GeV = d0Mass;
            mcbooster::Vector4R piPlus1 = (*(particles[0]))[evtNum];
            tpiPlus1_PE_GeV = piPlus1.get(0);
            tpiPlus1_Px_GeV = piPlus1.get(1);
            tpiPlus1_Py_GeV = piPlus1.get(2);
            tpiPlus1_Pz_GeV = piPlus1.get(3);
            mcbooster::Vector4R piMinus = (*(particles[1]))[evtNum];
            tpiMinus_PE_GeV = piMinus.get(0);
            tpiMinus_Px_GeV = piMinus.get(1);
            tpiMinus_Py_GeV = piMinus.get(2);
            tpiMinus_Pz_GeV = piMinus.get(3);
            mcbooster::Vector4R kMinus = (*(particles[2]))[evtNum];
            tkMinus_PE_GeV = kMinus.get(0);
            tkMinus_Px_GeV = kMinus.get(1);
            tkMinus_Py_GeV = kMinus.get(2);
            tkMinus_Pz_GeV = kMinus.get(3);
            mcbooster::Vector4R piPlus2 = (*(particles[3]))[evtNum];
            tpiPlus2_PE_GeV = piPlus2.get(0);
            tpiPlus2_Px_GeV = piPlus2.get(1);
            tpiPlus2_Py_GeV = piPlus2.get(2);
            tpiPlus2_Pz_GeV = piPlus2.get(3);
        } // end if writing 4 vecs

#ifdef ROOT_FOUND
        _tree->Fill();
#else
        _csvFile << tm12_2 << ", " << tm34_2 << ", " << tc12_2 << ", " << tc34_2 << ", " << tphi_2 << ", " << tdtime_2;
        if (_write4Vecs)
        {
            _csvFile << ", ";
            _csvFile << td0_Mass_GeV << ", ";
            _csvFile << tkMinus_PE_GeV << ", " << tkMinus_Px_GeV << ", " << tkMinus_Py_GeV << ", " << tkMinus_Pz_GeV << ", ";
            _csvFile << tpiPlus1_PE_GeV << ", " << tpiPlus1_Px_GeV << ", " << tpiPlus1_Py_GeV << ", " << tpiPlus1_Pz_GeV << ", ";
            _csvFile << tpiPlus2_PE_GeV << ", " << tpiPlus2_Px_GeV << ", " << tpiPlus2_Py_GeV << ", " << tpiPlus2_Pz_GeV << ", ";
            _csvFile << tpiMinus_PE_GeV << ", " << tpiMinus_Px_GeV << ", " << tpiMinus_Py_GeV << ", " << tpiMinus_Pz_GeV;
        }
        _csvFile << std::endl;
        _csvFile.flush();
#endif
    }
} // end namespace GooFit