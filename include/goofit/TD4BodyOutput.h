#pragma once

#ifdef ROOT_FOUND
#include <TFile.h>
#include <TTree.h>
#endif

#include <iostream>
#include <fstream>

#include <mcbooster/Vector4R.h>
#include <mcbooster/GContainers.h>

namespace GooFit {
class TD4BodyOutput final {
  public:
    TD4BodyOutput(const std::string &outFile, bool write4Vecs);
    virtual ~TD4BodyOutput();

    void writeEvent(double d0Mass,
                    const mcbooster::ParticlesSet_h &particles,
                    const mcbooster::VariableSet_h &variables,
                    unsigned int evtNum);

    void flushBuffer();
    bool flushNeeded();

  private:
    static const unsigned int _MAX_NUM_EVTS_BUFFER = 10000;
#ifdef ROOT_FOUND
    TFile _rootFile;
    TTree _tree;
#else
    static const std::string _DELIM;
    std::ofstream _csvFile;
    std::stringstream _ss;
#endif
    unsigned int _eventCount;
    bool _write4Vecs;
    double tm12_MeV_2, tm34_MeV_2, tc12_2, tc34_2, tphi_negPi_to_pi_2, tdtime_ps_2;
    double td0_Mass_MeV;
    double tkMinus_PE_MeV, tkMinus_Px_MeV, tkMinus_Py_MeV, tkMinus_Pz_MeV;
    double tpiPlus1_PE_MeV, tpiPlus1_Px_MeV, tpiPlus1_Py_MeV, tpiPlus1_Pz_MeV;
    double tpiPlus2_PE_MeV, tpiPlus2_Px_MeV, tpiPlus2_Py_MeV, tpiPlus2_Pz_MeV;
    double tpiMinus_PE_MeV, tpiMinus_Px_MeV, tpiMinus_Py_MeV, tpiMinus_Pz_MeV;
};
} // end namespace GooFit