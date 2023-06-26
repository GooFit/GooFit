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

  private:
#ifdef ROOT_FOUND
    TFile _rootFile;
    TTree _tree;
#else
    std::ofstream _csvFile;
#endif
    bool _write4Vecs;
    double tm12_2, tm34_2, tc12_2, tc34_2, tphi_2, tdtime_2;
    double td0_Mass_GeV;
    double tkMinus_PE_GeV, tkMinus_Px_GeV, tkMinus_Py_GeV, tkMinus_Pz_GeV;
    double tpiPlus1_PE_GeV, tpiPlus1_Px_GeV, tpiPlus1_Py_GeV, tpiPlus1_Pz_GeV;
    double tpiPlus2_PE_GeV, tpiPlus2_Px_GeV, tpiPlus2_Py_GeV, tpiPlus2_Pz_GeV;
    double tpiMinus_PE_GeV, tpiMinus_Px_GeV, tpiMinus_Py_GeV, tpiMinus_Pz_GeV;
};
} // end namespace GooFit