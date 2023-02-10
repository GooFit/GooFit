#pragma once

#include <string>
#include <tuple>

#include "goofit/PDFs/basic/PolynomialPdf.h"
#include "goofit/PDFs/physics/DP4Pdf.h"
#include "goofit/PDFs/physics/DalitzPlotHelpers.h"
#include "goofit/PDFs/physics/Amp4Body_TD.h"
#include "goofit/PDFs/physics/TruthResolution_Aux.h"
#include "goofit/UnbinnedDataSet.h"
#include "goofit/Variable.h"

namespace GooFit {

class ToyModel final {
  public:
    ToyModel(const fptype xMixingValue,
             const fptype yMixingValue,
             const std::vector<NormEvents_4Body_Base *> normEvents);

    ToyModel(const fptype xMixingValue, const fptype yMixingValue, const unsigned int modelMCEventsNorm);

    ~ToyModel();

    void setXMixingRangeForFit(const fptype error, const fptype lowerLimit, const fptype upperLimit);

    void setYMixingRangeForFit(const fptype error, const fptype lowerLimit, const fptype upperLimit);

    void setModelMaxWeight(const fptype wmax);

    void setGenerationOffset(const uint generationOffset);

    std::tuple<mcbooster::ParticlesSet_h, mcbooster::VariableSet_h, mcbooster::RealVector_h, mcbooster::BoolVector_h>
    generateSig(const int batchSize, const int seed);

    void addEventToCurrentDataToFit(
        fptype m12, fptype m34, fptype cos12, fptype cos34, fptype phi, fptype dt, fptype sigmat, int eventNum);

    void addEventToMCToPlot(
        fptype m12, fptype m34, fptype cos12, fptype cos34, fptype phi, fptype dt, fptype sigmat, int eventNum);

    std::vector<std::vector<fptype>> fitCurrentData(unsigned int sampleNum, const std::string &outFile);

    static const fptype D0_MASS;
    static const fptype PI_MASS;
    static const fptype K_MASS;
    static const fptype D0_MESON_RADIUS;
    static const fptype D0_TAU;
    static const fptype SQ_WS_TO_RS_RATE;

  private:
    Variable _tau;
    Variable _xmixing;
    Variable _ymixing;
    Variable _sqWStoRSrate;
    DecayInfo4t _decayInfo;

    Variable _rhoMass  = Variable("rho_mass", 0.77526);
    Variable _rhoWidth = Variable("rho_width", 0.1478);
    Variable _k892M    = Variable("K892M", 0.89581);
    Variable _k892W    = Variable("K892W", 0.0474);

    std::vector<SpinFactor *> _sf_K892_rho770_S{
        new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S, D0_MASS, 0, 1, 2, 3),
        new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S, D0_MASS, 3, 1, 2, 0)};

    std::vector<Lineshape *> LS_K892_rho770_S{new Lineshapes::RBW("K*(892)bar", _k892M, _k892W, 1, M_34, FF::BL2),
                                              new Lineshapes::RBW("rho(770)", _rhoMass, _rhoWidth, 1, M_12, FF::BL2),
                                              new Lineshapes::RBW("K*(892)bar", _k892M, _k892W, 1, M_13, FF::BL2),
                                              new Lineshapes::RBW("rho(770)", _rhoMass, _rhoWidth, 1, M_24, FF::BL2)};

    // CF model.
    Variable _k892_rho770_S_real = Variable("K892_rho770_S_real", 1.0);
    Variable _k892_rho770_S_imag = Variable("K892_rho770_S_imag", 0.0);
    std::vector<Amplitude *> _cf_amplitudes{new Amplitude(
        "K892_rho770_S", _k892_rho770_S_real, _k892_rho770_S_imag, LS_K892_rho770_S, _sf_K892_rho770_S, 2)};

    // DCS model.
    Variable _dcs_K892_rho770_S_real = Variable("WS_K892_rho770_S_real", 1.0);
    Variable _dcs_K892_rho770_S_imag = Variable("WS_K892_rho770_S_imag", 0.0);
    std::vector<Amplitude *> _dcs_amplitudes{new Amplitude(
        "K892_rho770_S", _dcs_K892_rho770_S_real, _dcs_K892_rho770_S_imag, LS_K892_rho770_S, _sf_K892_rho770_S, 2)};

    // Efficiency.
    Variable _constantOne = Variable("constantOne", 1);
    std::vector<Variable> _coefficients{_constantOne};
    Variable _constantZero = Variable("constantZero", 0);
    std::vector<Variable> _offsets{_constantZero, _constantZero};

    // Observables.
    Observable _m12          = Observable("m12", 0, 3);
    Observable _m34          = Observable("m34", 0, 3);
    Observable _cos12        = Observable("cos12", -1, 1);
    Observable _cos34        = Observable("cos34", -1, 1);
    Observable _phi          = Observable("phi", -3.5, 3.5);
    EventNumber _eventNumber = EventNumber("eventNumber");
    Observable _dtime        = Observable("dtime", 0, 10);
    Observable _sigmat       = Observable("sigmat", -3, 3);
    Observable _mistag       = Observable("mistag", 0, 1);

    std::vector<Observable> _vars = {_m12, _m34, _cos12, _cos34, _phi, _eventNumber, _dtime, _sigmat, _mistag};

    UnbinnedDataSet _currentDataToFit = UnbinnedDataSet(_vars);
    UnbinnedDataSet _mcToPlot         = UnbinnedDataSet(_vars);

    TruthResolution _dat;
    PolynomialPdf _eff;

    Amp4Body_TD *_dp;
};

} // namespace GooFit
