#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "ToyModel.h"

#include "goofit/FitManager.h"
#include "goofit/PDFs/basic/PolynomialPdf.h"
#include "goofit/PDFs/combine/AddPdf.h"
#include "goofit/PDFs/physics/Tddp4Pdf.h"
#include "goofit/PDFs/physics/TruthResolution_Aux.h"
#include "goofit/UnbinnedDataSet.h"

namespace GooFit {

const fptype ToyModel::D0_MASS          = 1.8645;
const fptype ToyModel::PI_MASS          = 0.13957018;
const fptype ToyModel::K_MASS           = 0.493677;
const fptype ToyModel::D0_MESON_RADIUS  = 5.0;
const fptype ToyModel::D0_TAU           = 0.4101;
const fptype ToyModel::SQ_WS_TO_RS_RATE = 1.0 / sqrt(300.0);

ToyModel::ToyModel(const fptype xMixingValue, const fptype yMixingValue, const unsigned int modelMCEventsNorm)
    : _tau("tau", ToyModel::D0_TAU)
    , _xmixing("xmixing", xMixingValue)
    , _ymixing("ymixing", yMixingValue)
    , _sqWStoRSrate("SqWStoRSrate", ToyModel::SQ_WS_TO_RS_RATE)
    , _decayInfo(_tau, _xmixing, _ymixing, _sqWStoRSrate)
    , _eff("constantEff", _vars, _coefficients, _offsets, 0) {
    _decayInfo.meson_radius = ToyModel::D0_MESON_RADIUS;
    _decayInfo.particle_masses
        = {ToyModel::D0_MASS, ToyModel::PI_MASS, ToyModel::PI_MASS, ToyModel::K_MASS, ToyModel::PI_MASS};
    _decayInfo.amplitudes_B.insert(
        std::end(_decayInfo.amplitudes_B), std::begin(_cf_amplitudes), std::end(_cf_amplitudes));
    _decayInfo.amplitudes.insert(
        std::end(_decayInfo.amplitudes), std::begin(_dcs_amplitudes), std::end(_dcs_amplitudes));
    std::string pdf_name = "test_TD";
    _dp = new Amp4Body_TD(pdf_name, _vars, _decayInfo, &_dat, &_eff, &_mistag, (long)0, modelMCEventsNorm);
}

ToyModel::ToyModel(const fptype xMixingValue,
                   const fptype yMixingValue,
                   const std::vector<NormEvents_4Body_Base *> normEvents)
    : _tau("tau", ToyModel::D0_TAU)
    , _xmixing("xmixing", xMixingValue)
    , _ymixing("ymixing", yMixingValue)
    , _sqWStoRSrate("SqWStoRSrate", ToyModel::SQ_WS_TO_RS_RATE)
    , _decayInfo(_tau, _xmixing, _ymixing, _sqWStoRSrate)
    , _eff("constantEff", _vars, _coefficients, _offsets, 0) {
    _decayInfo.meson_radius = ToyModel::D0_MESON_RADIUS;
    _decayInfo.particle_masses
        = {ToyModel::D0_MASS, ToyModel::PI_MASS, ToyModel::PI_MASS, ToyModel::K_MASS, ToyModel::PI_MASS};
    _decayInfo.amplitudes_B.insert(
        std::end(_decayInfo.amplitudes_B), std::begin(_cf_amplitudes), std::end(_cf_amplitudes));
    _decayInfo.amplitudes.insert(
        std::end(_decayInfo.amplitudes), std::begin(_dcs_amplitudes), std::end(_dcs_amplitudes));
    std::string pdf_name = "test_TD";
    _dp                  = new Amp4Body_TD(pdf_name, _vars, _decayInfo, &_dat, &_eff, &_mistag, normEvents);
}

void ToyModel::setXMixingRangeForFit(const fptype error, const fptype lowerLimit, const fptype upperLimit) {
    _decayInfo._xmixing.setFixed(false);
    _decayInfo._xmixing.setError(error);
    _decayInfo._xmixing.setLowerLimit(lowerLimit);
    _decayInfo._xmixing.setUpperLimit(upperLimit);
}

void ToyModel::setYMixingRangeForFit(const fptype error, const fptype lowerLimit, const fptype upperLimit) {
    _decayInfo._ymixing.setFixed(false);
    _decayInfo._ymixing.setError(error);
    _decayInfo._ymixing.setLowerLimit(lowerLimit);
    _decayInfo._ymixing.setUpperLimit(upperLimit);
}

void ToyModel::setModelMaxWeight(const fptype wmax) { _dp->setMaxWeight(wmax); }

std::tuple<mcbooster::ParticlesSet_h, mcbooster::VariableSet_h, mcbooster::RealVector_h, mcbooster::BoolVector_h>
ToyModel::generateSig(const int batchSize, const int seed) {
    return _dp->GenerateSig(batchSize, seed);
}

void ToyModel::addEventToCurrentDataToFit(
    fptype m12, fptype m34, fptype cos12, fptype cos34, fptype phi, fptype dt, fptype sigmat, int eventNum) {
    _m12.setValue(m12);
    _m34.setValue(m34);
    _cos12.setValue(cos12);
    _cos34.setValue(cos34);
    _phi.setValue(phi);
    _eventNumber.setValue(eventNum);
    _dtime.setValue(dt);
    _sigmat.setValue(sigmat);
    _mistag.setValue(0.0);
    _currentDataToFit.addEvent();
}

void ToyModel::addEventToMCToPlot(
    fptype m12, fptype m34, fptype cos12, fptype cos34, fptype phi, fptype dt, fptype sigmat, int eventNum) {
    _m12.setValue(m12);
    _m34.setValue(m34);
    _cos12.setValue(cos12);
    _cos34.setValue(cos34);
    _phi.setValue(phi);
    _eventNumber.setValue(eventNum);
    _dtime.setValue(dt);
    _sigmat.setValue(sigmat);
    _mistag.setValue(0.0);
    _mcToPlot.addEvent();
}

void ToyModel::setGenerationOffset(const uint generationOffset) { _dp->setGenerationOffset(generationOffset); }

std::vector<std::vector<fptype>> ToyModel::fitCurrentData(unsigned int sampleNum, const std::string &outFile) {
    // Build the PDF.
    printf("Building PDF...\n");
    Variable constant("constant1", 1.0);
    Variable constant2("constant2", 1.0);
    std::vector<Variable> backgrVars = {constant};
    PolynomialPdf backgr("backgr", _m12, backgrVars);
    AddPdf signal("signal", constant2, _dp, &backgr);

    // Set the current data.
    printf("Setting current data...\n");
    signal.setData(&_currentDataToFit);
    std::cout << "Setting data size of Amp4Body_TD PDF" << std::endl;
    _dp->setDataSize(_currentDataToFit.getNumEvents(), 9);

    // Perform the fit.
    printf("Fitting data (%i events)...\n", _currentDataToFit.getNumEvents());
    FitManager datapdf(&signal);
    ROOT::Minuit2::FunctionMinimum fitResults = datapdf.fit();

    // Get the results.
    printf("Retrieving fit results...\n");
    int convergedStatus = fitResults.IsValid() ? 1 : 0;
    int covStatus       = fitResults.UserState().CovarianceStatus();
    std::ofstream out;
    out.open(outFile.c_str(), std::ios::app);
    out.precision(10);
    out << sampleNum << " " << _decayInfo._xmixing.getValue() << " " << _decayInfo._xmixing.getError() << " "
        << _decayInfo._ymixing.getValue() << " " << _decayInfo._ymixing.getError() << " " << convergedStatus << " "
        << covStatus << std::endl;
    out.close();

    // Evaluate PDF for MC events and return.
    std::cout << "Evaluating on integration MC..." << std::endl;
    signal.setData(&_mcToPlot);
    _dp->setDataSize(_mcToPlot.getNumEvents(), 9);
    return signal.getCompProbsAtDataPoints();
}

ToyModel::~ToyModel() { delete _dp; }

} // namespace GooFit
