#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "ChristophModel.h"

// GooFit stuff
#include "goofit/FitManager.h"
#include "goofit/PDFs/basic/PolynomialPdf.h"
#include "goofit/PDFs/combine/AddPdf.h"
#include "goofit/PDFs/physics/Tddp4Pdf.h"
#include "goofit/PDFs/physics/TruthResolution_Aux.h"
#include "goofit/UnbinnedDataSet.h"

namespace GooFit {
const fptype ChristophModel::D0_MASS = 1.8645;
const fptype ChristophModel::PI_MASS = 0.13957018;
const fptype ChristophModel::K_MASS = 0.493677;
const fptype ChristophModel::D0_MESON_RADUIS = 5.0;
const fptype ChristophModel::D0_TAU = 0.4101;
const fptype ChristophModel::SQ_WS_TO_RS_RATE = 1.0 / sqrt(300.0);

ChristophModel::ChristophModel(const fptype xMixingValue,
                               const fptype yMixingValue,
                               const unsigned int modelMCEventsNorm,
                               bool special_integral)
    : _dk3piTau("tau", ChristophModel::D0_TAU),
      _dk3piXMixing("xmixing", xMixingValue),
      _dk3piYMixing("ymixing", yMixingValue),
      _dk3piSqWStoRSrate("SqWStoRSrate", ChristophModel::SQ_WS_TO_RS_RATE),
      _dk3piDecayInfo(_dk3piTau, _dk3piXMixing, _dk3piYMixing,
                      _dk3piSqWStoRSrate),
      _eff("constantEff", _modelVars, _coefficients, _offsets, 0) {
  _dk3piDecayInfo.meson_radius = ChristophModel::D0_MESON_RADUIS;
  _dk3piDecayInfo.particle_masses = {
      ChristophModel::D0_MASS, ChristophModel::PI_MASS, ChristophModel::PI_MASS,
      ChristophModel::K_MASS, ChristophModel::PI_MASS};
  _dk3piDecayInfo.amplitudes_B.insert(std::end(_dk3piDecayInfo.amplitudes_B),
                                      std::begin(_rs_amplitudes),
                                      std::end(_rs_amplitudes));
  _dk3piDecayInfo.amplitudes.insert(std::end(_dk3piDecayInfo.amplitudes),
                                    std::begin(_ws_amplitudes),
                                    std::end(_ws_amplitudes));

  _dp = new TDDP4("test_TD", _modelVars, _dk3piDecayInfo, &_dat, &_eff, 0,
                  modelMCEventsNorm,special_integral);
}

void ChristophModel::setXMixingRangeForFit(const fptype error,
                                           const fptype lowerLimit,
                                           const fptype upperLimit) {
  // constructing a Variable without giving a range, as is done when
  // initializing xmixing in the constructor,
  // sets fixed = true by default but we want to float xmixing if fitting
  _dk3piDecayInfo._xmixing.setFixed(false);
  _dk3piDecayInfo._xmixing.setError(error);
  _dk3piDecayInfo._xmixing.setLowerLimit(lowerLimit);
  _dk3piDecayInfo._xmixing.setUpperLimit(upperLimit);
}

void ChristophModel::setYMixingRangeForFit(const fptype error,
                                           const fptype lowerLimit,
                                           const fptype upperLimit) {
  // constructing a Variable without giving a range, as is done when
  // initializing ymixing in the constructor,
  // sets fixed = true by default but we want to float ymixing if fitting
  _dk3piDecayInfo._ymixing.setFixed(false);
  _dk3piDecayInfo._ymixing.setError(error);
  _dk3piDecayInfo._ymixing.setLowerLimit(lowerLimit);
  _dk3piDecayInfo._ymixing.setUpperLimit(upperLimit);
}

void ChristophModel::setModelMaxWeight(const fptype wmax) {
  _dp->setMaxWeight(wmax);
}

void ChristophModel::setGenerationOffset(const unsigned int generationOffset) {
  _dp->setGenerationOffset(generationOffset);
}

std::tuple<mcbooster::ParticlesSet_h, mcbooster::VariableSet_h,
           mcbooster::RealVector_h, mcbooster::BoolVector_h>
ChristophModel::generateSig(const int batchSize, const int seed) {
  return _dp->GenerateSig(batchSize, seed);
}

void ChristophModel::addEventToCurrentDataToFit(double m12, double m34,
                                                double cos12, double cos34,
                                                double phi, double dt,
                                                double sigmaT, int eventNum, double eff) {
  _model_m12.setValue(m12);
  _model_m34.setValue(m34);
  _model_cos12.setValue(cos12);
  _model_cos34.setValue(cos34);
  _model_phi.setValue(phi);
  _model_eventNumber.setValue(eventNum);
  _model_dtime.setValue(dt);
  _model_sigmat.setValue(sigmaT);
  _model_eff.setValue(eff);
  _currentDataToFit.addEvent();

  // std::cout<<"Added event with:"<<std::endl;
  // std::cout<<"    m12: "<<m12<<std::endl;
  // std::cout<<"    m34: "<<m34<<std::endl;
  // std::cout<<"  cos12: "<<cos12<<std::endl;
  // std::cout<<"  cos34: "<<cos34<<std::endl;
  // std::cout<<"    phi: "<<phi<<std::endl;
  // std::cout<<"     dt: "<<dt<<std::endl;
  // std::cout<<" sigmaT: "<<sigmaT<<std::endl;
  // std::cout<<"event #: "<<eventNum<<std::endl<<std::endl;

  // FIXME do we need to reset the Variable values after we are done adding the
  // event?
  // it seems the DataSet has to be backed by the same Variables as the model
}

void ChristophModel::fitCurrentData(unsigned int sampleNum,
                                    const std::string &outFile) {
  // build pdf
  std::cout << "Building PDF..." << std::endl;
  Variable constant("constant1", 1.0);
  Variable constant2("constant2", 1.0);
  std::vector<Variable> backgrVars = {constant};
  PolynomialPdf backgr("backgr", _model_m12, backgrVars);
  AddPdf signal("signal", constant2, _dp, &backgr);

  // set current data
  std::cout << "Setting current data..." << std::endl;
  signal.setData(&_currentDataToFit);
  _dp->setDataSize(_currentDataToFit.getNumEvents(), 9);

  // do fitting
  std::cout << "Fitting data (" << _currentDataToFit.getNumEvents()
            << " events)..." << std::endl;
  FitManager datapdf(&signal);
  // datapdf.setMaxCalls(10000);
  auto fitResults = datapdf.fit();

  // get results
  std::cout << "Retrieving fit results..." << std::endl;
  int convergedStatus =
      fitResults.IsValid() ? 1 : 0; // convert status to int to make it easier
                                    // to read/write fit log; 1 = converged, 0 =
                                    // did not converge
  // covStatus: 0 = not calculated, 1 = approx but not accurate, 2 = full matrix
  // but forced pos def, 3 = full accurate cov matrix according to minuit2 doc
  int covStatus = fitResults.UserState().CovarianceStatus();

  // write results to file
  // format is one line with dataSampleNum fittedXMixVal fittedXMixError
  // fittedYMixVal fittedYMixError convergedStatus covStatus
  std::cout << "Appending fit results to file " << outFile << "..."
            << std::endl;
  std::ofstream out;
  out.open(outFile.c_str(), std::ios::app);
  out.precision(10);
  out << sampleNum << " " << _dk3piDecayInfo._xmixing.getValue() << " "
      << _dk3piDecayInfo._xmixing.getError() << " "
      << _dk3piDecayInfo._ymixing.getValue() << " "
      << _dk3piDecayInfo._ymixing.getError() << " " << convergedStatus << " "
      << covStatus << std::endl;
  out.close();
}

ChristophModel::~ChristophModel() { delete _dp; }
} // end namespace GooFit
