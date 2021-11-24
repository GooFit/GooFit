#pragma once

#include <string>
#include <tuple>

#include "goofit/PDFs/basic/PolynomialPdf.h"
#include "goofit/PDFs/physics/DP4Pdf.h"
#include "goofit/PDFs/physics/DalitzPlotHelpers.h"
#include "goofit/PDFs/physics/Tddp4Pdf.h"
#include "goofit/PDFs/physics/TruthResolution_Aux.h"
#include "goofit/UnbinnedDataSet.h"
#include "goofit/Variable.h"

namespace GooFit {
class ChristophModel final {
public:
  ChristophModel(const fptype xMixingValue, const fptype yMixingValue,
                 const unsigned int modelMCEventsNorm);
  ~ChristophModel();

  void setXMixingRangeForFit(const fptype error, const fptype lowerLimit,
                             const fptype upperLimit);

  void setYMixingRangeForFit(const fptype error, const fptype lowerLimit,
                             const fptype upperLimit);

  void setModelMaxWeight(const fptype wmax);

  void setGenerationOffset(const unsigned int generationOffset);

  std::tuple<mcbooster::ParticlesSet_h, mcbooster::VariableSet_h,
             mcbooster::RealVector_h, mcbooster::BoolVector_h>
  generateSig(const int batchSize, const int seed);

  void addEventToCurrentDataToFit(double m12, double m34, double cos12,
                                  double cos34, double phi, double dt,
                                  double sigmaT, int eventNum);

  void fitCurrentData(unsigned int sampleNum, const std::string &outFile);

private:
  static const fptype D0_MASS;
  static const fptype PI_MASS;
  static const fptype K_MASS;
  static const fptype D0_MESON_RADUIS;
  static const fptype D0_TAU;
  static const fptype SQ_WS_TO_RS_RATE;

  Variable _dk3piTau;
  Variable _dk3piXMixing;
  Variable _dk3piYMixing;
  Variable _dk3piSqWStoRSrate;
  DecayInfo4t _dk3piDecayInfo;

  Variable _rhoMass = Variable("rho_mass", 0.77526);
  Variable _rhoWidth = Variable("rho_width", 0.1478);
  Variable _k892M = Variable("K892M", 0.89581);
  Variable _k892W = Variable("K892W", 0.0474);
  Variable _f600M = Variable("f600M", 0.519);
  Variable _f600W = Variable("f600W", 0.454);
  Variable _a1M = Variable("a1M", 1.237);
  Variable _a1W = Variable("a1W", 0.526);
  Variable _k1_1270M = Variable("K1_1270M", 1.28241);
  Variable _k1_1270W = Variable("K1_1270W", 0.06596);
  Variable _k0_1430M = Variable("K0_1430M", 1.425);
  Variable _k0_1430W = Variable("K0_1430W", 0.27);
  Variable _k1410M = Variable("K1410M", 1.414);
  Variable _k1410W = Variable("K1410W", 0.232);
  Variable _rho1450M = Variable("rho1450M", 1.465);
  Variable _rho1450W = Variable("rho1450W", 0.400);
  Variable _k1460M = Variable("K1460M", 1.351);
  Variable _k1460W = Variable("K1460W", 0.281);
  Variable _f0_1370M = Variable("f0_1370M", 1.350);
  Variable _f0_1370W = Variable("f0_1370W", 0.35);
  Variable _k1_1400M = Variable("K1_1400M", 1.403);
  Variable _k1_1400W = Variable("K1_1400W", 0.174);
  Variable _k2_1430M = Variable("K2_1430M", 1.4256);
  Variable _k2_1430W = Variable("K2_1430W", 0.0985);

  Variable _lass_a = Variable("lass_a", 2.07);
  Variable _lass_r = Variable("lass_r", 3.32);
  Variable _lass_pf = Variable("lass_pf", 0.0);
  Variable _lass_pr = Variable("lass_pr", 0.0);
  Variable _lass_F = Variable("lass_F", 1.0);
  std::vector<Variable> _lassVars{_lass_a, _lass_r, _lass_pf, _lass_pr,
                                  _lass_F};

  // Spin factors: we have two due to the bose symmetrization of the two pi+
  std::vector<SpinFactor *> _sf_K892_rho770_S{
      new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S,
                     D0_MESON_RADUIS, 0, 1, 2, 3),
      new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S,
                     D0_MESON_RADUIS, 3, 1, 2, 0)};

  // Lineshapes, also for both pi+ configurations
  std::vector<Lineshape *> LS_K892_rho770_S{
      new Lineshapes::RBW("K*(892)bar", _k892M, _k892W, 1, M_34, FF::BL2),
      new Lineshapes::RBW("rho(770)", _rhoMass, _rhoWidth, 1, M_12, FF::BL2),
      new Lineshapes::RBW("K*(892)bar", _k892M, _k892W, 1, M_13, FF::BL2),
      new Lineshapes::RBW("rho(770)", _rhoMass, _rhoWidth, 1, M_24, FF::BL2)};

  std::vector<SpinFactor *> _sf_K892_rho770_P{
      new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_P,
                     D0_MESON_RADUIS, 0, 1, 2, 3),
      new SpinFactor("SF", SF_4Body::FF_12_34_L1, D0_MESON_RADUIS, 0, 1, 2, 3),
      new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_P,
                     D0_MESON_RADUIS, 3, 1, 2, 0),
      new SpinFactor("SF", SF_4Body::FF_12_34_L1, D0_MESON_RADUIS, 3, 1, 2, 0)};

  std::vector<Lineshape *> _ls_K892_rho770_P{
      new Lineshapes::RBW("rho(770)", _rhoMass, _rhoWidth, 1, M_12, FF::BL2),
      new Lineshapes::RBW("K*(892)bar", _k892M, _k892W, 1, M_34, FF::BL2),
      new Lineshapes::RBW("rho(770)", _rhoMass, _rhoWidth, 1, M_24, FF::BL2),
      new Lineshapes::RBW("K*(892)bar", _k892M, _k892W, 1, M_13, FF::BL2)};

  std::vector<SpinFactor *> _sf_K892_rho770_D{
      new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_D,
                     D0_MESON_RADUIS, 0, 1, 2, 3),
      new SpinFactor("SF", SF_4Body::FF_12_34_L2, D0_MESON_RADUIS, 0, 1, 2, 3),
      new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_D,
                     D0_MESON_RADUIS, 3, 1, 2, 0),
      new SpinFactor("SF", SF_4Body::FF_12_34_L2, D0_MESON_RADUIS, 3, 1, 2, 0)};

  std::vector<Lineshape *> _ls_K892_rho770_D{
      new Lineshapes::RBW("rho(770)", _rhoMass, _rhoWidth, 1, M_12, FF::BL2),
      new Lineshapes::RBW("K*(892)bar", _k892M, _k892W, 1, M_34, FF::BL2),
      new Lineshapes::RBW("rho(770)", _rhoMass, _rhoWidth, 1, M_24, FF::BL2),
      new Lineshapes::RBW("K*(892)bar", _k892M, _k892W, 1, M_13, FF::BL2)};

  std::vector<SpinFactor *> _sf_K1410_rho770_S{
      new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S,
                     D0_MESON_RADUIS, 0, 1, 2, 3),
      new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S,
                     D0_MESON_RADUIS, 3, 1, 2, 0)};

  std::vector<Lineshape *> _ls_K1410_rho770_S{
      new Lineshapes::RBW("rho(770)", _rhoMass, _rhoWidth, 1, M_12, FF::BL2),
      new Lineshapes::RBW("K*(1410)", _k1410M, _k1410W, 1, M_34, FF::BL2),
      new Lineshapes::RBW("rho(770)", _rhoMass, _rhoWidth, 1, M_24, FF::BL2),
      new Lineshapes::RBW("K*(1410)", _k1410M, _k1410W, 1, M_13, FF::BL2)};

  std::vector<SpinFactor *> _sf_K1410_rho770_P{
      new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_P,
                     D0_MESON_RADUIS, 0, 1, 2, 3),
      new SpinFactor("SF", SF_4Body::FF_12_34_L1, D0_MESON_RADUIS, 0, 1, 2, 3),
      new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_P,
                     D0_MESON_RADUIS, 3, 1, 2, 0),
      new SpinFactor("SF", SF_4Body::FF_12_34_L1, D0_MESON_RADUIS, 3, 1, 2, 0)};

  std::vector<Lineshape *> _ls_K1410_rho770_P{
      new Lineshapes::RBW("rho(770)", _rhoMass, _rhoWidth, 1, M_12, FF::BL2),
      new Lineshapes::RBW("K*(1410)", _k1410M, _k1410W, 1, M_34, FF::BL2),
      new Lineshapes::RBW("rho(770)", _rhoMass, _rhoWidth, 1, M_24, FF::BL2),
      new Lineshapes::RBW("K*(1410)", _k1410M, _k1410W, 1, M_13, FF::BL2)};

  std::vector<SpinFactor *> _sf_K892_f0_600{
      new SpinFactor("SF", SF_4Body::DtoVS_VtoP1P2_StoP3P4, D0_MESON_RADUIS, 2,
                     3, 0, 1),
      new SpinFactor("SF", SF_4Body::FF_12_34_L1, D0_MESON_RADUIS, 2, 3, 0, 1),
      new SpinFactor("SF", SF_4Body::DtoVS_VtoP1P2_StoP3P4, D0_MESON_RADUIS, 2,
                     0, 3, 1),
      new SpinFactor("SF", SF_4Body::FF_12_34_L1, D0_MESON_RADUIS, 2, 0, 3, 1)};

  std::vector<Lineshape *> _ls_K892_f0_600{
      new Lineshapes::RBW("K*(892)bar", _k892M, _k892W, 1, M_34, FF::BL2),
      new Lineshapes::Bugg3("f600", _f600M, _f600W, 0, M_12, FF::BL2),
      new Lineshapes::RBW("K*(892)bar", _k892M, _k892W, 1, M_13, FF::BL2),
      new Lineshapes::Bugg3("f600", _f600M, _f600W, 0, M_24, FF::BL2)};

  std::vector<SpinFactor *> _sf_rho1450_K0_1430{
      new SpinFactor("SF", SF_4Body::DtoVS_VtoP1P2_StoP3P4, D0_MESON_RADUIS, 0,
                     1, 2, 3),
      new SpinFactor("SF", SF_4Body::FF_12_34_L1, D0_MESON_RADUIS, 0, 1, 2, 3),
      new SpinFactor("SF", SF_4Body::DtoVS_VtoP1P2_StoP3P4, D0_MESON_RADUIS, 3,
                     1, 2, 0),
      new SpinFactor("SF", SF_4Body::FF_12_34_L1, D0_MESON_RADUIS, 3, 1, 2, 0)};

  std::vector<Lineshape *> _ls_rho1450_K0_1430{
      new Lineshapes::RBW("rho(1450)", _rho1450M, _rho1450W, 1, M_12, FF::BL2),
      new Lineshapes::GLASS("K(0)*(1430)", _k0_1430M, _k0_1430W, 0, M_34,
                            FF::BL2, 1.5, _lassVars),
      new Lineshapes::RBW("rho(1450)", _rho1450M, _rho1450W, 1, M_24, FF::BL2),
      new Lineshapes::GLASS("K(0)*(1430)", _k0_1430M, _k0_1430W, 0, M_13,
                            FF::BL2, 1.5, _lassVars)};

  std::vector<SpinFactor *> _sf_K1460_K892{
      new SpinFactor("SF", SF_4Body::DtoPP1_PtoVP2_VtoP3P4, D0_MESON_RADUIS, 0,
                     1, 2, 3),
      new SpinFactor("SF", SF_4Body::DtoPP1_PtoVP2_VtoP3P4, D0_MESON_RADUIS, 3,
                     1, 2, 0)};

  std::vector<Lineshape *> _ls_K1460_K892{
      new Lineshapes::RBW("K1460", _k1460M, _k1460W, 1, M_34_2, FF::BL2),
      new Lineshapes::RBW("K*(892)bar", _k892M, _k892W, 1, M_34, FF::BL2),
      new Lineshapes::RBW("K1460", _k1460M, _k1460W, 1, M_13_2, FF::BL2),
      new Lineshapes::RBW("K*(892)bar", _k892M, _k892W, 1, M_13, FF::BL2)};

  std::vector<SpinFactor *> _sf_K1460_f0_1370{
      new SpinFactor("SF", SF_4Body::DtoPP1_PtoSP2_StoP3P4, D0_MESON_RADUIS, 0,
                     1, 2, 3),
      new SpinFactor("SF", SF_4Body::DtoPP1_PtoSP2_StoP3P4, D0_MESON_RADUIS, 3,
                     1, 2, 0)};

  std::vector<Lineshape *> _ls_K1460_f0_1370{
      new Lineshapes::RBW("K1460", _k1460M, _k1460W, 0, M_12_3, FF::BL2),
      new Lineshapes::RBW("f0_1370", _f0_1370M, _f0_1370W, 0, M_12, FF::BL2),
      new Lineshapes::RBW("K1460", _k1460M, _k1460W, 0, M_24_3, FF::BL2),
      new Lineshapes::RBW("f0_1370", _f0_1370M, _f0_1370W, 0, M_24, FF::BL2)};

  std::vector<SpinFactor *> _sf_K1_1270_K892{
      new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, D0_MESON_RADUIS, 0,
                     1, 2, 3),
      new SpinFactor("SF", SF_4Body::FF_123_4_L1, D0_MESON_RADUIS, 1, 2, 3, 0),
      new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, D0_MESON_RADUIS, 3,
                     1, 2, 0),
      new SpinFactor("SF", SF_4Body::FF_123_4_L1, D0_MESON_RADUIS, 1, 2, 0, 3)};

  std::vector<Lineshape *> _ls_K1_1270_K892{
      new Lineshapes::RBW("K1_1270", _k1_1270M, _k1_1270W, 0, M_34_2, FF::BL2),
      new Lineshapes::RBW("K*(892)bar", _k892M, _k892W, 1, M_34, FF::BL2),
      new Lineshapes::RBW("K1_1270", _k1_1270M, _k1_1270W, 0, M_13_2, FF::BL2),
      new Lineshapes::RBW("K*(892)bar", _k892M, _k892W, 1, M_13, FF::BL2)};

  std::vector<SpinFactor *> _sf_K1_1270_rho770{
      new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, D0_MESON_RADUIS, 3,
                     2, 0, 1),
      new SpinFactor("SF", SF_4Body::FF_123_4_L1, D0_MESON_RADUIS, 0, 1, 2, 3),
      new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, D0_MESON_RADUIS, 0,
                     2, 3, 1),
      new SpinFactor("SF", SF_4Body::FF_123_4_L1, D0_MESON_RADUIS, 1, 2, 3, 0)};

  std::vector<Lineshape *> _ls_K1_1270_rho770{
      new Lineshapes::RBW("K1_1270", _k1_1270M, _k1_1270W, 0, M_12_3, FF::BL2),
      new Lineshapes::RBW("rho(770)", _rhoMass, _rhoWidth, 1, M_12, FF::BL2),
      new Lineshapes::RBW("K1_1270", _k1_1270M, _k1_1270W, 0, M_24_3, FF::BL2),
      new Lineshapes::RBW("rho(770)", _rhoMass, _rhoWidth, 1, M_24, FF::BL2)};

  std::vector<SpinFactor *> _sf_K1_1270_K0_1430{
      new SpinFactor("SF", SF_4Body::DtoAP1_AtoSP2_StoP3P4, D0_MESON_RADUIS, 0,
                     1, 2, 3),
      new SpinFactor("SF", SF_4Body::FF_123_4_L1, D0_MESON_RADUIS, 1, 2, 3, 0),
      new SpinFactor("SF", SF_4Body::DtoAP1_AtoSP2_StoP3P4, D0_MESON_RADUIS, 3,
                     1, 2, 0),
      new SpinFactor("SF", SF_4Body::FF_123_4_L1, D0_MESON_RADUIS, 0, 1, 2, 3)};

  std::vector<Lineshape *> _ls_K1_1270_K0_1430{
      new Lineshapes::RBW("K(1)(1270)bar", _k1_1270M, _k1_1270W, 1, M_34_2,
                          FF::BL2),
      new Lineshapes::GLASS("K(0)*(1430)", _k0_1430M, _k0_1430W, 0, M_34,
                            FF::BL2, 1.5, _lassVars),
      new Lineshapes::RBW("K(1)(1270)bar2", _k1_1270M, _k1_1270W, 1, M_13_2,
                          FF::BL2),
      new Lineshapes::GLASS("K(0)*1430)", _k0_1430M, _k0_1430W, 0, M_13,
                            FF::BL2, 1.5, _lassVars)};

  std::vector<SpinFactor *> _sf_K1_1400_K892{
      new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, D0_MESON_RADUIS, 0,
                     1, 2, 3),
      new SpinFactor("SF", SF_4Body::FF_123_4_L1, D0_MESON_RADUIS, 1, 2, 3, 0),
      new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, D0_MESON_RADUIS, 3,
                     1, 2, 0),
      new SpinFactor("SF", SF_4Body::FF_123_4_L1, D0_MESON_RADUIS, 1, 2, 0, 3)};

  std::vector<Lineshape *> _ls_K1_1400_K892{
      new Lineshapes::RBW("K1_1400", _k1_1400M, _k1_1400W, 0, M_34_2, FF::BL2),
      new Lineshapes::RBW("K*(892)bar", _k892M, _k892W, 1, M_34, FF::BL2),
      new Lineshapes::RBW("K1_1400", _k1_1400M, _k1_1400W, 0, M_13_2, FF::BL2),
      new Lineshapes::RBW("K*(892)bar", _k892M, _k892W, 1, M_13, FF::BL2)};

  std::vector<SpinFactor *> _sf_K2_1430_K892{
      new SpinFactor("SF", SF_4Body::DtoTP1_TtoVP2_VtoP3P4, D0_MESON_RADUIS, 0,
                     1, 2, 3),
      new SpinFactor("SF", SF_4Body::FF_123_4_L2, D0_MESON_RADUIS, 1, 2, 3, 0),
      new SpinFactor("SF", SF_4Body::DtoTP1_TtoVP2_VtoP3P4, D0_MESON_RADUIS, 3,
                     1, 2, 0),
      new SpinFactor("SF", SF_4Body::FF_123_4_L2, D0_MESON_RADUIS, 1, 2, 0, 3)};

  std::vector<Lineshape *> _ls_K2_1430_K892{
      new Lineshapes::RBW("K2_1430", _k2_1430M, _k2_1430W, 2, M_34_2, FF::BL2),
      new Lineshapes::RBW("K*(892)bar", _k892M, _k892W, 1, M_34, FF::BL2),
      new Lineshapes::RBW("K2_1430", _k2_1430M, _k2_1430W, 2, M_13_2, FF::BL2),
      new Lineshapes::RBW("K*(892)bar", _k892M, _k892W, 1, M_13, FF::BL2)};

  std::vector<SpinFactor *> _sf_K2_1430_rho770{
      new SpinFactor("SF", SF_4Body::DtoTP1_TtoVP2_VtoP3P4, D0_MESON_RADUIS, 3,
                     2, 0, 1),
      new SpinFactor("SF", SF_4Body::FF_123_4_L2, D0_MESON_RADUIS, 0, 1, 2, 3),
      new SpinFactor("SF", SF_4Body::DtoTP1_TtoVP2_VtoP3P4, D0_MESON_RADUIS, 0,
                     2, 3, 1),
      new SpinFactor("SF", SF_4Body::FF_123_4_L2, D0_MESON_RADUIS, 3, 1, 2, 0)};

  std::vector<Lineshape *> _ls_K2_1430_rho770{
      new Lineshapes::RBW("K2_1430", _k2_1430M, _k2_1430W, 2, M_12_3, FF::BL2),
      new Lineshapes::RBW("rho(770)", _rhoMass, _rhoWidth, 1, M_12, FF::BL2),
      new Lineshapes::RBW("K2_1430", _k2_1430M, _k2_1430W, 2, M_24_3, FF::BL2),
      new Lineshapes::RBW("rho(770)", _rhoMass, _rhoWidth, 1, M_24, FF::BL2)};

  std::vector<SpinFactor *> _sf_a1_f0_600{
      new SpinFactor("SF", SF_4Body::DtoAP1_AtoSP2_StoP3P4, D0_MESON_RADUIS, 2,
                     3, 0, 1),
      new SpinFactor("SF", SF_4Body::FF_123_4_L1, D0_MESON_RADUIS, 0, 1, 3, 2),
      new SpinFactor("SF", SF_4Body::DtoAP1_AtoSP2_StoP3P4, D0_MESON_RADUIS, 2,
                     0, 3, 1),
      new SpinFactor("SF", SF_4Body::FF_123_4_L1, D0_MESON_RADUIS, 0, 1, 3, 2)};

  std::vector<Lineshape *> _ls_a1_f0_600{
      new Lineshapes::RBW("a(1)(1260)+", _a1M, _a1W, 1, M_12_4, FF::BL2, 5.71),
      new Lineshapes::Bugg3("f600", _f600M, _f600W, 0, M_12, FF::BL2),
      new Lineshapes::RBW("a(1)(1260)+", _a1M, _a1W, 1, M_24_1, FF::BL2, 5.71),
      new Lineshapes::Bugg3("f600", _f600M, _f600W, 0, M_24, FF::BL2)};

  std::vector<SpinFactor *> _sf_a1_rho770{
      new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, D0_MESON_RADUIS, 2,
                     3, 0, 1),
      new SpinFactor("SF", SF_4Body::FF_123_4_L1, D0_MESON_RADUIS, 0, 1, 3, 2),
      new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, D0_MESON_RADUIS, 2,
                     0, 3, 1),
      new SpinFactor("SF", SF_4Body::FF_123_4_L1, D0_MESON_RADUIS, 0, 1, 3, 2)};

  std::vector<Lineshape *> _ls_a1_rho770{
      new Lineshapes::RBW("a(1)(1260)+", _a1M, _a1W, 0, M_12_4, FF::BL2, 5.71),
      new Lineshapes::RBW("rho(770)", _rhoMass, _rhoWidth, 1, M_12, FF::BL2),
      new Lineshapes::RBW("a(1)(1260)+", _a1M, _a1W, 0, M_24_1, FF::BL2, 5.71),
      new Lineshapes::RBW("rho(770)", _rhoMass, _rhoWidth, 1, M_24, FF::BL2)};

  std::vector<SpinFactor *> _sf_a1_rho770_D{
      new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2Dwave_VtoP3P4,
                     D0_MESON_RADUIS, 2, 3, 0, 1),
      new SpinFactor("SF", SF_4Body::FF_123_4_L1, D0_MESON_RADUIS, 0, 1, 3, 2),
      new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2Dwave_VtoP3P4,
                     D0_MESON_RADUIS, 2, 0, 3, 1),
      new SpinFactor("SF", SF_4Body::FF_123_4_L1, D0_MESON_RADUIS, 0, 1, 3, 2)};

  std::vector<Lineshape *> _ls_a1_rho770_D{
      new Lineshapes::RBW("a(1)(1260)+", _a1M, _a1W, 2, M_12_4, FF::BL2, 5.71),
      new Lineshapes::RBW("rho(770)", _rhoMass, _rhoWidth, 1, M_12, FF::BL2),
      new Lineshapes::RBW("a(1)(1260)+", _a1M, _a1W, 2, M_24_1, FF::BL2, 5.71),
      new Lineshapes::RBW("rho(770)", _rhoMass, _rhoWidth, 1, M_24, FF::BL2)};

  std::vector<SpinFactor *> _sf_nonRes{
      new SpinFactor("SF", SF_4Body::ONE, D0_MESON_RADUIS, 2, 3, 0, 1),
      new SpinFactor("SF", SF_4Body::ONE, D0_MESON_RADUIS, 2, 0, 3, 1)};

  std::vector<Lineshape *> _ls_nonRes{
      new Lineshapes::One("nonRes", _a1M, _a1W, 0, M_12, FF::BL2),
      new Lineshapes::One("nonRes", _rhoMass, _rhoWidth, 0, M_34, FF::BL2),
      new Lineshapes::One("nonRes", _a1M, _a1W, 0, M_12, FF::BL2),
      new Lineshapes::One("nonRes", _rhoMass, _rhoWidth, 0, M_34, FF::BL2)};

  std::vector<SpinFactor *> _sf_NonResA_K892{
      new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2Dwave_VtoP3P4,
                     D0_MESON_RADUIS, 0, 1, 2, 3),
      new SpinFactor("SF", SF_4Body::FF_123_4_L1, D0_MESON_RADUIS, 1, 2, 3, 0),
      new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2Dwave_VtoP3P4,
                     D0_MESON_RADUIS, 3, 1, 2, 0),
      new SpinFactor("SF", SF_4Body::FF_123_4_L1, D0_MESON_RADUIS, 1, 2, 0, 3)};

  Variable _nr1 = Variable("NR1", 0.0);
  Variable _nr2 = Variable("NR2", 0.0);
  Variable _nr3 = Variable("NR3", 0.0);
  Variable _nr4 = Variable("NR4", 0.0);
  std::vector<Lineshape *> _ls_NonResA_K892{
      new Lineshapes::NonRes("K1_1400", _nr1, _nr2, 2, M_34_2, FF::BL2),
      new Lineshapes::RBW("K*(892)bar", _k892M, _k892W, 1, M_34, FF::BL2),
      new Lineshapes::NonRes("K1_1400", _nr3, _nr4, 2, M_13_2, FF::BL2),
      new Lineshapes::RBW("K*(892)bar", _k892M, _k892W, 1, M_13, FF::BL2)};

  // the very last parameter means that we have two permutations. so the first
  // half of the Lineshapes
  // and the first half of the spinfactors are amplitude 1, rest is amplitude 2
  // This means that it is important for symmetrized amplitueds that the
  // spinfactors and lineshapes are in the "right" order

  // RS Model
  Variable _k892_rho770_S_real = Variable("K892_rho770_S_real", 1.0);
  Variable _k892_rho770_S_imag = Variable("K892_rho770_S_imag", 0.0);
  Variable _k892_rho770_P_real = Variable("K892_rho770_P_real", 1.531);
  Variable _k892_rho770_P_imag = Variable("K892_rho770_P_imag", -0.886);
  Variable _k892_rho770_D_real = Variable("K892_rho770_D_real", 20.513);
  Variable _k892_rho770_D_imag = Variable("K892_rho770_D_imag", 20.398);
  Variable _k1410_rho770_P_real = Variable("K1410_rho770_P_real", 4.001);
  Variable _k1410_rho770_P_imag = Variable("K1410_rho770_P_imag", -2.620);
  Variable _k892_f0600_real = Variable("K892_f0600_real", -0.770);
  Variable _k892_f0600_imag = Variable("K892_f0600_imag", -1.530);
  Variable _rho1450_K0_1430_real = Variable("rho1450_K0_1430_real", -0.110);
  Variable _rho1450_K0_1430_imag = Variable("rho1450_K0_1430_imag", 1.814);
  Variable _k1460_K892_real = Variable("K1460_K892_real", -0.696);
  Variable _k1460_K892_imag = Variable("K1460_K892_imag", 0.326);
  Variable _k1460_f0_1370_real = Variable("K1460_f0_1370_real", -0.849);
  Variable _k1460_f0_1370_imag = Variable("K1460_f0_1370_imag", 0.972);
  Variable _k1_1270_K892_real = Variable("K1_1270_K892_real", 0.601);
  Variable _k1_1270_K892_imag = Variable("K1_1270_K892_imag", -0.182);
  Variable _k1_1270_rho770_real = Variable("K1_1270_rho770_real", -1.523);
  Variable _k1_1270_rho770_imag = Variable("K1_1270_rho770_imag", 1.244);
  Variable _k1_1270_K0_1430_real = Variable("K1_1270_K0_1430_real", 0.248);
  Variable _k1_1270_K0_1430_imag = Variable("K1_1270_K0_1430_imag", -0.088);
  Variable _k1_1400_K892_real = Variable("K1_1400_K892_real", -0.808);
  Variable _k1_1400_K892_imag = Variable("K1_1400_K892_imag", -0.358);
  Variable _nonResA_K892_real = Variable("NonResA_K892_real", -15.322);
  Variable _nonResA_K892_imag = Variable("NonResA_K892_imag", -12.089);
  Variable _k2_1430_K892_real = Variable("K2_1430_K892_real", 17.008);
  Variable _k2_1430_K892_imag = Variable("K2_1430_K892_imag", -5.014);
  Variable _k2_1430_rho770_real = Variable("K2_1430_rho770_real", 13.039);
  Variable _k2_1430_rho770_imag = Variable("K2_1430_rho770_imag", -1.935);
  Variable _a1_rho770_real = Variable("a1_rho770_real", -0.639);
  Variable _a1_rho770_imag = Variable("a1_rho770_imag", -6.801);
  Variable _a1_f0_600_real = Variable("a1_f0_600_real", -0.062);
  Variable _a1_f0_600_imag = Variable("a1_f0_600_imag", 2.872);
  Variable _a1_rho770_D_real = Variable("a1_rho770_D_real", -9.465);
  Variable _a1_rho770_D_imag = Variable("a1_rho770_D_imag", 15.390);
  Variable _nonRes_real = Variable("nonRes_real", -0.265);
  Variable _nonRes_imag = Variable("nonRes_imag", -0.003);
  std::vector<Amplitude *> _rs_amplitudes{
      new Amplitude("K892_rho770_S", _k892_rho770_S_real, _k892_rho770_S_imag,
                    LS_K892_rho770_S, _sf_K892_rho770_S, 2),
      new Amplitude("K892_rho770_P", _k892_rho770_P_real, _k892_rho770_P_imag,
                    _ls_K892_rho770_P, _sf_K892_rho770_P, 2),
      new Amplitude("K892_rho770_D", _k892_rho770_D_real, _k892_rho770_D_imag,
                    _ls_K892_rho770_D, _sf_K892_rho770_D, 2),
      new Amplitude("K1410_rho770", _k1410_rho770_P_real, _k1410_rho770_P_imag,
                    _ls_K1410_rho770_P, _sf_K1410_rho770_P, 2),
      new Amplitude("K892_f0600", _k892_f0600_real, _k892_f0600_imag,
                    _ls_K892_f0_600, _sf_K892_f0_600, 2),
      new Amplitude("rho1450_K0_1430", _rho1450_K0_1430_real,
                    _rho1450_K0_1430_imag, _ls_rho1450_K0_1430,
                    _sf_rho1450_K0_1430, 2),
      new Amplitude("K1460_K892", _k1460_K892_real, _k1460_K892_imag,
                    _ls_K1460_K892, _sf_K1460_K892, 2),
      new Amplitude("K1460_f0_1370", _k1460_f0_1370_real, _k1460_f0_1370_imag,
                    _ls_K1460_f0_1370, _sf_K1460_f0_1370, 2),
      new Amplitude("K1_1270_K892", _k1_1270_K892_real, _k1_1270_K892_imag,
                    _ls_K1_1270_K892, _sf_K1_1270_K892, 2),
      new Amplitude("K1_1270_rho770", _k1_1270_rho770_real,
                    _k1_1270_rho770_imag, _ls_K1_1270_rho770,
                    _sf_K1_1270_rho770, 2),
      new Amplitude("K1_1270_K0_1430", _k1_1270_K0_1430_real,
                    _k1_1270_K0_1430_imag, _ls_K1_1270_K0_1430,
                    _sf_K1_1270_K0_1430, 2),
      new Amplitude("K1_1400_K892", _k1_1400_K892_real, _k1_1400_K892_imag,
                    _ls_K1_1400_K892, _sf_K1_1400_K892, 2),
      new Amplitude("NonResA_K892", _nonResA_K892_real, _nonResA_K892_imag,
                    _ls_NonResA_K892, _sf_NonResA_K892, 2),
      new Amplitude("K2_1430_K892", _k2_1430_K892_real, _k2_1430_K892_imag,
                    _ls_K2_1430_K892, _sf_K2_1430_K892, 2),
      new Amplitude("K2_1430_rho770", _k2_1430_rho770_real,
                    _k2_1430_rho770_imag, _ls_K2_1430_rho770,
                    _sf_K2_1430_rho770, 2),
      new Amplitude("a1_rho770", _a1_rho770_real, _a1_rho770_imag,
                    _ls_a1_rho770, _sf_a1_rho770, 2),
      new Amplitude("a1_f0_600", _a1_f0_600_real, _a1_f0_600_imag,
                    _ls_a1_f0_600, _sf_a1_f0_600, 2),
      new Amplitude("a1_rho770_D", _a1_rho770_D_real, _a1_rho770_D_imag,
                    _ls_a1_rho770_D, _sf_a1_rho770_D, 2),
      new Amplitude("nonRes", _nonRes_real, _nonRes_imag, _ls_nonRes,
                    _sf_nonRes, 2)};

  // WS Model
  Variable _ws_K892_rho770_S_real = Variable("WS_K892_rho770_S_real", 1.0);
  Variable _ws_K892_rho770_S_imag = Variable("WS_K892_rho770_S_imag", 0.0);
  Variable _ws_K892_rho770_P_real = Variable("WS_K892_rho770_P_real", -0.109);
  Variable _ws_K892_rho770_P_imag = Variable("WS_K892_rho770_P_imag", 1.653);
  Variable _ws_K892_rho770_D_real = Variable("WS_K892_rho770_D_real", 25.463);
  Variable _ws_K892_rho770_D_imag = Variable("WS_K892_rho770_D_imag", 2.662);
  Variable _ws_rho1450_K0_1430_real =
      Variable("WS_rho1450_K0_1430_real", 2.353);
  Variable _ws_rho1450_K0_1430_imag =
      Variable("WS_rho1450_K0_1430_imag", -0.234);
  Variable _ws_K1_1270_K892_real = Variable("WS_K1_1270_K892_real", -0.035);
  Variable _ws_K1_1270_K892_imag = Variable("WS_K1_1270_K892_imag", -1.405);
  Variable _ws_K1_1270_rho770_real = Variable("WS_K1_1270_rho770_real", 2.42);
  Variable _ws_K1_1270_rho770_imag = Variable("WS_K1_1270_rho770_imag", 2.471);
  Variable _ws_K1_1270_K0_1430_real =
      Variable("WS_K1_1270_K0_1430_real", -1.990);
  Variable _ws_K1_1270_K0_1430_imag =
      Variable("WS_K1_1270_K0_1430_imag", -2.105);
  Variable _ws_K1_1400_K892_real = Variable("WS_K1_1400_K892_real", -3.347);
  Variable _ws_K1_1400_K892_imag = Variable("WS_K1_1400_K892_imag", -2.237);
  Variable _ws_nonRes_real = Variable("WS_nonRes_real", -0.456 / 6);
  Variable _ws_nonRes_imag = Variable("WS_nonRes_imag", -0.041 / 6);
  std::vector<Amplitude *> _ws_amplitudes{
      new Amplitude("WS_K892_rho770_S", _ws_K892_rho770_S_real,
                    _ws_K892_rho770_S_imag, LS_K892_rho770_S, _sf_K892_rho770_S,
                    2),
      new Amplitude("WS_K892_rho770_P", _ws_K892_rho770_P_real,
                    _ws_K892_rho770_P_imag, _ls_K892_rho770_P,
                    _sf_K892_rho770_P, 2),
      new Amplitude("WS_K892_rho770_D", _ws_K892_rho770_D_real,
                    _ws_K892_rho770_D_imag, _ls_K892_rho770_D,
                    _sf_K892_rho770_D, 2),
      new Amplitude("WS_rho1450_K0_1430", _ws_rho1450_K0_1430_real,
                    _ws_rho1450_K0_1430_imag, _ls_rho1450_K0_1430,
                    _sf_rho1450_K0_1430, 2),
      new Amplitude("WS_K1_1270_K892", _ws_K1_1270_K892_real,
                    _ws_K1_1270_K892_imag, _ls_K1_1270_K892, _sf_K1_1270_K892,
                    2),
      new Amplitude("WS_K1_1270_rho770", _ws_K1_1270_rho770_real,
                    _ws_K1_1270_rho770_imag, _ls_K1_1270_rho770,
                    _sf_K1_1270_rho770, 2),
      new Amplitude("WS_K1_1270_K0_1430", _ws_K1_1270_K0_1430_real,
                    _ws_K1_1270_K0_1430_imag, _ls_K1_1270_K0_1430,
                    _sf_K1_1270_K0_1430, 2),
      new Amplitude("WS_K1_1400_K892", _ws_K1_1400_K892_real,
                    _ws_K1_1400_K892_imag, _ls_K1_1400_K892, _sf_K1_1400_K892,
                    2),
      new Amplitude("WS_nonRes", _ws_nonRes_real, _ws_nonRes_imag, _ls_nonRes,
                    _sf_nonRes, 2)};

  Variable _constantOne = Variable("constantOne", 1);
  std::vector<Variable> _coefficients{_constantOne};

  Variable _constantZero = Variable("constantZero", 0);
  std::vector<Variable> _offsets{_constantZero, _constantZero};

  Observable _model_m12 = Observable("m12", 0, 3);
  Observable _model_m34 = Observable("m34", 0, 3);
  Observable _model_cos12 = Observable("cos12", -1, 1);
  Observable _model_cos34 = Observable("m12", -1, 1);
  Observable _model_phi = Observable("phi", -3.5, 3.5);
  EventNumber _model_eventNumber = EventNumber("eventNumber");
  Observable _model_dtime = Observable("dtime", 0, 10);
  Observable _model_sigmat = Observable("sigmat", -3, 3);
  std::vector<Observable> _modelVars = {
      _model_m12, _model_m34,         _model_cos12, _model_cos34,
      _model_phi, _model_eventNumber, _model_dtime, _model_sigmat};

  UnbinnedDataSet _currentDataToFit = UnbinnedDataSet(_modelVars);

  TruthResolution _dat;
  PolynomialPdf _eff;
  TDDP4 *_dp;
};
} // end namespace GooFit
