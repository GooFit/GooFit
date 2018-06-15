#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

from goofit import *
import numpy as np
import os

print_goofit_info()

GDIR = os.path.dirname(os.path.abspath(__file__))
for i in range(2):
    GDIR = os.path.dirname(GDIR) # up to goofit dir

# Constants used in more than one PDF component.
_mD0       = 1.8645
piPlusMass = 0.13957018
KmMass     = .493677

m12 = Observable("m12", 0, 3)
m34 = Observable("m34", 0, 3)
cos12 = Observable("cos12", -1, 1)
cos34 = Observable("m12", -1, 1)
phi = Observable("phi", -3.5, 3.5)
eventNumber = EventNumber("eventNumber", 0, INT_MAX)

currData = UnbinnedDataSet(m12, m34, cos12, cos34, phi, eventNumber)

data = np.loadtxt(os.path.join(GDIR, "examples/DP4/ToyMC.txt"))
currData.from_numpy(data.T)

DK3P_DI = DecayInfo4()
DK3P_DI.meson_radius    = 1.5
DK3P_DI.particle_masses = (_mD0, piPlusMass, piPlusMass, KmMass, piPlusMass)

RhoMass = Variable("rho_mass", 0.77526)
RhoWidth = Variable("rho_width", 0.1478)
KstarM = Variable("KstarM", 0.89581)
KstarW = Variable("KstarW", 0.0474)
f600M = Variable("f600M", 0.519)
f600W = Variable("f600W", 0.454)
a1M = Variable("a1M", 1.23)
a1W = Variable("a1W", 0.42)
K1M = Variable("K1M", 1.272)
K1W = Variable("K1W", 0.09)
K1430M = Variable("K1430M", 1.414)
K1430W = Variable("K1430W", 0.29)

# Spin factors: we have two due to the bose symmetrization of the two pi+
SFKRS = (SpinFactor("SF", SF_4Body.DtoV1V2_V1toP1P2_V2toP3P4_S, _mD0, 0, 1, 2, 3),
         SpinFactor("SF", SF_4Body.DtoV1V2_V1toP1P2_V2toP3P4_S, _mD0, 3, 1, 2, 0))

SFKRP = (SpinFactor("SF", SF_4Body.DtoV1V2_V1toP1P2_V2toP3P4_P, _mD0, 0, 1, 2, 3),
         SpinFactor("SF", SF_4Body.DtoV1V2_V1toP1P2_V2toP3P4_P, _mD0, 3, 1, 2, 0))

SFKRD = (SpinFactor("SF", SF_4Body.DtoV1V2_V1toP1P2_V2toP3P4_D, _mD0, 0, 1, 2, 3),
         SpinFactor("SF", SF_4Body.DtoV1V2_V1toP1P2_V2toP3P4_D, _mD0, 3, 1, 2, 0))

SFKF = (SpinFactor("SF", SF_4Body.DtoVS_VtoP1P2_StoP3P4, _mD0, 2, 3, 0, 1),
        SpinFactor("SF", SF_4Body.DtoVS_VtoP1P2_StoP3P4, _mD0, 2, 0, 3, 1))

SFKK = (SpinFactor("SF", SF_4Body.DtoAP1_AtoSP2_StoP3P4, _mD0, 0, 1, 3, 2),
        SpinFactor("SF", SF_4Body.DtoAP1_AtoSP2_StoP3P4, _mD0, 3, 1, 0, 2))

SFK1R = (SpinFactor("SF", SF_4Body.DtoAP1_AtoVP2_VtoP3P4, _mD0, 3, 2, 0, 1),
         SpinFactor("SF", SF_4Body.DtoAP1_AtoVP2_VtoP3P4, _mD0, 0, 2, 3, 1))

SFA1R = (SpinFactor("SF", SF_4Body.DtoAP1_AtoVP2_VtoP3P4, _mD0, 2, 3, 0, 1),
         SpinFactor("SF", SF_4Body.DtoAP1_AtoVP2_VtoP3P4, _mD0, 2, 0, 3, 1))

SFA1RD = (SpinFactor("SF", SF_4Body.DtoAP1_AtoVP2Dwave_VtoP3P4, _mD0, 2, 3, 0, 1),
          SpinFactor("SF", SF_4Body.DtoAP1_AtoVP2Dwave_VtoP3P4, _mD0, 2, 0, 3, 1))

# Lineshapes, also for both pi+ configurations
LSKRS = (Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_12),
         Lineshapes.RBW("K*(892)bar", KstarM, KstarW, 1, M_34),
         Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_24),
         Lineshapes.RBW("K*(892)bar", KstarM, KstarW, 1, M_13))

LSKRP = (Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_12),
         Lineshapes.RBW("K*(892)bar", KstarM, KstarW, 1, M_34),
         Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_24),
         Lineshapes.RBW("K*(892)bar", KstarM, KstarW, 1, M_13))

LSKRD = (Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_12),
         Lineshapes.RBW("K*(892)bar", KstarM, KstarW, 1, M_34),
         Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_24),
         Lineshapes.RBW("K*(892)bar", KstarM, KstarW, 1, M_13))

LSKF = (Lineshapes.RBW("K*(892)bar", KstarM, KstarW, 1, M_34),
        Lineshapes.Bugg("f600", f600M, f600W, 0, M_12),
        Lineshapes.RBW("K*(892)bar", KstarM, KstarW, 1, M_13),
        Lineshapes.Bugg("f600", f600M, f600W, 0, M_24))

LSKK = (Lineshapes.SBW("K(1)(1270)bar", K1M, K1W, 1, M_34_2),
        Lineshapes.LASS("K(0)*(1430)bar", K1430M, K1430W, 0, M_34),
        Lineshapes.SBW("K(1)(1270)bar2", K1M, K1W, 1, M_13_2),
        Lineshapes.LASS("K(0)*(1430)bar2", K1430M, K1430W, 0, M_13))

LSK1R = (Lineshapes.SBW("K(1)(1270)bar", K1M, K1W, 0, M_12_3),
         Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_12),
         Lineshapes.SBW("K(1)(1270)bar", K1M, K1W, 0, M_24_3),
         Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_24))

LSA1R = (Lineshapes.SBW("a(1)(1260)+", a1M, a1W, 0, M_12_4),
         Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_12),
         Lineshapes.SBW("a(1)(1260)+", a1M, a1W, 0, M_24_1),
         Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_24))

LSA1RD = (Lineshapes.SBW("a(1)(1260)+", a1M, a1W, 2, M_12_4),
          Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_12),
          Lineshapes.SBW("a(1)(1260)+", a1M, a1W, 2, M_24_1),
          Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_24))

# the very last parameter means that we have two permutations. so the first half of the Lineshapes
# and the first half of the spinfactors are amplitude 1, rest is amplitude
# This means that it is important for symmetrized amplitueds that the spinfactors and lineshapes are in the "right"
# order

# Amplitudes with the correct fixed values.
# Amplitude* Bose_symmetrized_AMP_S = Amplitude( "K*(892)rho(770)_S", Variable("amp_real1", -0.115177), new
# Variable("amp_imag1", 0.153976), LSKRS, SFKRS, 2)
# Amplitude* Bose_symmetrized_AMP_P = Amplitude( "K*(892)rho(770)_P", Variable("amp_real2", -0.0298697),
# Variable("amp_imag2", -0.0722874), LSKRP, SFKRP, 2)
# Amplitude* Bose_symmetrized_AMP_D = Amplitude( "K*(892)rho(770)_D", Variable("amp_real3", -0.452212), new
# Variable("amp_imag3", 0.426521), LSKRD, SFKRD, 2)

# Amplitudes with floating slightly different values to be fitted.

Bose_symmetrized_AMP_S = Amplitude("K*(892)rho(770)_S",
                                 Variable("amp_real1", -0.1, 0.001, 0, 0),
                                 Variable("amp_imag1", 0.1, 0.001, 0, 0),
                                 LSKRS,
                                 SFKRS,
                                 2)

Bose_symmetrized_AMP_P = Amplitude("K*(892)rho(770)_P",
                                 Variable("amp_real2", -0.02, 0.001, 0, 0),
                                 Variable("amp_imag2", -0.07, 0.001, 0, 0),
                                 LSKRP,
                                 SFKRP,
                                 2)
Bose_symmetrized_AMP_D = Amplitude("K*(892)rho(770)_D",
                                 Variable("amp_real3", -0.4, 0.001, 0, 0),
                                 Variable("amp_imag3", 0.4, 0.001, 0, 0),
                                 LSKRD,
                                 SFKRD,
                                 2)

Bose_symmetrized_KF = Amplitude("KF",
                              Variable("amp_real4", 0.0120787),
                              Variable("amp_imag4", -0.0332525),
                              LSKF,
                              SFKF,
                              2)

Bose_symmetrized_KK = Amplitude("LSKK",
                              Variable("amp_real5", 0.0109033),
                              Variable("amp_imag5", -0.00186219),
                              LSKK,
                              SFKK,
                              2)

Bose_symmetrized_K1R = Amplitude("LSK1R",
                               Variable("amp_real6", -0.10728),
                               Variable("amp_imag6", -0.130213),
                               LSK1R,
                               SFK1R,
                               2)

Bose_symmetrized_A1R = Amplitude("LSA1R",
                               Variable("amp_real7", 1.0),
                               Variable("amp_imag7", 0.0),
                               LSA1R,
                               SFA1R,
                               2)

Bose_symmetrized_A1RD = Amplitude("LSA1RD",
                                Variable("amp_real8", -0.94921),
                                Variable("amp_imag8", -1.73407),
                                LSA1RD,
                                SFA1RD,
                                2)

DK3P_DI.amplitudes = (Bose_symmetrized_KF,
                      Bose_symmetrized_AMP_S,
                      Bose_symmetrized_AMP_P,
                      Bose_symmetrized_AMP_D,
                      Bose_symmetrized_KK,
                      Bose_symmetrized_K1R,
                      Bose_symmetrized_A1R,
                      Bose_symmetrized_A1RD)

for res in LSKF:
    res.setParameterConstantness(True)

for res in LSKRS:
    res.setParameterConstantness(True)

for res in LSKRP:
    res.setParameterConstantness(True)

for res in LSKRD:
    res.setParameterConstantness(True)

for res in LSKK:
    res.setParameterConstantness(True)

for res in LSK1R:
    res.setParameterConstantness(True)

for res in LSA1R:
    res.setParameterConstantness(True)

for res in LSA1RD:
    res.setParameterConstantness(True)

constantOne = Variable("constantOne", 1)
constantZero = Variable("constantZero", 0)

observables  = (m12, m34, cos12, cos34, phi, eventNumber)
coefficients = (constantOne,)
offsets      = (constantZero, constantZero)

eff = PolynomialPdf("constantEff", observables, coefficients, offsets, 0)
dp = Amp4Body("test", observables, DK3P_DI, eff, 1000000)

constant = (Variable("constant", 0.1),)
constant2 = Variable("constant2", 1.0)
backgr = PolynomialPdf("backgr", m12, (constant))
signal = AddPdf("signal", constant2, dp, backgr)

signal.setData(currData)

dp.setDataSize(currData.getNumEvents(), 6)

datapdf = FitManager(signal)
datapdf.fit()

