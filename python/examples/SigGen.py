#!/usr/bin/env python
# -*- coding: utf-8 -*-

from goofit import *

# Constants used in more than one PDF component.
_mD0       = 1.8645
piPlusMass = 0.13957018
KmMass     =.493677

DK3P_DI = DecayInfo_DP()
DK3P_DI.meson_radius = 1.5
DK3P_DI.particle_masses = (_mD0,piPlusMass,piPlusMass,KmMass,piPlusMass)


RhoMass  = Variable("rho_mass", 0.77526, 0.01, 0.7, 0.8)
RhoWidth = Variable("rho_width", 0.1478, 0.01, 0.1, 0.2)
KstarM   = Variable("KstarM", 0.89581, 0.01, 0.9, 0.1)
KstarW   = Variable("KstarW", 0.0474, 0.01, 0.1, 0.2)
f600M    = Variable("f600M", 0.519, 0.01, 0.75, 0.85)
f600W    = Variable("f600W", 0.454, 0.01, 0.75, 0.85)
a1M      = Variable("a1M", 1.23, 0.01, 1.2, 1.3)
a1W      = Variable("a1W", 0.42, 0.01, 0.37, 0.47)
K1M      = Variable("K1M", 1.272, 0.01, 1.2, 1.3)
K1W      = Variable("K1W", 0.09, 0.01, 0.08, 0.1)
K1430M   = Variable("K1430M", 1.414, 0.01, 1.4, 1.5)
K1430W   = Variable("K1430W",.29, 0.01, 0.25, 0.35)

# Spin factors: we have two due to the bose symmetrization of the two pi+
SFKRS  = (SpinFactor("SF", SF_4Body.DtoV1V2_V1toP1P2_V2toP3P4_S, 0, 1, 2, 3),
SpinFactor("SF", SF_4Body.DtoV1V2_V1toP1P2_V2toP3P4_S, 3, 1, 2, 0))

SFKRP  = (SpinFactor("SF", SF_4Body.DtoV1V2_V1toP1P2_V2toP3P4_P, 0, 1, 2, 3),
SpinFactor("SF", SF_4Body.DtoV1V2_V1toP1P2_V2toP3P4_P, 3, 1, 2, 0))

SFKRD  = (SpinFactor("SF", SF_4Body.DtoV1V2_V1toP1P2_V2toP3P4_D, 0, 1, 2, 3),
SpinFactor("SF", SF_4Body.DtoV1V2_V1toP1P2_V2toP3P4_D, 3, 1, 2, 0))

SFKF   = (SpinFactor("SF", SF_4Body.DtoVS_VtoP1P2_StoP3P4, 2, 3, 0, 1),
SpinFactor("SF", SF_4Body.DtoVS_VtoP1P2_StoP3P4, 2, 0, 3, 1))

SFKK   = (SpinFactor("SF", SF_4Body.DtoAP1_AtoSP2_StoP3P4, 0, 1, 3, 2),
SpinFactor("SF", SF_4Body.DtoAP1_AtoSP2_StoP3P4, 3, 1, 0, 2))

SFK1R  = (SpinFactor("SF", SF_4Body.DtoAP1_AtoVP2_VtoP3P4, 3, 2, 0, 1),
SpinFactor("SF", SF_4Body.DtoAP1_AtoVP2_VtoP3P4, 0, 2, 3, 1))

SFA1R  = (SpinFactor("SF", SF_4Body.DtoAP1_AtoVP2_VtoP3P4, 2, 3, 0, 1),
SpinFactor("SF", SF_4Body.DtoAP1_AtoVP2_VtoP3P4, 2, 0, 3, 1))

SFA1RD = (SpinFactor("SF", SF_4Body.DtoAP1_AtoVP2Dwave_VtoP3P4, 2, 3, 0, 1),
SpinFactor("SF", SF_4Body.DtoAP1_AtoVP2Dwave_VtoP3P4, 2, 0, 3, 1))


# Lineshapes, also for both pi+ configurations
LSKRS = (Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_12, LS.BW, FF.BL2),Lineshape("K*(892)bar", KstarM, KstarW, 1, M_34, LS.BW, FF.BL2),Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_24, LS.BW, FF.BL2),Lineshape("K*(892)bar", KstarM, KstarW, 1, M_13, LS.BW, FF.BL2))
LSKRP = (Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_12, LS.BW, FF.BL2),Lineshape("K*(892)bar", KstarM, KstarW, 1, M_34, LS.BW, FF.BL2),Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_24, LS.BW, FF.BL2),Lineshape("K*(892)bar", KstarM, KstarW, 1, M_13, LS.BW, FF.BL2))
LSKRD = (Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_12, LS.BW, FF.BL2),Lineshape("K*(892)bar", KstarM, KstarW, 1, M_34, LS.BW, FF.BL2),Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_24, LS.BW, FF.BL2),Lineshape("K*(892)bar", KstarM, KstarW, 1, M_13, LS.BW, FF.BL2))
LSKF = (Lineshape("K*(892)bar", KstarM, KstarW, 1, M_34, LS.BW),Lineshape("f600", f600M, f600W, 0, M_12, LS.Bugg),Lineshape("K*(892)bar", KstarM, KstarW, 1, M_13, LS.BW),Lineshape("f600", f600M, f600W, 0, M_24, LS.Bugg))
LSKK = (Lineshape("K(1)(1270)bar", K1M, K1W, 1, M_34_2, LS.SBW),Lineshape("K(0)*(1430)bar", K1430M, K1430W, 0, M_34, LS.Lass),Lineshape("K(1)(1270)bar2", K1M, K1W, 1, M_13_2, LS.SBW),Lineshape("K(0)*(1430)bar2", K1430M, K1430W, 0, M_13, LS.Lass))
LSK1R = (Lineshape("K(1)(1270)bar", K1M, K1W, 0, M_12_3, LS.SBW),Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_12, LS.BW),Lineshape("K(1)(1270)bar", K1M, K1W, 0, M_24_3, LS.SBW),Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_24, LS.BW))
LSA1R = (Lineshape("a(1)(1260)+", a1M, a1W, 0, M_12_4, LS.SBW),Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_12, LS.BW),Lineshape("a(1)(1260)+", a1M, a1W, 0, M_24_1, LS.SBW),Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_24, LS.BW))
LSA1RD = (Lineshape("a(1)(1260)+", a1M, a1W, 2, M_12_4, LS.SBW),Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_12, LS.BW),Lineshape("a(1)(1260)+", a1M, a1W, 2, M_24_1, LS.SBW),Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_24, LS.BW))

# the very last parameter means that we have two permutations. so the first half of the Lineshapes
# and the first half of the spinfactors are amplitude 1, rest is amplitude 2
# This means that it is important for symmetrized amplitueds that the spinfactors and lineshapes are in the "right"
# order
Bose_symmetrized_AMP_S = Amplitude("K*(892)rho(770)_S",
Variable("amp_real1", -0.115177),
Variable("amp_imag1", 0.153976),
LSKRS,
SFKRS,
2)
Bose_symmetrized_AMP_P = Amplitude("K*(892)rho(770)_P",
Variable("amp_real2", -0.0298697),
Variable("amp_imag2", -0.0722874),
LSKRP,
SFKRP,
2)
Bose_symmetrized_AMP_D = Amplitude("K*(892)rho(770)_D",
Variable("amp_real3", -0.452212),
Variable("amp_imag3", 0.426521),
LSKRD,
SFKRD,
2)

Bose_symmetrized_KF = Amplitude(
"KF", Variable("amp_real3", 0.0120787), Variable("amp_imag3", -0.0332525), LSKF, SFKF, 2)
Bose_symmetrized_KK = Amplitude(
"LSKK", Variable("amp_real3", 0.0109033), Variable("amp_imag3", -0.00186219), LSKK, SFKK, 2)
Bose_symmetrized_K1R = Amplitude(
"LSK1R", Variable("amp_real3", -0.10728), Variable("amp_imag3", -0.130213), LSK1R, SFK1R, 2)
Bose_symmetrized_A1R = Amplitude(
"LSA1R", Variable("amp_real3", 1.0), Variable("amp_imag3", 0.0), LSA1R, SFA1R, 2)
Bose_symmetrized_A1RD = Amplitude(
"LSA1RD", Variable("amp_real3", -0.94921), Variable("amp_imag3", -1.73407), LSA1RD, SFA1RD, 2)

DK3P_DI.amplitudes = (Bose_symmetrized_KF,Bose_symmetrized_AMP_S,Bose_symmetrized_AMP_P,Bose_symmetrized_AMP_D,Bose_symmetrized_KK,Bose_symmetrized_K1R,Bose_symmetrized_A1R,Bose_symmetrized_A1RD)

m12                 = Variable("m12", 0, 3)
m34                 = Variable("m34", 0, 3)
cos12               = Variable("cos12", -1, 1)
cos34               = Variable("m12", -1, 1)
phi                 = Variable("phi", -3.5, 3.5)
eventNumber         = CountingVariable("eventNumber", 0, INT_MAX)
constantOne         = Variable("constantOne", 1)
constantZero        = Variable("constantZero", 0)


observables = (m12,m34,cos12,cos34,phi,eventNumber)
offsets = (constantZero, constantZero,)
coefficients = (constantOne,)

eff = PolynomialPdf("constantEff", observables, coefficients, offsets, 0)
dp  = DPPdf("test", observables, DK3P_DI, eff, 5)

#file = TFile("SigGen.root", "RECREATE")
#tree = TTree("events", "events")


'''
tree.Branch("m12", tm12, "m12/D")
tree.Branch("m34", tm34, "m34/D")
tree.Branch("c12", tc12, "c12/D")
tree.Branch("c34", tc34, "c34/D")
tree.Branch("phi", tphi, "phi/D")
'''
k=0
while k < 4:
    numEvents = 1000000
    dp.setGenerationOffset(k * numEvents)
    tuple = dp.GenerateSig(numEvents)
    k+=1


#tree.Write()
#file.Close()
#return 0
