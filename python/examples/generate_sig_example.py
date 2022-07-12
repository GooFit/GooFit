#!/usr/bin/env python


import matplotlib.pyplot as plt
import numpy as np

from goofit import *

print_goofit_info()

# Constants used in more than one PDF component.
_mD0 = 1.8645
piPlusMass = 0.13957018
KmMass = 0.493677

# DK3P_DI = DecayInfo4()

DK3P_DI = DecayInfo4t(
    Variable("tau", 0.4101, 0.001, 0.300, 0.500),
    Variable("xmixing", 0.0),
    Variable("ymixing", 0.0),
    Variable("SqWStoRSrate", 1.0 / np.sqrt(300.0)),
)

DK3P_DI.meson_radius = 1.5
DK3P_DI.particle_masses = (_mD0, piPlusMass, piPlusMass, KmMass, piPlusMass)


RhoMass = Variable("rho_mass", 0.77526, 0.01, 0.7, 0.8)
RhoWidth = Variable("rho_width", 0.1478, 0.01, 0.1, 0.2)
KstarM = Variable("KstarM", 0.89581, 0.01, 0.9, 0.1)
KstarW = Variable("KstarW", 0.0474, 0.01, 0.1, 0.2)
f600M = Variable("f600M", 0.519, 0.01, 0.75, 0.85)
f600W = Variable("f600W", 0.454, 0.01, 0.75, 0.85)
a1M = Variable("a1M", 1.23, 0.01, 1.2, 1.3)
a1W = Variable("a1W", 0.42, 0.01, 0.37, 0.47)
K1M = Variable("K1M", 1.272, 0.01, 1.2, 1.3)
K1W = Variable("K1W", 0.09, 0.01, 0.08, 0.1)
K1430M = Variable("K1430M", 1.414, 0.01, 1.4, 1.5)
K1430W = Variable("K1430W", 0.29, 0.01, 0.25, 0.35)

# Spin factors: we have two due to the bose symmetrization of the two pi+
SFKRS = (
    SpinFactor("SF", SF_4Body.DtoV1V2_V1toP1P2_V2toP3P4_S, _mD0, 0, 1, 2, 3),
    SpinFactor("SF", SF_4Body.DtoV1V2_V1toP1P2_V2toP3P4_S, _mD0, 3, 1, 2, 0),
)

SFKRP = (
    SpinFactor("SF", SF_4Body.DtoV1V2_V1toP1P2_V2toP3P4_P, _mD0, 0, 1, 2, 3),
    SpinFactor("SF", SF_4Body.DtoV1V2_V1toP1P2_V2toP3P4_P, _mD0, 3, 1, 2, 0),
)

SFKRD = (
    SpinFactor("SF", SF_4Body.DtoV1V2_V1toP1P2_V2toP3P4_D, _mD0, 0, 1, 2, 3),
    SpinFactor("SF", SF_4Body.DtoV1V2_V1toP1P2_V2toP3P4_D, _mD0, 3, 1, 2, 0),
)

SFKF = (
    SpinFactor("SF", SF_4Body.DtoVS_VtoP1P2_StoP3P4, _mD0, 2, 3, 0, 1),
    SpinFactor("SF", SF_4Body.DtoVS_VtoP1P2_StoP3P4, _mD0, 2, 0, 3, 1),
)

SFKK = (
    SpinFactor("SF", SF_4Body.DtoAP1_AtoSP2_StoP3P4, _mD0, 0, 1, 3, 2),
    SpinFactor("SF", SF_4Body.DtoAP1_AtoSP2_StoP3P4, _mD0, 3, 1, 0, 2),
)

SFK1R = (
    SpinFactor("SF", SF_4Body.DtoAP1_AtoVP2_VtoP3P4, _mD0, 3, 2, 0, 1),
    SpinFactor("SF", SF_4Body.DtoAP1_AtoVP2_VtoP3P4, _mD0, 0, 2, 3, 1),
)

SFA1R = (
    SpinFactor("SF", SF_4Body.DtoAP1_AtoVP2_VtoP3P4, _mD0, 2, 3, 0, 1),
    SpinFactor("SF", SF_4Body.DtoAP1_AtoVP2_VtoP3P4, _mD0, 2, 0, 3, 1),
)

SFA1RD = (
    SpinFactor("SF", SF_4Body.DtoAP1_AtoVP2Dwave_VtoP3P4, _mD0, 2, 3, 0, 1),
    SpinFactor("SF", SF_4Body.DtoAP1_AtoVP2Dwave_VtoP3P4, _mD0, 2, 0, 3, 1),
)


# Lineshapes, also for both pi+ configurations
LSKRS = (
    Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_12),
    Lineshapes.RBW("K*(892)bar", KstarM, KstarW, 1, M_34),
    Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_24),
    Lineshapes.RBW("K*(892)bar", KstarM, KstarW, 1, M_13),
)

LSKRP = (
    Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_12),
    Lineshapes.RBW("K*(892)bar", KstarM, KstarW, 1, M_34),
    Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_24),
    Lineshapes.RBW("K*(892)bar", KstarM, KstarW, 1, M_13),
)

LSKRD = (
    Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_12),
    Lineshapes.RBW("K*(892)bar", KstarM, KstarW, 1, M_34),
    Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_24),
    Lineshapes.RBW("K*(892)bar", KstarM, KstarW, 1, M_13),
)

LSKF = (
    Lineshapes.RBW("K*(892)bar", KstarM, KstarW, 1, M_34),
    Lineshapes.Bugg("f600", f600M, f600W, 0, M_12),
    Lineshapes.RBW("K*(892)bar", KstarM, KstarW, 1, M_13),
    Lineshapes.Bugg("f600", f600M, f600W, 0, M_24),
)

LSKK = (
    Lineshapes.SBW("K(1)(1270)bar", K1M, K1W, 1, M_34_2),
    Lineshapes.LASS("K(0)*(1430)bar", K1430M, K1430W, 0, M_34),
    Lineshapes.SBW("K(1)(1270)bar2", K1M, K1W, 1, M_13_2),
    Lineshapes.LASS("K(0)*(1430)bar2", K1430M, K1430W, 0, M_13),
)

LSK1R = (
    Lineshapes.SBW("K(1)(1270)bar", K1M, K1W, 0, M_12_3),
    Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_12),
    Lineshapes.SBW("K(1)(1270)bar", K1M, K1W, 0, M_24_3),
    Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_24),
)

LSA1R = (
    Lineshapes.SBW("a(1)(1260)+", a1M, a1W, 0, M_12_4),
    Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_12),
    Lineshapes.SBW("a(1)(1260)+", a1M, a1W, 0, M_24_1),
    Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_24),
)

LSA1RD = (
    Lineshapes.SBW("a(1)(1260)+", a1M, a1W, 2, M_12_4),
    Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_12),
    Lineshapes.SBW("a(1)(1260)+", a1M, a1W, 2, M_24_1),
    Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_24),
)

# the very last parameter means that we have two permutations. so the first half of the Lineshapes
# and the first half of the spinfactors are amplitude 1, rest is amplitude 2
# This means that it is important for symmetrized amplitueds that the spinfactors and lineshapes are in the "right"
# order
Bose_symmetrized_AMP_S = Amplitude(
    "K*(892)rho(770)_S",
    Variable("amp_real1", -0.115177),
    Variable("amp_imag1", 0.153976),
    LSKRS,
    SFKRS,
    2,
)
Bose_symmetrized_AMP_P = Amplitude(
    "K*(892)rho(770)_P",
    Variable("amp_real2", -0.0298697),
    Variable("amp_imag2", -0.0722874),
    LSKRP,
    SFKRP,
    2,
)
Bose_symmetrized_AMP_D = Amplitude(
    "K*(892)rho(770)_D",
    Variable("amp_real3", -0.452212),
    Variable("amp_imag3", 0.426521),
    LSKRD,
    SFKRD,
    2,
)

Bose_symmetrized_KF = Amplitude(
    "KF",
    Variable("amp_real3", 0.0120787),
    Variable("amp_imag3", -0.0332525),
    LSKF,
    SFKF,
    2,
)

Bose_symmetrized_KK = Amplitude(
    "LSKK",
    Variable("amp_real3", 0.0109033),
    Variable("amp_imag3", -0.00186219),
    LSKK,
    SFKK,
    2,
)

Bose_symmetrized_K1R = Amplitude(
    "LSK1R",
    Variable("amp_real3", -0.10728),
    Variable("amp_imag3", -0.130213),
    LSK1R,
    SFK1R,
    2,
)

Bose_symmetrized_A1R = Amplitude(
    "LSA1R", Variable("amp_real3", 1.0), Variable("amp_imag3", 0.0), LSA1R, SFA1R, 2
)

Bose_symmetrized_A1RD = Amplitude(
    "LSA1RD",
    Variable("amp_real3", -0.94921),
    Variable("amp_imag3", -1.73407),
    LSA1RD,
    SFA1RD,
    2,
)

DK3P_DI.amplitudes = (
    Bose_symmetrized_KF,
    Bose_symmetrized_AMP_S,
    Bose_symmetrized_AMP_P,
    Bose_symmetrized_AMP_D,
    Bose_symmetrized_KK,
    Bose_symmetrized_K1R,
    Bose_symmetrized_A1R,
    Bose_symmetrized_A1RD,
)

m12 = Observable("m12", 0, 3)
m34 = Observable("m34", 0, 3)
cos12 = Observable("cos12", -1, 1)
cos34 = Observable("m12", -1, 1)
phi = Observable("phi", -3.5, 3.5)
eventNumber = EventNumber("eventNumber", 0, int(10e7))
dtime = Observable("dtime", 0, 100)
sigmat = Observable("sigmat", -3, 3)
constantOne = Variable("constantOne", 1)
constantZero = Variable("constantZero", 0)

# time independent observable order
# observables = (m12, m34, cos12, cos34, phi, eventNumber)
# time dependent observable order
observables = (m12, m34, cos12, cos34, phi, eventNumber, dtime, sigmat)
offsets = (
    constantZero,
    constantZero,
)
coefficients = (constantOne,)

eff = PolynomialPdf("constantEff", observables, coefficients, offsets, 0)
# dp = Amp4Body("test", observables, DK3P_DI, eff, 5)
mistag = None
res = TruthResolution()

dp = Amp4Body_TD("test", observables, DK3P_DI, res, eff, mistag, 0, 1)
accepted = []
accepted_m12 = []
accepted_m34 = []
accepted_c12 = []
accepted_c34 = []
accepted_phi = []
accepted_dtime = []

for k in range(4):
    numEvents = 50000
    dp.setGenerationOffset(k * numEvents)
    print("Using accept-reject method", end=" ")
    particles, variables, weights, flags = dp.GenerateSig(numEvents)
    # Flags == 1 means the data was accepted
    accepted_m12.append((variables[0])[flags])
    accepted_m34.append((variables[1])[flags])
    accepted_c12.append((variables[2])[flags])
    accepted_c34.append((variables[3])[flags])
    accepted_phi.append((variables[4])[flags])
    accepted_dtime.append((variables[4])[flags])
    print(accepted_m12[-1:-5])
    accepted = int(np.sum(flags))
    print("would leave you with", len(accepted_m12[0]), "out of", numEvents, "events")

print(accepted_m12[0])
# plt.figure()
# plt.hist(accepted_m12[0], bins=100)
# # plt.savefig('amp4body_amplitude_plots/siggen_example_m12.png')
# plt.savefig("amp4body_td_amplitude_plots/siggen_example_m12.png")


def plot_multi_body_resonance(amplitude, index):
    accepted_m12_amp = []
    accepted_m34_amp = []
    accepted_c12_amp = []
    accepted_c34_amp = []
    accepted_phi_amp = []
    accepted_dtime_amp = []
    # DK3P_DI_amp = DecayInfo4()
    # create a new Decay Info object

    DK3P_DI_amp = DecayInfo4t(
        Variable("tau", 0.4101, 0.001, 0.300, 0.500),
        Variable("xmixing", 0.0),
        Variable("ymixing", 0.0),
        Variable("SqWStoRSrate", 1.0 / np.sqrt(300.0)),
    )

    DK3P_DI_amp.meson_radius = 1.5
    DK3P_DI_amp.particle_masses = (_mD0, piPlusMass, piPlusMass, KmMass, piPlusMass)
    DK3P_DI_amp.amplitudes = [(amplitude)]
    # create an Amp4Body_TD with a single amplitude
    mistag = None
    res = TruthResolution()
    # dp_amp = Amp4Body("test_amp", observables, DK3P_DI_amp, eff, 5)
    dp_amp = Amp4Body_TD("test_amp", observables, DK3P_DI_amp, res, eff, mistag, 0, 1)
    print("Generating toys with a single amplitude")
    for k in range(4):
        numEvents = 50000
        dp_amp.setGenerationOffset(k * numEvents)
        print("Using accept-reject method", end=" ")
        particles, variables, weights, flags = dp_amp.GenerateSig(numEvents)
        accepted = int(np.sum(flags))
        print("would leave you with", accepted, "out of", numEvents, "events")
        accepted_m12_amp.append((variables[0])[flags])
        accepted_m34_amp.append((variables[1])[flags])
        accepted_c12_amp.append((variables[2])[flags])
        accepted_c34_amp.append((variables[3])[flags])
        accepted_phi_amp.append((variables[4])[flags])
        accepted_dtime_amp.append((variables[5])[flags])
    # plot output
    plt.figure()
    plt.hist(accepted_m12_amp, bins=100)
    # plt.savefig(f'amp4body_amplitude_plots/amplitude_{i}_accepted_m12.png')
    plt.savefig(f"amp4body_td_amplitude_plots/amplitude_{index}_accepted_m12.png")


for i in range(len(DK3P_DI.amplitudes)):
    plot_multi_body_resonance(DK3P_DI.amplitudes[i], i)
