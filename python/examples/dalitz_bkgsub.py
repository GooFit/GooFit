#!/usr/bin/env python

import sys

import matplotlib
import pandas as pd
import numpy as np
from copy import deepcopy

from goofit import *

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# This is really unnecessary...
from dalitz import getToyData

print_goofit_info()

_mD0 = 1.86484
_mD02 = _mD0 * _mD0
_mD02inv = 1.0 / _mD02
piPlusMass = 0.13957018
piZeroMass = 0.1349766

def make_sig_pdf(m12, m13, eventNumber, eff=None, fitMasses=False, fixAmps=False):
    # Constants used in more than one PDF component.
    # motherM = Variable("motherM", _mD0)
    # chargeM = Variable("chargeM", piPlusMass)
    # neutrlM = Variable("neutrlM", piZeroMass)
    # massSum = Variable("massSum", _mD0 *_mD0 + 2 * piPlusMass * piPlusMass + piZeroMass * piZeroMass) # = 3.53481

    constantOne = Variable("constantOne", 1)
    constantZero = Variable("constantZero", 0)

    fixedRhoMass = Variable("rho_mass", 0.7758, 0.01, 0.7, 0.8)
    fixedRhoWidth = Variable("rho_width", 0.1503, 0.01, 0.1, 0.2)

    dtop0pp = DecayInfo3()
    dtop0pp.motherMass = _mD0
    dtop0pp.daug1Mass = piZeroMass
    dtop0pp.daug2Mass = piPlusMass
    dtop0pp.daug3Mass = piPlusMass
    dtop0pp.meson_radius = 1.5

    rhop = Resonances.RBW(
        "rhop",
        Variable("rhop_amp_real", 1),
        Variable("rhop_amp_imag", 0),
        fixedRhoMass,
        fixedRhoWidth,
        1,
        PAIR_12,
    )

    rhom = Resonances.RBW(
        "rhom",
        Variable("rhom_amp_real", 0.714)
        if fixAmps
        else Variable("rhom_amp_real", 0.714, 0.001, 0, 0),
        Variable("rhom_amp_imag", -0.025)
        if fixAmps
        else Variable("rhom_amp_imag", -0.025, 0.1, 0, 0),
        fixedRhoMass,
        fixedRhoWidth,
        1,
        PAIR_13,
    )

    rho0 = Resonances.RBW(
        "rho0",
        Variable("rho0_amp_real", 0.565)
        if fixAmps
        else Variable("rho0_amp_real", 0.565, 0.001, 0, 0),
        Variable("rho0_amp_imag", 0.164)
        if fixAmps
        else Variable("rho0_amp_imag", 0.164, 0.1, 0, 0),
        fixedRhoMass,
        fixedRhoWidth,
        1,
        PAIR_23,
    )

    sharedMass = Variable("rhop_1450_mass", 1.465, 0.01, 1.0, 2.0)
    shareWidth = Variable("rhop_1450_width", 0.400, 0.01, 0.01, 5.0)

    rhop_1450 = Resonances.RBW(
        "rhop_1450",
        Variable("rhop_1450_amp_real", -0.174)
        if fixAmps
        else Variable("rhop_1450_amp_real", -0.174, 0.001, 0, 0),
        Variable("rhop_1450_amp_imag", -0.117)
        if fixAmps
        else Variable("rhop_1450_amp_imag", -0.117, 0.1, 0, 0),
        sharedMass,
        shareWidth,
        1,
        PAIR_12,
    )

    rho0_1450 = Resonances.RBW(
        "rho0_1450",
        Variable("rho0_1450_amp_real", 0.325)
        if fixAmps
        else Variable("rho0_1450_amp_real", 0.325, 0.001, 0, 0),
        Variable("rho0_1450_amp_imag", 0.057)
        if fixAmps
        else Variable("rho0_1450_amp_imag", 0.057, 0.1, 0, 0),
        sharedMass,
        shareWidth,
        1,
        PAIR_23,
    )

    rhom_1450 = Resonances.RBW(
        "rhom_1450",
        Variable("rhom_1450_amp_real", 0.788)
        if fixAmps
        else Variable("rhom_1450_amp_real", 0.788, 0.001, 0, 0),
        Variable("rhom_1450_amp_imag", 0.226)
        if fixAmps
        else Variable("rhom_1450_amp_imag", 0.226, 0.1, 0, 0),
        sharedMass,
        shareWidth,
        1,
        PAIR_13,
    )

    sharedMass = Variable("rhop_1700_mass", 1.720, 0.01, 1.6, 1.9)
    shareWidth = Variable("rhop_1700_width", 0.250, 0.01, 0.1, 1.0)

    rhop_1700 = Resonances.RBW(
        "rhop_1700",
        Variable("rhop_1700_amp_real", 2.151)
        if fixAmps
        else Variable("rhop_1700_amp_real", 2.151, 0.001, 0, 0),
        Variable("rhop_1700_amp_imag", -0.658)
        if fixAmps
        else Variable("rhop_1700_amp_imag", -0.658, 0.1, 0, 0),
        sharedMass,
        shareWidth,
        1,
        PAIR_12,
    )

    rho0_1700 = Resonances.RBW(
        "rho0_1700",
        Variable("rho0_1700_amp_real", 2.400)
        if fixAmps
        else Variable("rho0_1700_amp_real", 2.400, 0.001, 0, 0),
        Variable("rho0_1700_amp_imag", -0.734)
        if fixAmps
        else Variable("rho0_1700_amp_imag", -0.734, 0.1, 0, 0),
        sharedMass,
        shareWidth,
        1,
        PAIR_23,
    )

    rhom_1700 = Resonances.RBW(
        "rhom_1700",
        Variable("rhom_1700_amp_real", 1.286)
        if fixAmps
        else Variable("rhom_1700_amp_real", 1.286, 0.001, 0, 0),
        Variable("rhom_1700_amp_imag", -1.532)
        if fixAmps
        else Variable("rhom_1700_amp_imag", -1.532, 0.1, 0, 0),
        sharedMass,
        shareWidth,
        1,
        PAIR_13,
    )

    f0_980 = Resonances.RBW(
        "f0_980",
        Variable("f0_980_amp_real", 0.008 * (-_mD02))
        if fixAmps
        else Variable("f0_980_amp_real", 0.008 * (-_mD02), 0.001, 0, 0),
        Variable("f0_980_amp_imag", -0.013 * (-_mD02))
        if fixAmps
        else Variable("f0_980_amp_imag", -0.013 * (-_mD02), 0.1, 0, 0),
        Variable("f0_980_mass", 0.980, 0.01, 0.8, 1.2),
        Variable("f0_980_width", 0.044, 0.001, 0.001, 0.08),
        0,
        PAIR_23,
    )

    f0_1370 = Resonances.RBW(
        "f0_1370",
        Variable("f0_1370_amp_real", -0.058 * (-_mD02))
        if fixAmps
        else Variable("f0_1370_amp_real", -0.058 * (-_mD02), 0.001, 0, 0),
        Variable("f0_1370_amp_imag", 0.026 * (-_mD02))
        if fixAmps
        else Variable("f0_1370_amp_imag", 0.026 * (-_mD02), 0.1, 0, 0),
        Variable("f0_1370_mass", 1.434, 0.01, 1.2, 1.6),
        Variable("f0_1370_width", 0.173, 0.01, 0.01, 0.4),
        0,
        PAIR_23,
    )

    f0_1500 = Resonances.RBW(
        "f0_1500",
        Variable("f0_1500_amp_real", 0.057 * (-_mD02))
        if fixAmps
        else Variable("f0_1500_amp_real", 0.057 * (-_mD02), 0.001, 0, 0),
        Variable("f0_1500_amp_imag", 0.012 * (-_mD02))
        if fixAmps
        else Variable("f0_1500_amp_imag", 0.012 * (-_mD02), 0.1, 0, 0),
        Variable("f0_1500_mass", 1.507, 0.01, 1.3, 1.7),
        Variable("f0_1500_width", 0.109, 0.01, 0.01, 0.3),
        0,
        PAIR_23,
    )

    f0_1710 = Resonances.RBW(
        "f0_1710",
        Variable("f0_1710_amp_real", 0.070 * (-_mD02))
        if fixAmps
        else Variable("f0_1710_amp_real", 0.070 * (-_mD02), 0.001, 0, 0),
        Variable("f0_1710_amp_imag", 0.087 * (-_mD02))
        if fixAmps
        else Variable("f0_1710_amp_imag", 0.087 * (-_mD02), 0.1, 0, 0),
        Variable("f0_1710_mass", 1.714, 0.01, 1.5, 2.9),
        Variable("f0_1710_width", 0.140, 0.01, 0.01, 0.5),
        0,
        PAIR_23,
    )

    f2_1270 = Resonances.RBW(
        "f2_1270",
        Variable("f2_1270_amp_real", -1.027 * (-_mD02inv))
        if fixAmps
        else Variable("f2_1270_amp_real", -1.027 * (-_mD02inv), 0.001, 0, 0),
        Variable("f2_1270_amp_imag", -0.162 * (-_mD02inv))
        if fixAmps
        else Variable("f2_1270_amp_imag", -0.162 * (-_mD02inv), 0.1, 0, 0),
        Variable("f2_1270_mass", 1.2754, 0.01, 1.0, 1.5),
        Variable("f2_1270_width", 0.1851, 0.01, 0.01, 0.4),
        2,
        PAIR_23,
    )

    f0_600 = Resonances.RBW(
        "f0_600",
        Variable("f0_600_amp_real", 0.068 * (-_mD02))
        if fixAmps
        else Variable("f0_600_amp_real", 0.068 * (-_mD02), 0.001, 0, 0),
        Variable("f0_600_amp_imag", 0.010 * (-_mD02))
        if fixAmps
        else Variable("f0_600_amp_imag", 0.010 * (-_mD02), 0.1, 0, 0),
        Variable("f0_600_mass", 0.500, 0.01, 0.3, 0.7),
        Variable("f0_600_width", 0.400, 0.01, 0.2, 0.6),
        0,
        PAIR_23,
    )

    dtop0pp.resonances = (
        rhop,
        rho0,
        rhom,
        rhop_1450,
        rho0_1450,
        rhom_1450,
        rhop_1700,
        rho0_1700,
        rhom_1700,
        f0_980,
        f0_1370,
        f0_1500,
        f0_1710,
        f2_1270,
        f0_600,
    )

    if not fitMasses:
        for res in dtop0pp.resonances:
            res.setParameterConstantness(True)

    if not eff:
        # By default create a constant efficiency.
        observables = (
            m12,
            m13,
        )
        offsets = (
            constantZero,
            constantZero,
        )
        coefficients = (constantOne,)
        eff = PolynomialPdf("constantEff", observables, coefficients, offsets, 0)

    return Amp3Body("signalPDF", m12, m13, eventNumber, dtop0pp, eff)


def make_bkg_pdf(m12, m13, eventNumber):
    bkg_di = DecayInfo3()
    bkg_di.motherMass = _mD0
    bkg_di.daug1Mass = piZeroMass
    bkg_di.daug2Mass = piPlusMass
    bkg_di.daug3Mass = piPlusMass
    bkg_di.meson_radius = 1.5

    # Define the amplitude (just phasespace).
    nonr = Resonances.NonRes("bkg_nonr", Variable("nrConstantOne", 1), Variable("nrConstantZero", 0))
    bkg_di.resonances = (nonr,)

    # Define the efficiency (constant).
    observables = (m12, m13)
    offsets = (Variable("constantZero1", 0), Variable("constantZero2", 0))
    coefficients = (Variable("constantOne1", 1),)
    eff = PolynomialPdf("constantEff", observables, coefficients, offsets, 0)

    return Amp3Body("bkgPDF", m12, m13, eventNumber, bkg_di, eff)


def main():
    filename = GOOFIT_SOURCE_DIR + "/examples/dalitz/dalitz_toyMC_000.txt"
    m12, m13, eventNumber, signal_data = getToyData(filename)

    signal_pdf = make_sig_pdf(m12, m13, eventNumber, fixAmps=True)
    bkg_pdf = make_bkg_pdf(m12, m13, eventNumber)
    bkg_fit_pdf = PolynomialPdf("bkgFitPDF", (m12, m13), (Variable("one", 1),), (Variable("zero1", 0), Variable("zero2", 0)), 0)

    nSig = 1000000
    nBkg = 25000
    print('Generating signal...')
    sig_particles, sig_variables, sig_weights, sig_flags = signal_pdf.GenerateSig(nSig)
    signal_pdf.setGenerationOffset(nSig)
    acc_sig = int(np.sum(sig_flags))
    sig_data = np.vstack([sig_variables[0][sig_flags], sig_variables[2][sig_flags], np.arange(acc_sig, dtype='uint32')])
    print(f'Accepted {acc_sig} out of {nSig} events')

    # print('Generating background...')
    # bkg_pdf.setGenerationOffset(nSig)
    # bkg_particles, bkg_variables, bkg_weights, bkg_flags = bkg_pdf.GenerateSig(nBkg)
    # bkg_pdf.setGenerationOffset(nSig + nBkg)
    # acc_bkg = int(np.sum(bkg_flags))
    # bkg_data = np.vstack([bkg_variables[0][bkg_flags], bkg_variables[2][bkg_flags], np.arange(acc_bkg, dtype='uint32') + acc_sig])
    # print(f'Accepted {acc_bkg} out of {nBkg} events')

    print('Creating dataset with signal and background...')
    data = UnbinnedDataSet(m12, m13, eventNumber)
    data.from_matrix(sig_data)
    # For now, fit a sample with 0 background.
    # data.from_matrix(bkg_data)

    print('Creating total PDF...')
    # signal_fraction = Variable("signal_fraction", 1.*acc_sig/(acc_sig + acc_bkg), 0., 1.)
    signal_fraction = Variable("signal_fraction", 1., 0., 1.)
    total_pdf = AddPdf("totalPDF", signal_fraction, signal_pdf, bkg_pdf)

    print('Fitting the total PDF...')
    total_pdf.setData(data)
    fitman = FitManager(total_pdf)
    fitman.fit()

if __name__ == '__main__':
    sys.exit(main())