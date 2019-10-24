#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

from goofit import *
import pandas as pd
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

print_goofit_info()

# Basic constants
_mD0       = 1.86484
_mD02      = _mD0 * _mD0
_mD02inv   = 1. / _mD02
piPlusMass = 0.13957018
piZeroMass = 0.1349766


def cpuGetM23(massPZ,massPM):
    print("cpuGetM23")
    return (_mD02 + piZeroMass * piZeroMass + piPlusMass * piPlusMass + piPlusMass * piPlusMass - massPZ - massPM)


def getToyData(toyFileName):
    print("Reading {} ...".format(toyFileName))
    out = pd.read_csv(toyFileName,
            delim_whitespace=True,
            comment='=',
	    skip_blank_lines=True,
            names='event m23 m12 m13 e23 e12 e13 tau sigma_t mD0 deltaM probSig dst_charge run event2 sig frac1 frac2 frac3 frac4'.split())
    out['vals'] = out.index
    print("Read in", len(out.index), "lines of data")

    m12 = Observable("m12", 0, 3)
    m13 = Observable("m13", 0, 3)
    m12.numbins = 240
    m13.numbins = 240
    eventNumber = EventNumber("eventNumber")

    data = UnbinnedDataSet(m12, m13, eventNumber)
    data.from_matrix(out[['m12','m13','vals']].T)

    # Add figure/subplot
    fig, ax = plt.subplots()
    ax.hist2d(out.m12, out.m13, bins=[100,100])
    fig.savefig("dalitz_data_plot.png")
    print("Original data plot:", "dalitz_data_plot.png")

    return m12, m13, eventNumber, data

# Not used, for example
def makeKzeroVeto(m12, m13, motherM, neutrlM, chargeM):
    kVetoInfo = VetoInfo()
    kVetoInfo.cyclic_index = PAIR_23
    kVetoInfo.minimum      = Variable("veto_min", 0.475 * 0.475)
    kVetoInfo.maximum      = Variable("veto_max", 0.505 * 0.505)
    vetos = (kVetoInfo)
    kzero_veto = DalitzVetoPdf("kzero_veto", m12, m13, motherM, neutrlM, chargeM, chargeM, vetos)
    return kzero_veto


def makeSignalPdf(m12, m13, eventNumber, eff = None, fitMasses = False):
    # Constants used in more than one PDF component.
    # motherM = Variable("motherM", _mD0)
    # chargeM = Variable("chargeM", piPlusMass)
    # neutrlM = Variable("neutrlM", piZeroMass)
    # massSum = Variable("massSum", _mD0 *_mD0 + 2 * piPlusMass * piPlusMass + piZeroMass * piZeroMass) # = 3.53481

    constantOne  = Variable("constantOne", 1)
    constantZero = Variable("constantZero", 0)

    fixedRhoMass = Variable("rho_mass", 0.7758, 0.01, 0.7, 0.8)
    fixedRhoWidth = Variable("rho_width", 0.1503, 0.01, 0.1, 0.2)

    dtop0pp = DecayInfo3()
    dtop0pp.motherMass   = _mD0
    dtop0pp.daug1Mass    = piZeroMass
    dtop0pp.daug2Mass    = piPlusMass
    dtop0pp.daug3Mass    = piPlusMass
    dtop0pp.meson_radius = 1.5

    rhop = Resonances.RBW("rhop",
                        Variable("rhop_amp_real", 1),
                        Variable("rhop_amp_imag", 0),
                        fixedRhoMass,
                        fixedRhoWidth,
                        1,
                        PAIR_12)

    fixAmps = False

    rhom = Resonances.RBW("rhom",
                        Variable("rhom_amp_real", 0.714) if fixAmps else Variable("rhom_amp_real", 0.714, 0.001, 0, 0),
                        Variable("rhom_amp_imag", -0.025) if fixAmps else Variable("rhom_amp_imag", -0.025, 0.1, 0, 0),
                        fixedRhoMass,
                        fixedRhoWidth,
                        1,
                        PAIR_13)

    rho0 = Resonances.RBW("rho0",
                        Variable("rho0_amp_real", 0.565) if fixAmps else Variable("rho0_amp_real", 0.565, 0.001, 0, 0),
                        Variable("rho0_amp_imag", 0.164) if fixAmps else Variable("rho0_amp_imag", 0.164, 0.1, 0, 0),
                        fixedRhoMass,
                        fixedRhoWidth,
                        1,
                        PAIR_23)

    sharedMass = Variable("rhop_1450_mass", 1.465, 0.01, 1.0, 2.0)
    shareWidth = Variable("rhop_1450_width", 0.400, 0.01, 0.01, 5.0)

    rhop_1450 = Resonances.RBW("rhop_1450",
                             Variable("rhop_1450_amp_real", -0.174) if fixAmps else Variable("rhop_1450_amp_real", -0.174, 0.001, 0, 0),
                             Variable("rhop_1450_amp_imag", -0.117) if fixAmps else Variable("rhop_1450_amp_imag", -0.117, 0.1, 0, 0),
                             sharedMass,
                             shareWidth,
                             1,
                             PAIR_12)

    rho0_1450 = Resonances.RBW("rho0_1450",
                             Variable("rho0_1450_amp_real", 0.325) if fixAmps else Variable("rho0_1450_amp_real", 0.325, 0.001, 0, 0),
                             Variable("rho0_1450_amp_imag", 0.057) if fixAmps else Variable("rho0_1450_amp_imag", 0.057, 0.1, 0, 0),
                             sharedMass,
                             shareWidth,
                             1,
                             PAIR_23)

    rhom_1450 = Resonances.RBW("rhom_1450",
                             Variable("rhom_1450_amp_real", 0.788) if fixAmps else Variable("rhom_1450_amp_real", 0.788, 0.001, 0, 0),
                             Variable("rhom_1450_amp_imag", 0.226) if fixAmps else Variable("rhom_1450_amp_imag", 0.226, 0.1, 0, 0),
                             sharedMass,
                             shareWidth,
                             1,
                             PAIR_13)

    sharedMass = Variable("rhop_1700_mass", 1.720, 0.01, 1.6, 1.9)
    shareWidth = Variable("rhop_1700_width", 0.250, 0.01, 0.1, 1.0)

    rhop_1700 = Resonances.RBW("rhop_1700",
                             Variable("rhop_1700_amp_real", 2.151) if fixAmps else Variable("rhop_1700_amp_real", 2.151, 0.001, 0, 0),
                             Variable("rhop_1700_amp_imag", -0.658) if fixAmps else Variable("rhop_1700_amp_imag", -0.658, 0.1, 0, 0),
                             sharedMass,
                             shareWidth,
                             1,
                             PAIR_12)

    rho0_1700 = Resonances.RBW("rho0_1700",
                             Variable("rho0_1700_amp_real", 2.400) if fixAmps else Variable("rho0_1700_amp_real", 2.400, 0.001, 0, 0),
                             Variable("rho0_1700_amp_imag", -0.734) if fixAmps else Variable("rho0_1700_amp_imag", -0.734, 0.1, 0, 0),
                             sharedMass,
                             shareWidth,
                             1,
                             PAIR_23)

    rhom_1700 = Resonances.RBW("rhom_1700",
                             Variable("rhom_1700_amp_real", 1.286) if fixAmps else Variable("rhom_1700_amp_real", 1.286, 0.001, 0, 0),
                             Variable("rhom_1700_amp_imag", -1.532) if fixAmps else Variable("rhom_1700_amp_imag", -1.532, 0.1, 0, 0),
                             sharedMass,
                             shareWidth,
                             1,
                             PAIR_13)

    f0_980    = Resonances.RBW("f0_980",
                            Variable("f0_980_amp_real", 0.008 * (-_mD02))if fixAmps else Variable("f0_980_amp_real", 0.008 * (-_mD02), 0.001, 0, 0),
                            Variable("f0_980_amp_imag", -0.013 * (-_mD02))if fixAmps else Variable("f0_980_amp_imag", -0.013 * (-_mD02), 0.1, 0, 0),
                            Variable("f0_980_mass", 0.980, 0.01, 0.8, 1.2),
                            Variable("f0_980_width", 0.044, 0.001, 0.001, 0.08),
                            0,
                            PAIR_23)

    f0_1370   = Resonances.RBW("f0_1370",
                            Variable("f0_1370_amp_real", -0.058 * (-_mD02))if fixAmps else Variable("f0_1370_amp_real", -0.058 * (-_mD02), 0.001, 0, 0),
                            Variable("f0_1370_amp_imag", 0.026 * (-_mD02))if fixAmps else Variable("f0_1370_amp_imag", 0.026 * (-_mD02), 0.1, 0, 0),
                            Variable("f0_1370_mass", 1.434, 0.01, 1.2, 1.6),
                            Variable("f0_1370_width", 0.173, 0.01, 0.01, 0.4),
                            0,
                            PAIR_23)

    f0_1500   = Resonances.RBW("f0_1500",
                            Variable("f0_1500_amp_real", 0.057 * (-_mD02))if fixAmps else Variable("f0_1500_amp_real", 0.057 * (-_mD02), 0.001, 0, 0),
                            Variable("f0_1500_amp_imag", 0.012 * (-_mD02))if fixAmps else Variable("f0_1500_amp_imag", 0.012 * (-_mD02), 0.1, 0, 0),
                            Variable("f0_1500_mass", 1.507, 0.01, 1.3, 1.7),
                            Variable("f0_1500_width", 0.109, 0.01, 0.01, 0.3),
                            0,
                            PAIR_23)

    f0_1710   = Resonances.RBW("f0_1710",
                            Variable("f0_1710_amp_real", 0.070 * (-_mD02))if fixAmps else Variable("f0_1710_amp_real", 0.070 * (-_mD02), 0.001, 0, 0),
                            Variable("f0_1710_amp_imag", 0.087 * (-_mD02))if fixAmps else Variable("f0_1710_amp_imag", 0.087 * (-_mD02), 0.1, 0, 0),
                            Variable("f0_1710_mass", 1.714, 0.01, 1.5, 2.9),
                            Variable("f0_1710_width", 0.140, 0.01, 0.01, 0.5),
                            0,
                            PAIR_23)

    f2_1270 = Resonances.RBW("f2_1270",
                        Variable("f2_1270_amp_real", -1.027 * (-_mD02inv)) if fixAmps else Variable("f2_1270_amp_real", -1.027 * (-_mD02inv), 0.001, 0, 0),
                        Variable("f2_1270_amp_imag", -0.162 * (-_mD02inv))if fixAmps else Variable("f2_1270_amp_imag", -0.162 * (-_mD02inv), 0.1, 0, 0),
                        Variable("f2_1270_mass", 1.2754, 0.01, 1.0, 1.5),
                        Variable("f2_1270_width", 0.1851, 0.01, 0.01, 0.4),
                        2,
                        PAIR_23)

    f0_600 = Resonances.RBW("f0_600",
                        Variable("f0_600_amp_real", 0.068 * (-_mD02))if fixAmps else Variable("f0_600_amp_real", 0.068 * (-_mD02), 0.001, 0, 0),
                        Variable("f0_600_amp_imag", 0.010 * (-_mD02))if fixAmps else Variable("f0_600_amp_imag", 0.010 * (-_mD02), 0.1, 0, 0),
                        Variable("f0_600_mass", 0.500, 0.01, 0.3, 0.7),
                        Variable("f0_600_width", 0.400, 0.01, 0.2, 0.6),
                        0,
                        PAIR_23)

    nonr = Resonances.NonRes("nonr",
                        Variable("nonr_amp_real", 0.5595 * (-1))if fixAmps else Variable("nonr_amp_real", 0.5595 * (-1), 0.001, 0, 0),
                        Variable("nonr_amp_imag", -0.108761 * (-1))if fixAmps else Variable("nonr_amp_imag", -0.108761 * (-1), 0.1, 0, 0))

    dtop0pp.resonances = (nonr,rhop,rho0,rhom,rhop_1450,rho0_1450,rhom_1450,rhop_1700,rho0_1700,rhom_1700,f0_980,f0_1370,f0_1500,f0_1710,f2_1270,f0_600)

    if not fitMasses:
        for res in dtop0pp.resonances:
            res.setParameterConstantness(True)

    if not eff:
        # By default create a constant efficiency.
        observables = (m12,m13,)
        offsets = (constantZero,constantZero,)
        coefficients = (constantOne,)
        eff = PolynomialPdf("constantEff", observables, coefficients, offsets, 0)

    d = Amp3Body("signalPDF", m12, m13, eventNumber, dtop0pp, eff)
    return d

def runToyFit(toyFileName):
    m12, m13, eventNumber, data = getToyData(toyFileName)
    signal = makeSignalPdf(m12, m13, eventNumber)
    signal.setData(data)
    signal.setDataSize(data.getNumEvents())
    fitman = FitManager(signal)

    fitman.fit()
    return fitman, signal, m12, m13

def main():
    filename = GOOFIT_SOURCE_DIR + "/examples/dalitz/dalitz_toyMC_000.txt"
    fitman, signal, m12, m23 =  runToyFit(filename)

    # A wrapper to avoid segfaulting when accessing a complex component
    prodpdf = ProdPdf("prodpdf", [signal])

    # Add nice tool for making data or plotting
    dplotter = DalitzPlotter(prodpdf, signal)

    arr = dplotter.make2D()
    extent = dplotter.getExtent()
    plt.imshow(arr, extent=extent, origin='lower')

    plt.savefig("dalitz_pdf_plot.png")
    print("Fit PDF plot:", "dalitz_pdf_plot.png")

    # Double int's needed here since int may return a long and sys.exit needs a real int for Python2
    return int(int(fitman))

if __name__ == '__main__':
    sys.exit(main())
