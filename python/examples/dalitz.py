#not /usr/bin/env python
# -*- coding: utf-8 -*-

from goofit import *
import io

data = UnbinnedDataSet()
m12  = Variable("m12",0)
m13  = Variable("m13",0)
eventNumber = CountingVariable("eventNumber",0,0)
fitMasses = False
fixedRhoMass = Variable("rho_mass", 0.7758, 0.01, 0.7, 0.8)
fixedRhoWidth = Variable("rho_width", 0.1503, 0.01, 0.1, 0.2)

_mD0       = 1.86484
_mD02      = _mD0 * _mD0
_mD02inv   = 1. / _mD02
piPlusMass = 0.13957018
piZeroMass = 0.1349766

#Constants used in more than one PDF component.
motherM = Variable("motherM", _mD0)
chargeM = Variable("chargeM", piPlusMass)
neutrlM = Variable("neutrlM", piZeroMass)
massSum = Variable("massSum", _mD0 *_mD0 + 2 * piPlusMass * piPlusMass + piZeroMass * piZeroMass) # = 3.53481
constantOne  = Variable("constantOne", 1)
constantZero = Variable("constantZero", 0)

kzero_veto = 0

def cpuGetM23(massPZ,massPM):
    print("cpuGetM23")
    return (_mD02 + piZeroMass * piZeroMass + piPlusMass * piPlusMass + piPlusMass * piPlusMass - massPZ - massPM)


def getToyData(toyFileName):
    print("getToyData")
    data = UnbinnedDataSet(m12,m13,eventNumber)
    print("yay")
    buffer = ""
    print("yay0")
    with open(toyFileName) as openfileobject:
        print("yay1")
        for reader in openfileobject:
            reader >> buffer

            if buffer == "====":
                break

            print(buffer)

    print("yay")
    dummy = 0

    with open(toyFileName) as openfileobject:
        for reader in openfileobject:
            reader >> dummy
            reader >> dummy # m23, m(pi+ pi-), called m12 in processToyRoot convention.
            reader >> m12  # Already swapped according to D* charge. m12 = m(pi+pi0)
            reader >> m13

            # Errors on Dalitz variables
            reader >> dummy
            reader >> dummy
            reader >> dummy

            reader >> dummy # Decay time
            reader >> dummy # sigma_t

            reader >> dummy # Md0
            reader >> dummy # deltaM
            reader >> dummy # ProbSig
            reader >> dummy # Dst charge
            reader >> dummy # Run
            reader >> dummy # Event
            reader >> dummy # Signal and four bkg fractions.
            reader >> dummy
            reader >> dummy
            reader >> dummy
            reader >> dummy

            eventNumber.setValue(data.getNumEvents())
            data.addEvent()




def makeKzeroVeto():
    print("makeKzeroVeto")
    if kzero_veto:
        return kzero_veto

    kVetoInfo = VetoInfo()
    kVetoInfo.cyclic_index = PAIR_23
    kVetoInfo.minimum      = Variable("veto_min", 0.475 * 0.475)
    kVetoInfo.maximum      = Variable("veto_max", 0.505 * 0.505)
    vetos = (kVetoInfo)
    kzero_veto = DalitzVetoPdf("kzero_veto", m12, m13, motherM, neutrlM, chargeM, chargeM, vetos)
    return kzero_veto


def makeSignalPdf(eff = 0):
    print("makeSignalPdf")
    dtop0pp = DecayInfo()
    dtop0pp.motherMass   = _mD0
    dtop0pp.daug1Mass    = piZeroMass
    dtop0pp.daug2Mass    = piPlusMass
    dtop0pp.daug3Mass    = piPlusMass
    dtop0pp.meson_radius = 1.5

    rhop = ResonancePdf("rhop",
                        RBW,
                        Variable("rhop_amp_real", 1),
                        Variable("rhop_amp_imag", 0),
                        fixedRhoMass,
                        fixedRhoWidth,
                        1,
                        PAIR_12)

    fixAmps = False

    rhom = ResonancePdf("rhom",
                        RBW,
                        Variable("rhom_amp_real", 0.714) if fixAmps else Variable("rhom_amp_real", 0.714, 0.001, 0, 0),
                        Variable("rhom_amp_imag", -0.025) if fixAmps else Variable("rhom_amp_imag", -0.025, 0.1, 0, 0),
                        fixedRhoMass,
                        fixedRhoWidth,
                        1,
                        PAIR_13)

    rho0 = ResonancePdf("rho0",
                        RBW,
                        Variable("rho0_amp_real", 0.565) if fixAmps else Variable("rho0_amp_real", 0.565, 0.001, 0, 0),
                        Variable("rho0_amp_imag", 0.164) if fixAmps else Variable("rho0_amp_imag", 0.164, 0.1, 0, 0),
                        fixedRhoMass,
                        fixedRhoWidth,
                        1,
                        PAIR_23)

    sharedMass = Variable("rhop_1450_mass", 1.465, 0.01, 1.0, 2.0)
    shareWidth = Variable("rhop_1450_width", 0.400, 0.01, 0.01, 5.0)

    rhop_1450 = ResonancePdf("rhop_1450",
                             RBW,
                             Variable("rhop_1450_amp_real", -0.174) if fixAmps else Variable("rhop_1450_amp_real", -0.174, 0.001, 0, 0),
                             Variable("rhop_1450_amp_imag", -0.117) if fixAmps else Variable("rhop_1450_amp_imag", -0.117, 0.1, 0, 0),
                             sharedMass,
                             shareWidth,
                             1,
                             PAIR_12)

    rho0_1450 = ResonancePdf("rho0_1450",
                             RBW,
                             Variable("rho0_1450_amp_real", 0.325) if fixAmps else Variable("rho0_1450_amp_real", 0.325, 0.001, 0, 0),
                             Variable("rho0_1450_amp_imag", 0.057) if fixAmps else Variable("rho0_1450_amp_imag", 0.057, 0.1, 0, 0),
                             sharedMass,
                             shareWidth,
                             1,
                             PAIR_23)

    rhom_1450 = ResonancePdf("rhom_1450",
                             RBW,
                             Variable("rhom_1450_amp_real", 0.788) if fixAmps else Variable("rhom_1450_amp_real", 0.788, 0.001, 0, 0),
                             Variable("rhom_1450_amp_imag", 0.226) if fixAmps else Variable("rhom_1450_amp_imag", 0.226, 0.1, 0, 0),
                             sharedMass,
                             shareWidth,
                             1,
                             PAIR_13)

    sharedMass = Variable("rhop_1700_mass", 1.720, 0.01, 1.6, 1.9)
    shareWidth = Variable("rhop_1700_width", 0.250, 0.01, 0.1, 1.0)

    rhop_1700 = ResonancePdf("rhop_1700",
                             RBW,
                             Variable("rhop_1700_amp_real", 2.151) if fixAmps else Variable("rhop_1700_amp_real", 2.151, 0.001, 0, 0),
                             Variable("rhop_1700_amp_imag", -0.658) if fixAmps else Variable("rhop_1700_amp_imag", -0.658, 0.1, 0, 0),
                             sharedMass,
                             shareWidth,
                             1,
                             PAIR_12)

    rho0_1700 = ResonancePdf("rho0_1700",
                             RBW,
                             Variable("rho0_1700_amp_real", 2.400) if fixAmps else Variable("rho0_1700_amp_real", 2.400, 0.001, 0, 0),
                             Variable("rho0_1700_amp_imag", -0.734) if fixAmps else Variable("rho0_1700_amp_imag", -0.734, 0.1, 0, 0),
                             sharedMass,
                             shareWidth,
                             1,
                             PAIR_23)

    rhom_1700 = ResonancePdf("rhom_1700",
                             RBW,
                             Variable("rhom_1700_amp_real", 1.286) if fixAmps else Variable("rhom_1700_amp_real", 1.286, 0.001, 0, 0),
                             Variable("rhom_1700_amp_imag", -1.532) if fixAmps else Variable("rhom_1700_amp_imag", -1.532, 0.1, 0, 0),
                             sharedMass,
                             shareWidth,
                             1,
                             PAIR_13)

    f0_980    = ResonancePdf("f0_980",
                            RBW,Variable("f0_980_amp_real", 0.008 * (-_mD02))if fixAmps else Variable("f0_980_amp_real", 0.008 * (-_mD02), 0.001, 0, 0),
                            Variable("f0_980_amp_imag", -0.013 * (-_mD02))if fixAmps else Variable("f0_980_amp_imag", -0.013 * (-_mD02), 0.1, 0, 0),
                            Variable("f0_980_mass", 0.980, 0.01, 0.8, 1.2),
                            Variable("f0_980_width", 0.044, 0.001, 0.001, 0.08),
                            0,
                            PAIR_23)

    f0_1370   = ResonancePdf("f0_1370",
                            RBW,
                            Variable("f0_1370_amp_real", -0.058 * (-_mD02))if fixAmps else Variable("f0_1370_amp_real", -0.058 * (-_mD02), 0.001, 0, 0),
                            Variable("f0_1370_amp_imag", 0.026 * (-_mD02))if fixAmps else Variable("f0_1370_amp_imag", 0.026 * (-_mD02), 0.1, 0, 0),
                            Variable("f0_1370_mass", 1.434, 0.01, 1.2, 1.6),
                            Variable("f0_1370_width", 0.173, 0.01, 0.01, 0.4),
                            0,
                            PAIR_23)

    f0_1500   = ResonancePdf("f0_1500",
                            RBW,
                            Variable("f0_1500_amp_real", 0.057 * (-_mD02))if fixAmps else Variable("f0_1500_amp_real", 0.057 * (-_mD02), 0.001, 0, 0),
                            Variable("f0_1500_amp_imag", 0.012 * (-_mD02))if fixAmps else Variable("f0_1500_amp_imag", 0.012 * (-_mD02), 0.1, 0, 0),
                            Variable("f0_1500_mass", 1.507, 0.01, 1.3, 1.7),
                            Variable("f0_1500_width", 0.109, 0.01, 0.01, 0.3),
                            0,
                            PAIR_23)

    f0_1710   = ResonancePdf("f0_1710",
                            RBW,
                            Variable("f0_1710_amp_real", 0.070 * (-_mD02))if fixAmps else Variable("f0_1710_amp_real", 0.070 * (-_mD02), 0.001, 0, 0),
                            Variable("f0_1710_amp_imag", 0.087 * (-_mD02))if fixAmps else Variable("f0_1710_amp_imag", 0.087 * (-_mD02), 0.1, 0, 0),
                            Variable("f0_1710_mass", 1.714, 0.01, 1.5, 2.9),
                            Variable("f0_1710_width", 0.140, 0.01, 0.01, 0.5),
                            0,
                            PAIR_23)

    f2_1270 = ResonancePdf("f2_1270",
                        RBW,
                        Variable("f2_1270_amp_real", -1.027 * (-_mD02inv)) if fixAmps else Variable("f2_1270_amp_real", -1.027 * (-_mD02inv), 0.001, 0, 0),
                        Variable("f2_1270_amp_imag", -0.162 * (-_mD02inv))if fixAmps else Variable("f2_1270_amp_imag", -0.162 * (-_mD02inv), 0.1, 0, 0),
                        Variable("f2_1270_mass", 1.2754, 0.01, 1.0, 1.5),
                        Variable("f2_1270_width", 0.1851, 0.01, 0.01, 0.4),
                        2,
                        PAIR_23)

    f0_600 = ResonancePdf("f0_600",
                        RBW,
                        Variable("f0_600_amp_real", 0.068 * (-_mD02))if fixAmps else Variable("f0_600_amp_real", 0.068 * (-_mD02), 0.001, 0, 0),
                        Variable("f0_600_amp_imag", 0.010 * (-_mD02))if fixAmps else Variable("f0_600_amp_imag", 0.010 * (-_mD02), 0.1, 0, 0),
                        Variable("f0_600_mass", 0.500, 0.01, 0.3, 0.7),
                        Variable("f0_600_width", 0.400, 0.01, 0.2, 0.6),
                        0,
                        PAIR_23)

    nonr = ResonancePdf("nonr",
                        NONRES,
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


    return DalitzPlotPdf("signalPDF", m12, m13, eventNumber, dtop0pp, eff)


def runToyFit(toyFileName):
    print("runToyFit")
    m12 = Variable("m12", 0, 3)
    m13 = Variable("m13", 0, 3)
    m12.numbins = 240
    m13.numbins = 240
    eventNumber = CountingVariable("eventNumber", 0, INT_MAX)
    getToyData(toyFileName)

    signal = makeSignalPdf()
    signal.setData(data)
    signal.setDataSize(data.getNumEvents())
    datapdf = FitManager(signal)
    datapdf.fit()
    return datapdf


def main():
    print("main")
    filename = "dalitz_toyMC_000.txt"
    f= open(filename,"w+")
    return runToyFit(filename)

main()

