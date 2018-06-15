#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

from goofit import *
import numpy as np

print_goofit_info()

_mD0       = 1.8645
piPlusMass = 0.13957018
KmMass     = .493677
def main():
    DK3P_DI = DecayInfo4t(
        Variable("tau", 0.4101, 0.001, 0.300, 0.500),
        Variable("xmixing", 0.005, 0.001, 0, 0),
        Variable("ymixing", 0.01, 0.001, 0, 0),
        Variable("SqWStoRSrate", 1.0 / np.sqrt(300.0))
        )
    DK3P_DI.meson_radius = 1.5
    DK3P_DI.particle_masses = (_mD0,piPlusMass,piPlusMass,KmMass,piPlusMass)


    RhoMass  = Variable("rho_mass", 0.77526, 0.01, 0.7, 0.8)
    RhoWidth = Variable("rho_width", 0.1478, 0.01, 0.1, 0.2)
    KstarM   = Variable("KstarM", 0.89581, 0.01, 0.9, 0.1)
    KstarW   = Variable("KstarW", 0.0474, 0.01, 0.1, 0.2)


    SFKRS  = (SpinFactor("SF", SF_4Body .DtoV1V2_V1toP1P2_V2toP3P4_S, _mD0, 0, 1, 2, 3),
              SpinFactor("SF", SF_4Body .DtoV1V2_V1toP1P2_V2toP3P4_S, _mD0, 3, 1, 2, 0))

    SFKRP  = (SpinFactor("SF", SF_4Body .DtoV1V2_V1toP1P2_V2toP3P4_P, _mD0, 0, 1, 2, 3),
              SpinFactor("SF", SF_4Body .DtoV1V2_V1toP1P2_V2toP3P4_P, _mD0, 3, 1, 2, 0))

    SFKRD  = (SpinFactor("SF", SF_4Body .DtoV1V2_V1toP1P2_V2toP3P4_D, _mD0, 0, 1, 2, 3),
              SpinFactor("SF", SF_4Body .DtoV1V2_V1toP1P2_V2toP3P4_D, _mD0, 3, 1, 2, 0))

    SFKF   = (SpinFactor("SF", SF_4Body .DtoVS_VtoP1P2_StoP3P4, _mD0, 2, 3, 0, 1),
              SpinFactor("SF", SF_4Body .DtoVS_VtoP1P2_StoP3P4, _mD0, 2, 0, 3, 1))

    SFKK   = (SpinFactor("SF", SF_4Body .DtoAP1_AtoSP2_StoP3P4, _mD0, 0, 1, 3, 2),
              SpinFactor("SF", SF_4Body .DtoAP1_AtoSP2_StoP3P4, _mD0, 3, 1, 0, 2))

    SFK1R  = (SpinFactor("SF", SF_4Body .DtoAP1_AtoVP2_VtoP3P4, _mD0, 3, 2, 0, 1),
              SpinFactor("SF", SF_4Body .DtoAP1_AtoVP2_VtoP3P4, _mD0, 0, 2, 3, 1))

    SFA1R  = (SpinFactor("SF", SF_4Body .DtoAP1_AtoVP2_VtoP3P4, _mD0, 2, 3, 0, 1),
              SpinFactor("SF", SF_4Body .DtoAP1_AtoVP2_VtoP3P4, _mD0, 2, 0, 3, 1))

    SFA1RD = (SpinFactor("SF", SF_4Body .DtoAP1_AtoVP2Dwave_VtoP3P4, _mD0, 2, 3, 0, 1),
              SpinFactor("SF", SF_4Body .DtoAP1_AtoVP2Dwave_VtoP3P4, _mD0, 2, 0, 3, 1))

    LSKRS = (Lineshapes.RBW("rho(770)",RhoMass, RhoWidth, 1, M_12),
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

    Bose_symmetrized_AMP_S = Amplitude("K*(892)rho(770)_S",
                                        Variable("amp_real1", 1.0),
                                        Variable("amp_imag1", 0.0),
                                        LSKRS,
                                        SFKRS,
                                        2)

    Bose_symmetrized_AMP_P = Amplitude("K*(892)rho(770)_P",
                                        Variable("amp_real2", 0.526),
                                        Variable("amp_imag2", -0.626),
                                        LSKRP,
                                        SFKRP,
                                        2)

    Bose_symmetrized_AMP_D = Amplitude("K*(892)rho(770)_D",
                                        Variable("amp_real3", 26.537),
                                        Variable("amp_imag3", 12.284),
                                        LSKRD,
                                        SFKRD,
                                        2)

    Bose_symmetrized_AMP_S_B = Amplitude("B_K*(892)rho(770)_S",
                                          Variable("amp_real1", 1.0),
                                          Variable("amp_imag1", 0),
                                          LSKRS,
                                          SFKRS,
                                          2)

    Bose_symmetrized_AMP_P_B = Amplitude("B_K*(892)rho(770)_P",
                                          Variable("amp_real2", -0.145),
                                          Variable("amp_imag2", 0.86),
                                          LSKRP,
                                          SFKRP,
                                          2)

    Bose_symmetrized_AMP_D_B = Amplitude("B_K*(892)rho(770)_D",
                                          Variable("amp_real3", 24.343),
                                          Variable("amp_imag3", 5.329),
                                          LSKRD,
                                          SFKRD,
                                          2)

    DK3P_DI.amplitudes_B = (Bose_symmetrized_AMP_S,
                            Bose_symmetrized_AMP_P,
                            Bose_symmetrized_AMP_D,
                            Bose_symmetrized_AMP_S_B,
                            Bose_symmetrized_AMP_P_B,
                            Bose_symmetrized_AMP_D_B)



    m12                 = Observable("m12", 0, 3)
    m34                 = Observable("m34", 0, 3)
    cos12               = Observable("cos12", -1, 1)
    cos34               = Observable("m12", -1, 1)
    phi                 = Observable("phi", -3.5, 3.5)
    eventNumber         = EventNumber("eventNumber")
    dtime               = Observable("dtime", 0, 10)
    sigmat              = Observable("sigmat", -3, 3)
    constantOne         = Variable("constantOne", 1)
    constantZero        = Variable("constantZero", 0)

    observables = (m12,m34,cos12,cos34,phi,eventNumber,dtime,sigmat)
    offsets = (constantZero, constantZero,)
    coefficients = (constantOne,)


    res = TruthResolution()
    eff = PolynomialPdf("constantEff", observables, coefficients, offsets, 0)
    dp  = Amp4Body_TD("test", observables, DK3P_DI, res, eff, None, 1)


    for k in range(10):
        numEvents = 800000

        dp.setGenerationOffset(k * numEvents)
        particles, variables, weights, flags = dp.GenerateSig(numEvents)

    return 0

if __name__ == "__main__":
    assert main() == 0
