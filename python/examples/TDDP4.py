#!/usr/bin/env python
# -*- coding: utf-8 -*-

from goofit import *

import math

_mD0       = 1.8645
piPlusMass = 0.13957018
KmMass     = .493677
'''
TString output = "test_10_15.output"
app.add_option("-o,--output,output", output, "File to output", true).check(GooFit .NonexistentPath)

trials = 100
app.add_option("-t,--trials,output", trials, "Number of trials", true)
'''
DK3P_DI = DecayInfo_DP()
DK3P_DI.meson_radius = 1.5
DK3P_DI.particle_masses = (_mD0,piPlusMass,piPlusMass,KmMass,piPlusMass)


RhoMass  = Variable("rho_mass", 0.77526, 0.01, 0.7, 0.8)
RhoWidth = Variable("rho_width", 0.1478, 0.01, 0.1, 0.2)
KstarM   = Variable("KstarM", 0.89581, 0.01, 0.9, 0.1)
KstarW   = Variable("KstarW", 0.0474, 0.01, 0.1, 0.2)


SFKRS  = (SpinFactor("SF", SF_4Body .DtoV1V2_V1toP1P2_V2toP3P4_S, 0, 1, 2, 3),
          SpinFactor("SF", SF_4Body .DtoV1V2_V1toP1P2_V2toP3P4_S, 3, 1, 2, 0))

SFKRP  = (SpinFactor("SF", SF_4Body .DtoV1V2_V1toP1P2_V2toP3P4_P, 0, 1, 2, 3),
          SpinFactor("SF", SF_4Body .DtoV1V2_V1toP1P2_V2toP3P4_P, 3, 1, 2, 0))

SFKRD  = (SpinFactor("SF", SF_4Body .DtoV1V2_V1toP1P2_V2toP3P4_D, 0, 1, 2, 3),
          SpinFactor("SF", SF_4Body .DtoV1V2_V1toP1P2_V2toP3P4_D, 3, 1, 2, 0))

SFKF   = (SpinFactor("SF", SF_4Body .DtoVS_VtoP1P2_StoP3P4, 2, 3, 0, 1),
          SpinFactor("SF", SF_4Body .DtoVS_VtoP1P2_StoP3P4, 2, 0, 3, 1))

SFKK   = (SpinFactor("SF", SF_4Body .DtoAP1_AtoSP2_StoP3P4, 0, 1, 3, 2),
          SpinFactor("SF", SF_4Body .DtoAP1_AtoSP2_StoP3P4, 3, 1, 0, 2))

SFK1R  = (SpinFactor("SF", SF_4Body .DtoAP1_AtoVP2_VtoP3P4, 3, 2, 0, 1),
          SpinFactor("SF", SF_4Body .DtoAP1_AtoVP2_VtoP3P4, 0, 2, 3, 1))

SFA1R  = (SpinFactor("SF", SF_4Body .DtoAP1_AtoVP2_VtoP3P4, 2, 3, 0, 1),
          SpinFactor("SF", SF_4Body .DtoAP1_AtoVP2_VtoP3P4, 2, 0, 3, 1))

SFA1RD = (SpinFactor("SF", SF_4Body .DtoAP1_AtoVP2Dwave_VtoP3P4, 2, 3, 0, 1),
          SpinFactor("SF", SF_4Body .DtoAP1_AtoVP2Dwave_VtoP3P4, 2, 0, 3, 1))

LSKRS = (Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_12, LS .BW, FF .BL2),
         Lineshape("K*(892)bar", KstarM, KstarW, 1, M_34, LS .BW, FF .BL2),
         Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_24, LS .BW, FF .BL2),Lineshape("K*(892)bar", KstarM, KstarW, 1, M_13, LS .BW, FF .BL2))
LSKRP = (Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_12, LS .BW, FF .BL2),Lineshape("K*(892)bar", KstarM, KstarW, 1, M_34, LS .BW, FF .BL2),Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_24, LS .BW, FF .BL2),Lineshape("K*(892)bar", KstarM, KstarW, 1, M_13, LS .BW, FF .BL2))
LSKRD = (Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_12, LS .BW, FF .BL2),Lineshape("K*(892)bar", KstarM, KstarW, 1, M_34, LS .BW, FF .BL2),Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_24, LS .BW, FF .BL2),Lineshape("K*(892)bar", KstarM, KstarW, 1, M_13, LS .BW, FF .BL2))

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


DK3P_DI._tau          = Variable("tau", 0.4101, 0.001, 0.300, 0.500)
DK3P_DI._xmixing      = Variable("xmixing", 0.005, 0.001, 0, 0)
DK3P_DI._ymixing      = Variable("ymixing", 0.01, 0.001, 0, 0)
DK3P_DI._SqWStoRSrate = Variable("SqWStoRSrate", 1.0 / math.sqrt(300.0))

m12                 = Variable("m12", 0, 3)
m34                 = Variable("m34", 0, 3)
cos12               = Variable("cos12", -1, 1)
cos34               = Variable("m12", -1, 1)
phi                 = Variable("phi", -3.5, 3.5)
eventNumber         = CountingVariable("eventNumber", 0, INT_MAX)
dtime               = Variable("dtime", 0, 10)
sigmat              = Variable("sigmat", -3, 3)
constantOne         = Variable("constantOne", 1)
constantZero        = Variable("constantZero", 0)

observables = (m12,m34,cos12,cos34,phi,eventNumber,dtime,sigmat)
offsets = (constantZero, constantZero,)
coefficients = (constantOne,)


dat = TruthResolution()
eff = PolynomialPdf("constantEff", observables, coefficients, offsets, 0)
dp  = TDDP4("test", observables, DK3P_DI, dat, eff, None, 1)

'''
file = TFile(output, "RECREATE")
tree = TTree("events", "events")



tree.Branch("m12", tm12, "m12/D")
tree.Branch("m34", tm34, "m34/D")
tree.Branch("c12", tc12, "c12/D")
tree.Branch("c34", tc34, "c34/D")
tree.Branch("phi", tphi, "phi/D")
tree.Branch("dtime", tdtime, "dtime/D")
tree.Branch("D0_E", D0_E, "D0_E/D")
tree.Branch("D0_Px", D0_Px, "D0_Px/D")
tree.Branch("D0_Py", D0_Py, "D0_Py/D")
tree.Branch("D0_Pz", D0_Pz, "D0_Pz/D")
tree.Branch("D0_pdg", D0_pdg, "D0_pdg/I")
tree.Branch("Kplus_E", Kplus_E, "Kplus_E/D")
tree.Branch("Kplus_Px", Kplus_Px, "Kplus_Px/D")
tree.Branch("Kplus_Py", Kplus_Py, "Kplus_Py/D")
tree.Branch("Kplus_Pz", Kplus_Pz, "Kplus_Pz/D")
tree.Branch("Kplus_pdg", Kplus_pdg, "Kplus_pdg/I")
tree.Branch("Piminus1_E", Piminus1_E, "Piminus1_E/D")
tree.Branch("Piminus1_Px", Piminus1_Px, "Piminus1_Px/D")
tree.Branch("Piminus1_Py", Piminus1_Py, "Piminus1_Py/D")
tree.Branch("Piminus1_Pz", Piminus1_Pz, "Piminus1_Pz/D")
tree.Branch("Piminus1_pdg", Piminus1_pdg, "Piminus1_pdg/I")
tree.Branch("Piminus2_E", Piminus2_E, "Piminus2_E/D")
tree.Branch("Piminus2_Px", Piminus2_Px, "Piminus2_Px/D")
tree.Branch("Piminus2_Py", Piminus2_Py, "Piminus2_Py/D")
tree.Branch("Piminus2_Pz", Piminus2_Pz, "Piminus2_Pz/D")
tree.Branch("Piminus2_pdg", Piminus2_pdg, "Piminus2_pdg/I")
tree.Branch("Piplus_E", Piplus_E, "Piplus_E/D")
tree.Branch("Piplus_Px", Piplus_Px, "Piplus_Px/D")
tree.Branch("Piplus_Py", Piplus_Py, "Piplus_Py/D")
tree.Branch("Piplus_Pz", Piplus_Pz, "Piplus_Pz/D")
tree.Branch("Piplus_pdg", Piplus_pdg, "Piplus_pdg/I")
'''
for k in range(5):
    numEvents = 800000
    dp.setGenerationOffset(k * numEvents)
    tuple = dp.GenerateSig(numEvents)

    particles = tuple[0]
    variables = tuple[1]
    weights   = tuple[2]
    flags     = tuple[3]
    '''
    accepted   = thrust .count_if(flags.begin(), flags.end(), thrust .identity<bool>())
    GOOFIT_INFO(
        "Run #{}: Using accept-reject method would leave you with {} out of {} events", k, accepted, numEvents)
        '''
    k= k+1


    while i<weights.size():
        if flags[i] == 1:
            tm12         = (variables[0])[i]
            tm34         = (variables[1])[i]
            tc12         = (variables[2])[i]
            tc34         = (variables[3])[i]
            tphi         = (variables[4])[i]
            tdtime       = (variables[5])[i]
            D0_E         = 1864
            D0_Px        = 0.0
            D0_Py        = 0.0
            D0_Pz        = 0.0
            D0_pdg       = 421
            Kplus_E      = 1000 * (particles[2])[i].get(0)
            Kplus_Px     = 1000 * (particles[2])[i].get(1)
            Kplus_Py     = 1000 * (particles[2])[i].get(2)
            Kplus_Pz     = 1000 * (particles[2])[i].get(3)
            Kplus_pdg    = -321
            Piminus1_E   = 1000 * (particles[3])[i].get(0)
            Piminus1_Px  = 1000 * (particles[3])[i].get(1)
            Piminus1_Py  = 1000 * (particles[3])[i].get(2)
            Piminus1_Pz  = 1000 * (particles[3])[i].get(3)
            Piminus1_pdg = 211
            Piminus2_E   = 1000 * (particles[0])[i].get(0)
            Piminus2_Px  = 1000 * (particles[0])[i].get(1)
            Piminus2_Py  = 1000 * (particles[0])[i].get(2)
            Piminus2_Pz  = 1000 * (particles[0])[i].get(3)
            Piminus2_pdg = 211
            Piplus_E     = 1000 * (particles[1])[i].get(0)
            Piplus_Px    = 1000 * (particles[1])[i].get(1)
            Piplus_Py    = 1000 * (particles[1])[i].get(2)
            Piplus_Pz    = 1000 * (particles[1])[i].get(3)
            Piplus_pdg   = -211

            tree.Fill()




'''
    delete variables[0]
    delete variables[1]
    delete variables[2]
    delete variables[3]
    delete variables[4]
    delete variables[5]

    delete particles[0]
    delete particles[1]
    delete particles[2]
    delete particles[3]
'''

tree.Write()
file.Close()
#return 0

