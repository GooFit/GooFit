 #script to generate a WS model sample which we will then use to fit with the numerical renormalization method and what is currently done

from __future__ import print_function, division

from goofit import *
import numpy as np
import ROOT
import pandas as pd
import matplotlib.pyplot as plt
import joblib 
import uproot

from hep_ml.reweight import  BinsReweighter, GBReweighter,FoldingReweighter
from hep_ml import reweight

print_goofit_info()

_mD0       = 1.8645
piPlusMass = 0.13957018
KmMass     = .493677

#plotting labels
m12label = r'm$_{\pi^{\mp}\pi^{\pm}}$ [GeV/c$^{2}$]'
m34label = r'm$_{K^{\mp}\pi^{\pm}}$ [GeV/c$^{2}$]'
c12label = r'$cos(\theta_{1})$ [rad]'
c34label = r'$cos(\theta_{2})$ [rad]'
philabel = r'$\phi$ [rad]' 
dTimeLabel = r'$D^{0}$ Decay Time [ns]'

rs_frame = ROOT.RDataFrame("DecayTree",["/hepgpu6-data1/johncob/K3Pi/k3pi_scripts/k3pi_workflow_notebooks/root_files/rs_selected_single_candidates.root"])
rs_selected_frame = rs_frame.Filter('Dst_ReFit_D0_M_best >= 1855 && Dst_ReFit_D0_M_best <= 1875').Filter('deltam_ReFit >= 144.5 && deltam_ReFit <= 146.5')

data_dtime = rs_selected_frame.AsNumpy(columns=["D0_Loki_BPVLTIME"])["D0_Loki_BPVLTIME"]


label_dict = {"m12":m12label,"m34":m34label,"c12":c12label,"c34":c34label,"phi":philabel,"dtime":dTimeLabel}


def plot_generated_events(norm_events,reweighted=False):
    plt.subplots(2,3,figsize=(20,12))
    column_vars = ["m12","m34","c12","c34","phi","dtime"]

    if reweighted:
        for id,column in enumerate(column_vars,1):
            ax = plt.subplot(2,3,id)
            xlim = np.percentile(np.hstack([norm_events[column]]), [0.01, 99.99])
            n1,bins1,patches1 = plt.hist(norm_events[column],weights=(1.0/ norm_events["bdt_weight"])*norm_events["importance_weight"], bins=100,range=xlim,density=True,histtype='step')
            plt.xlabel(label_dict[column])
        plt.savefig("plots/norm_events_distributions_reweighted.png")
    else:
        for id,column in enumerate(column_vars,1):
            ax = plt.subplot(2,3,id)
            xlim = np.percentile(np.hstack([norm_events[column]]), [0.01, 99.99])
            n1,bins1,patches1 = plt.hist(norm_events[column], bins=100,range=xlim,density=True,histtype='step')
            plt.xlabel(label_dict[column])
        plt.savefig("plots/generated_WS.png")

#order of variables array: m12, m34, c12, c34, phi, dtime 
def process_generated_sample(variables,flags):
    m12 = []
    m34 = []
    c12 = []
    c34 = []
    phi = []
    dtime = []
    for i in range(len(flags)):
        if flags[i] == True:
            m12.append(variables[0][i])
            m34.append(variables[1][i])
            c12.append(variables[2][i])
            c34.append(variables[3][i])
            phi.append(variables[4][i])
            dtime.append(variables[5][i])
    return pd.DataFrame({"m12":m12,"m34":m34,"c12":c12,"c34":c34,"phi":phi,"dtime":dtime})

def main():
    DK3P_DI = DecayInfo4t(
        Variable("tau", 0.4101, 0.001, 0.300, 0.500),
        Variable("xmixing", 0.005, 0.001, 0, 0),
        Variable("ymixing", 0.01, 0.001, 0, 0),
        Variable("SqWStoRSrate", 1.0 / np.sqrt(300.0))
        )
    DK3P_DI.meson_radius = 1.5
    DK3P_DI.particle_masses = (_mD0,piPlusMass,piPlusMass,KmMass,piPlusMass)

    
    RhoMass  = Variable("rho_mass", 0.77526)
    RhoWidth = Variable("rho_width", 0.1478)
    KstarM   = Variable("KstarM", 0.89581)
    KstarW   = Variable("KstarW", 0.0474)


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

    LSKRS = (Lineshapes.RBW("rho(770)",RhoMass, RhoWidth, 1, M_12,FF.BL2),
             Lineshapes.RBW("K*(892)bar", KstarM, KstarW, 1, M_34,FF.BL2),
             Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_24,FF.BL2),
             Lineshapes.RBW("K*(892)bar", KstarM, KstarW, 1, M_13,FF.BL2))

    LSKRP = (Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_12,FF.BL2),
             Lineshapes.RBW("K*(892)bar", KstarM, KstarW, 1, M_34,FF.BL2),
             Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_24,FF.BL2),
             Lineshapes.RBW("K*(892)bar", KstarM, KstarW, 1, M_13,FF.BL2))

    LSKRD = (Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_12,FF.BL2),
             Lineshapes.RBW("K*(892)bar", KstarM, KstarW, 1, M_34,FF.BL2),
             Lineshapes.RBW("rho(770)", RhoMass, RhoWidth, 1, M_24,FF.BL2),
             Lineshapes.RBW("K*(892)bar", KstarM, KstarW, 1, M_13,FF.BL2))

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
    DK3P_DI.amplitudes = (Bose_symmetrized_AMP_S,
                            Bose_symmetrized_AMP_P,
                            Bose_symmetrized_AMP_D)

    DK3P_DI.amplitudes_B = ( Bose_symmetrized_AMP_S_B,
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
    dp  = Amp4Body_TD("test", observables, DK3P_DI, res, eff, None, 5000000)

    n_toys = 1
    numEvents = 70000
    particles, variables, weights, flags = dp.GenerateSig(numEvents)
    generated_df = process_generated_sample(variables,flags)
    test_figure = plt.figure()
    plt.hist(generated_df['m12'],bins=100)
    test_figure.savefig("./plots/generated_m12.png")
    print(generated_df['m12'].values)
    #plot_generated_events
    return 0

if __name__ == "__main__":
    assert main() == 0