#fit a toy with realistic acceptance and decay time weights incorporated (weights set to 1)
#!/usr/bin/env python
# -*- coding: utf-8 -*-

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


columns = ['m12','m34','c12','c34','phi','dtime']

#plotting labels
m12label = r'm$_{\pi^{\mp}\pi^{\pm}}$ [GeV/c$^{2}$]'
m34label = r'm$_{K^{\mp}\pi^{\pm}}$ [GeV/c$^{2}$]'
c12label = r'$cos(\theta_{1})$ [rad]'
c34label = r'$cos(\theta_{2})$ [rad]'
philabel = r'$\phi$ [rad]' 
dTimeLabel = r'$D^{0}$ Decay Time [ns]'

#read in generated MC data that we will use to test the fit
mc_df = pd.read_csv('amp4body_cpp_generated_events.txt',names=columns,delimiter=" ")
mc_df['eff_weight'] = np.ones(len(mc_df))
print(mc_df.head())

label_dict = {"m12":m12label,"m34":m34label,"c12":c12label,"c34":c34label,"phi":philabel,"dtime":dTimeLabel}
imp_samp_time=1.
#this function calculates the corresponding importance sampled weight using the triangular function generated decay times
def two_parts(xs, turn=0.6):
    ys = np.ones(len(xs))
    ys[xs<turn] = xs[xs<turn]*np.exp(-turn/imp_samp_time)/turn
    ys[xs>=turn] = np.exp(-xs[xs>=turn]/imp_samp_time)
    norm = 0.5*np.exp(-turn/imp_samp_time)/turn*turn**2 - 0.5*np.exp(-turn/imp_samp_time)/turn*0.1725**2
    norm += (-imp_samp_time*np.exp(-3.26/imp_samp_time)+imp_samp_time*np.exp(-turn/imp_samp_time))
    return ys/norm

#this generates the decay time distribution for the normalisation events
def two_parts_generate(turn=0.6, size=1):
    int_left = 0.5*np.exp(-turn/imp_samp_time)/turn*turn**2 - 0.5*np.exp(-turn/imp_samp_time)/turn*0.1725**2
    int_right = (-imp_samp_time*np.exp(-3.26/imp_samp_time)+imp_samp_time*np.exp(-turn/imp_samp_time))
    print(int_left, int_right)
    fl = int_left/(int_left+int_right)
    print(fl)
    left = np.random.triangular(0, turn, turn, size=int(fl*size*3.0))
    left = left[left>0.1725][:int(fl*size)]
    right = np.random.exponential(imp_samp_time, size=int((size-int(fl*size))*1.4)) + turn
    return np.append(left, right[right<3.26][:size-int(fl*size)])

def get_bdt_weights(df):
    bdt_model = joblib.load("/hepgpu6-data1/johncob/K3Pi/k3pi_scripts/k3pi_workflow_notebooks/bdt_models/rs_efficiency_model_kfolded_6d_ord_full_stats.sav")
    training_vars_ord = ["m12_ord","m34_ord","cos12_ord","cos34_ord","phi_ord","log_dtime"]
    bdt_df = df.rename(columns={"m12":"m12_ord","m34":"m34_ord","c12":"cos12_ord","c34":"cos34_ord","phi":"phi_ord"})
    bdt_weights = bdt_model.predict_weights(bdt_df[training_vars_ord])
    return bdt_weights


def main():
    DK3P_DI = DecayInfo4t(
        Variable("tau", 0.4101, 0.001, 0.300, 0.500),
        Variable("xmixing", 0.005, 0.001, 0., 0.),
        Variable("ymixing", 0.01, 0.001, 0, 0.),
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
                                          Variable("amp_real1_B", 1.0),
                                          Variable("amp_imag1_B", 0),
                                          LSKRS,
                                          SFKRS,
                                          2)

    Bose_symmetrized_AMP_P_B = Amplitude("B_K*(892)rho(770)_P",
                                          Variable("amp_real2_B", -0.145),
                                          Variable("amp_imag2_B", 0.86),
                                          LSKRP,
                                          SFKRP,
                                          2)

    Bose_symmetrized_AMP_D_B = Amplitude("B_K*(892)rho(770)_D",
                                          Variable("amp_real3_B", 24.343),
                                          Variable("amp_imag3_B", 5.329),
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
    cos34               = Observable("cos34", -1, 1)
    phi                 = Observable("phi", -3.5, 3.5)
    eventNumber         = EventNumber("eventNumber")
    dtime               = Observable("dtime", 0, 10)
    sigmat              = Observable("sigmat", -3, 3)
    data_eff_weight     = Observable("data_eff_weight",-1000,1000)
    constantOne         = Variable("constantOne", 1)
    constantZero        = Variable("constantZero", 0)

    observables = (m12,m34,cos12,cos34,phi,eventNumber,dtime,sigmat,data_eff_weight)
    offsets = (constantZero, constantZero,)
    coefficients = (constantOne,)

    #create the dataset to fit
    mc_df['eventNumber'] = [i for i in range(len(mc_df))]
    dataset = UnbinnedDataSet(m12,m34,cos12,cos34,phi,dtime,data_eff_weight,eventNumber)
    dataset.from_numpy(mc_df.to_numpy().T)

    for i in range(10):
        print(dataset.getValue(m12,i),dataset.getValue(m34,i),dataset.getValue(cos12,i),dataset.getValue(cos34,i),dataset.getValue(phi,i),dataset.getValue(dtime,i),dataset.getValue(data_eff_weight,i),dataset.getValue(eventNumber,i))

    res = TruthResolution()
    eff = PolynomialPdf("constantEff", observables, coefficients, offsets, 0)
    dp  = Amp4Body_TD("test", observables, DK3P_DI, res, eff, None, 5000000)
    print(f"Configuring dataset of {dataset.getNumEvents()} events")
    #configure amp4body parameters for the fit
    dp.setDataSize(dataset.getNumEvents(), 9)
    fit_manager = FitManager(dp)
    print("Fitting...")
    fitResults = fit_manager.fit()
    convergedStatus = fitResults.IsValid()
    print(convergedStatus)

    #do the integration normally without modifying the renormalization events

    #plot_decay_time_distribution(norm_events_df)
    #plot_normalisation_events(norm_events_df,reweighted=True)

    """
    #required for realistic fit model
    dp.set_special_integral(False)
    m12_arr = dp.get_norm_m12()
    m34_arr = dp.get_norm_m34()
    c12_arr = dp.get_norm_c12()
    c34_arr = dp.get_norm_c34()
    phi_arr = dp.get_norm_phi()
    dtime_arr = dp.get_norm_dtime()
    importance_weights = dp.get_norm_importance_weights() #importance weights for decay time importance sampling
    eff_weights = dp.get_norm_eff()
    columns = {"m12":m12_arr,"m34":m34_arr,"c12":c12_arr,"c34":c34_arr,"phi":phi_arr,"dtime":dtime_arr,"importance_weight":importance_weights,"eff_weights":eff_weights}
    norm_events_df = pd.DataFrame(columns)

    print("generating weights and decay times")
    norm_events_df['dtime'] = two_parts_generate(size=len(norm_events_df))
    norm_events_df['importance_weight'] = two_parts(norm_events_df['dtime'])
    norm_events_df['log_dtime'] = np.log(norm_events_df['dtime'])
    norm_events_df["bdt_weight"] = get_bdt_weights(norm_events_df)
    #copy back generated data
    print("copying back modified arrays")
    dp.set_norm_dtime(norm_events_df.dtime.values)
    dp.set_norm_eff(norm_events_df.bdt_weight.values)
    dp.set_norm_importance_weights(norm_events_df.importance_weight.values)
    """

    return 0

if __name__ == "__main__":
    assert main() == 0
