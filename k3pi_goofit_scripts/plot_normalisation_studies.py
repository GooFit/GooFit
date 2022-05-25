#plotting script to produce plots of fits for individual fit parameters and their pulls (e.g. Gaussian fits or plots like those of 6.14/15 in Dominik's thesis)
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import subprocess as sb
import os
import glob
import ROOT

files = glob.glob("../build-cuda/k3pi_goofit_scripts/ws_signal_only_fit_results/*.txt")
true_x = 0.0045
true_y = 0.0062
def get_xy_values():
    files = glob.glob("../build-cuda/k3pi_goofit_scripts/ws_signal_only_fit_results/signal_only_70k_sample_x_0.0045_y_0.0062_fit_results_seed_*.txt")
    x_values = []
    x_err = []
    y_values = []
    y_err = []
    fit_files = []
    for f in files:
        if "floated_DCS" not in f:
            fit_files.append(f)
            file_values = np.loadtxt(f,delimiter=' ')
            if file_values[5] == 1 and file_values[6] == 3:
                x_values.append(file_values[1])
                x_err.append(file_values[2])
                y_values.append(file_values[3])
                y_err.append(file_values[4])
    return (x_values,x_err,y_values,y_err)

def get_xy_dcs_values():
    files = glob.glob("../build-cuda/k3pi_goofit_scripts/ws_signal_only_fit_results/signal_only_70k_sample_x_0.0045_y_0.0062_fit_results_seed_*.txt")
    x_values = []
    x_err = []
    y_values = []
    y_err = []
    fit_files = []
    for f in files:
        if "floated_DCS" in f:
            fit_files.append(f)
            file_values = np.loadtxt(f,delimiter=' ',dtype=str)
            if int(file_values[5]) == 1 and int(file_values[6]) == 3:
                x_values.append(float(file_values[1]))
                x_err.append(float(file_values[2]))
                y_values.append(float(file_values[3]))
                y_err.append(float(file_values[4]))
    return (x_values,x_err,y_values,y_err)

def get_xy_ensemble_values():
    files = glob.glob("../build-cuda/k3pi_goofit_scripts/ws_signal_only_fit_results/signal_only_70k_sample_x_0.0045_y_0.0062_fit_results_sample_num_*_same_norm_seed.txt")
    files_default = glob.glob("../build-cuda/k3pi_goofit_scripts/ws_signal_only_fit_results/signal_only_70k_sample_x_0.0045_y_0.0062_fit_results_sample_num_*_same_norm_seed_default_integral.txt")
    x_values = []
    y_values = []
    x_values_default = []
    y_values_default = []
    fit_files = []
    for f,f_default in zip(files,files_default):
        if "floated_DCS" not in f:
            fit_files.append(f)
            file_values = np.loadtxt(f,delimiter=' ',dtype=str)
            if int(file_values[5]) == 1 and int(file_values[6]) == 3:
                x_values.append(float(file_values[1]))
                y_values.append(float(file_values[3]))
            fit_values = np.loadtxt(f_default,delimiter=' ',dtype=str)
            if int(file_values[5]) == 1 and int(file_values[6]) == 3:
                x_values_default.append(float(file_values[1]))
                y_values_default.append(float(file_values[3]))
    return (x_values,x_values_default,y_values,y_values_default)

def get_xy_ensemble_dcs_values():
    files = glob.glob("../build-cuda/k3pi_goofit_scripts/ws_signal_only_fit_results/signal_only_70k_sample_x_0.0045_y_0.0062_fit_results_sample_num_*_same_norm_seed_floated_dcs.txt")
    files_default = glob.glob("../build-cuda/k3pi_goofit_scripts/ws_signal_only_fit_results/signal_only_70k_sample_x_0.0045_y_0.0062_fit_results_sample_num_*_same_norm_seed_floated_dcs_default_integral.txt")
    x_values = []
    y_values = []
    x_values_default = []
    y_values_default = []
    fit_files = []
    for f,f_default in zip(files,files_default):
        with open(f) as f_open:
            file_values = [x.split() for x in f_open.readlines()]
            file_values = [item for sublist in file_values for item in sublist]
            if int(file_values[5]) == 1 and int(file_values[6]) == 3:
                x_values.append(float(file_values[1]))
                y_values.append(float(file_values[3]))
            f_open.close()
        with open(f_default) as f_open:
            file_values = [x.split() for x in f_open.readlines()]
            file_values = [item for sublist in file_values for item in sublist]
            if int(file_values[5]) == 1 and int(file_values[6]) == 3:
                x_values_default.append(float(file_values[1]))
                y_values_default.append(float(file_values[3]))
            f_open.close()
    return (x_values,x_values_default,y_values,y_values_default)


def fit_array(arr,varname,dcs=False,ensemble=False,default_int=False):
    low_lim,high_lim =( 0.95*min(arr),1.05*max(arr))
    var = ROOT.RooRealVar(varname,varname,0.9*min(arr),1.1*max(arr))
    mean = ROOT.RooRealVar(varname + " mean",varname + " Gaussian mean",(low_lim + high_lim)/2, low_lim,high_lim)
    sigma = ROOT.RooRealVar(varname + " sigma",varname + " Gaussian width",0.001,0.,2.)
    gauss = ROOT.RooGaussian(varname+"_gauss",varname + " gaussian PDF",var,mean,sigma)
    dataset = ROOT.RooDataSet(varname +"_dataset",varname + "_data",ROOT.RooArgSet(var))
    for val in arr:
        var.setVal(val)
        dataset.add(ROOT.RooArgSet(var))
    gauss.fitTo(dataset)
    c = ROOT.TCanvas()
    varframe = var.frame()
    dataset.plotOn(varframe)
    gauss.plotOn(varframe)
    varframe.GetXaxis().SetTitle(varname)
    varframe.SetTitle("")
    pt = ROOT.TPaveText(0.05,.1,0.6,0.8)
    pt.AddText(f"#mu = {mean.getValV()}")
    pt.AddText(f"#sigma = {sigma.getValV()}")

    varframe.Draw()
    if varname == "fitted x":
        xline = ROOT.TLine(0.0045,0,0.0045,15)
        xline.SetLineColor(ROOT.kRed)
        xline.SetLineWidth(5)
        xline.Draw()
    elif varname == "fitted y":
        xline = ROOT.TLine(0.0062,0,0.0062,15)
        xline.SetLineColor(ROOT.kRed)
        xline.SetLineWidth(5)
        xline.Draw()
    pt.Draw()
    canvas_name = varname +"_fit_xy"
    if dcs:
        canvas_name += "_dcs_floated"
    if ensemble:
        canvas_name += "_ensemble"
    if default_int:
        canvas_name += "_default_integral"
    canvas_name += ".png"
    c.SaveAs(canvas_name)
    print(f"Fit results for {varname} with dcs floated:{dcs}")
    print(f"Mean:{mean.getValV()}")
    print(f"Sigma:{sigma.getValV()}")


#results for fitting unique toys with the same set of normalisation events
x_values,x_values_default,y_values,y_values_default = get_xy_ensemble_values()

x_values_dcs,x_values_default_dcs,y_values_dcs,y_values_default_dcs =  get_xy_ensemble_dcs_values()
fit_array(x_values_dcs,"fitted x",dcs=True,ensemble=True)
fit_array(y_values_dcs,"fitted y",dcs=True,ensemble=True)
fit_array(x_values_default_dcs,"fitted x",dcs=True,ensemble=True,default_int=True)
fit_array(y_values_default_dcs,"fitted y",dcs=True,ensemble=True,default_int=True)

fit_array(x_values,"fitted x",dcs=False,ensemble=True)
fit_array(y_values,"fitted y",dcs=False,ensemble=True)
fit_array(x_values_default,"fitted x",dcs=False,ensemble=True,default_int=True)
fit_array(y_values_default,"fitted y",dcs=False,ensemble=True,default_int=True)

"""
#results for fitting the same toy multiple times. Used for estimating error on 6D integration
x_values, x_err, y_values, y_err = get_xy_values()
x_values_dcs,x_err_dcs,y_values_dcs,y_err_dcs = get_xy_dcs_values()
print(x_values_dcs)
fit_array(x_values_dcs,"fitted x",dcs=True)
fit_array(y_values_dcs,"fitted y",dcs=True)
fit_array(x_values,"fitted x",dcs=False)
fit_array(y_values,"fitted y",dcs=False)
"""