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



def fit_array(arr,varname):
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
    c.SaveAs(varname+"_fit_xy_floated.png")



x_values, x_err, y_values, y_err = get_xy_values()
fit_array(x_values,"fitted x")
#fit_array(x_err,"fitted x error")
fit_array(y_values,"fitted y")
#fit_array(y_err,"fitted y error")

"""
x = ROOT.RooRealVar("x","x",0.0044,0.005)
mean = ROOT.RooRealVar("mean","x Gaussian mean",0.0045, 0.004,0.005)
sigma = ROOT.RooRealVar("sigma","x Gaussian width",0.001,0.,1.)
gauss = ROOT.RooGaussian("gauss","gaussian PDF",x,mean,sigma)
x_dataset = ROOT.RooDataSet("x_data","x_data",ROOT.RooArgSet(x))
for x_val in x_values:
    x.setVal(x_val)
    x_dataset.add(ROOT.RooArgSet(x))
gauss.fitTo(x_dataset)

c = ROOT.TCanvas()
xframe = x.frame()
x_dataset.plotOn(xframe)
gauss.plotOn(xframe)
xframe.GetXaxis().SetTitle("Fitted x value")
xframe.SetTitle("")
xframe.Draw("E0")
xline = ROOT.TLine(0.0045,0,0.0045,15)
xline.SetLineColor(ROOT.kRed)
xline.SetLineWidth(5)
xline.Draw()
c.SaveAs("fitted_x_values_xy_floated.png")
"""