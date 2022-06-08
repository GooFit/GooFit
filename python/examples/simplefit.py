# -*- coding: utf-8 -*-
import argparse
import sys
from math import atan2, cos, pi, sin, sqrt

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colors
from root_pandas import read_root
from scipy.optimize import curve_fit
from scipy.stats import norm

from goofit import *

# from fitter_tools import *


print_goofit_info()

m12 = Observable("m12", 0, 3.5)
m13 = Observable("m13", 0, 3.5)


def getData(filename, tree, catIndex):
    df = read_root(
        filename,
        tree,
        columns=["m12", "m13", "sig", "__index__", "B_DTFKS_D0_CTAU", "B_DTFKS_D0_M"],
    )
    df = df.query("B_DTFKS_D0_CTAU/0.3 < 0.4101")
    df["category"] = catIndex
    # df = df.query('B_DTFKS_D0_M > 1834 & B_DTFKS_D0_M < 1894')
    df.reset_index(drop=True, inplace=True)
    df["index"] = range(df.shape[0])

    print(df)
    array = df[["m12", "m13", "sig", "category", "index"]].to_numpy()
    dataset = UnbinnedDataSet(m12, m13, sigprob, category, eventNumber)
    dataset.from_matrix(array.T)

    return dataset


def getDatasets(years, modes):
    datasets = {}
    filename = "root://eoslhcb.cern.ch//eos/lhcb/user/m/mhilton/KSPiPi-ntuples/tuples-BDT/sWeight_{mode}_{year}.root"
    catIndex = 0
    for year in years:
        for mode in modes:
            print(f"using category {catIndex} for {mode} {year}")
            dataset = getData(
                filename.format(mode=mode, year=year), "tree_sWeights", catIndex
            )
            datasets[f"{year}_{mode}"] = dataset
            catIndex = catIndex + 1
    for key, value in datasets.items():
        print(f"Number of events in {key} dataset: ", value.getNumEvents())
    return datasets


def make_toy_data(dp, m12, m13, sigprob, eventNumber):
    print("make toy data")
    data = UnbinnedDataSet(m12, m13, sigprob, eventNumber)

    print(m12.getNumBins())
    print(m13.getNumBins())

    xbins, ybins = [], []
    for i in range(m12.getNumBins()):
        m12.setValue(m12.getLowerLimit() + m12.getBinSize() * (i + 0.5))
        for j in range(m13.getNumBins()):
            m13.setValue(m13.getLowerLimit() + m13.getBinSize() * (j + 0.5))
            if inDalitz(m12.getValue(), m13.getValue()):
                xbins.append(i)
                ybins.append(j)
                sigprob.setValue(0.9)
                data.addEvent()
                eventNumber.setValue(eventNumber.getValue() + 1)
    dp.setData(data)

    print("Normalise")
    dp.normalize()
    pdfvals = dp.getCompProbsAtDataPoints()  # Memory error
    print("pdfvals", pdfvals)

    return data


def makeFit(years, modes, plotprefix="", random_seed=-1):

    if random_seed > -1:
        plotprefix = f"{plotprefix}rand{random_seed}_"
    constantOne = Variable("One", 1)
    constantZero = Variable("Zero", 0)
    dalitz_smoothing = Variable("dalitz_smoothing", 0.5)
    # eff = PolynomialPdf('constantEff', [m12, m13], [constantOne], [constantZero,constantZero], 0)

    m12.setNumBins(4 * 175)
    m13.setNumBins(4 * 175)

    # modes = ['SingleTag_D0ToKsPiPiLL', 'SingleTag_D0ToKsPiPiDD']
    datasets = getDatasets(years, modes)

    arrays = [a.to_matrix() for a in datasets.values()]
    tot_matrix = np.concatenate(arrays, axis=1)
    # tot_matrix[3] = np.arange(tot_matrix.shape[1])
    print(tot_matrix)

    data_combined = UnbinnedDataSet(m12, m13, sigprob, category, eventNumber)
    data_combined.from_matrix(tot_matrix)  # check this works
    print(
        "Number of events in combined dataset: {}".format(data_combined.getNumEvents())
    )
    # try out resetin cachecounter
    Amp3Body.resetCacheCounter()

    globaldecinfo = makeDecayInfo()
    signal_pdfs = {}
    bkg_pdfs = {}
    total_pdfs = {}
    filename = "root://eoslhcb.cern.ch//eos/lhcb/user/m/mhilton/KSPiPi-ntuples/tuples-BDT/sWeight_{mode}_{year}.root"
    for year in years:
        for mode in modes:
            efficiency = makeEfficiency(mode, year)
            yearmode = f"{year}_{mode}"
            signal_pdf = Amp3Body(
                f"signalPDF_{yearmode}",
                m12,
                m13,
                eventNumber,
                globaldecinfo,
                efficiency,
            )
            signal_pdfs[yearmode] = signal_pdf

            bghist_dalitz = BinnedDataSet(m12, m13)
            fillBGhist(filename.format(mode=mode, year=year), bghist_dalitz)
            BGPdf = SmoothHistogramPdf(
                f"BGPdf_{yearmode}", bghist_dalitz, dalitz_smoothing
            )
            bkg_pdfs[yearmode] = BGPdf

            total_pdfs[yearmode] = EventWeightedAddPdf(
                f"total_{yearmode}", [sigprob], [signal_pdf, BGPdf]
            )
    print(total_pdfs)
    stepFunction = BinTransformPdf(
        "stepFunction", [category], [-0.5], [1], [len(datasets)]
    )
    finalPDF_comb = MappedPdf("finalPDF_comb", stepFunction, list(total_pdfs.values()))

    for param in finalPDF_comb.getParameters():
        if "mass" in param.getName() or "width" in param.getName():
            if fixMass:
                param.setFixed(True)
        if "beta" in param.getName() or "f_prod" in param.getName():
            if fixKmatrix:
                param.setFixed(True)  # Fix K-matrix parameters
        if "amp_real" in param.getName() or "amp_imag" in param.getName():
            param.setFixed(False)
        if "prod1" in param.getName() or "prod5" in param.getName():
            param.setFixed(True)
    params_df = pd.read_csv("Fit_params_SingleTag_D0ToKsPiPiLL.csv")
    for column in params_df.columns:
        for param in finalPDF_comb.getParameters():
            if column == param.getName():
                print(
                    f"setting parameter {param.getName()} from {param.getValue()} to {params_df[column][0]}"
                )
                param.setValue(params_df[column][0])
                print(f"its values is now {param.getValue()}")

    if random_seed > -1:
        print("resampling starting values of free parameters")
        for param in finalPDF_comb.getParameters():
            if not param.IsFixed():
                param_val = param.getValue()
                param_err = param.getError()
                np.random.seed((random_seed * 537) + 142)
                resampled_param = np.random.normal(param_val, param_err)
                print(
                    f"resampling parameter {param.getName()} from {param.getValue()} to {resampled_param}"
                )
                param.setValue(resampled_param)

    rho_770_amp_imag.setFixed(True)
    rho_770_amp_real.setFixed(True)

    m12.setNumBins(8 * 175)
    m13.setNumBins(8 * 175)
    offset = 0
    for year in years:
        for mode in modes:
            print("offset", offset)
            signal_pdfs[f"{year}_{mode}"].setDataSize(
                datasets[f"{year}_{mode}"].getNumEvents(), 5, offset
            )
            offset = datasets[f"{year}_{mode}"].getNumEvents() + offset

    finalPDF_comb.setData(data_combined)

    fitman = FitManager(finalPDF_comb)
    fitman.setMaxCalls(320000)
    print("Running fit...")
    func_min = fitman.fit()

    params = {}
    for param in finalPDF_comb.getParameters():
        # params[param.getName()] = '{:.4f} +/- {:.4f}'.format(param.getValue(), param.getError())
        params[param.getName()] = param.getValue()
        params[param.getName() + "_err"] = param.getError()

    df = pd.DataFrame(params, index=[0])
    df.to_csv(f"{plotprefix}Fitted_params_timeintegrated_{mode}_{year}.csv")

    # reset cache counter before plotting
    Amp3Body.resetCacheCounter()
    for year in years:
        for mode in modes:
            print(f"Fit fractions {year} {mode}:")
            print(
                [
                    signal_pdfs[f"{year}_{mode}"].fit_fractions()[i][i]
                    for i in range(len(globaldecinfo.resonances))
                ]
            )
            print(
                f"Sum of fit fractions {year} {mode}: ",
                np.diag(signal_pdfs[f"{year}_{mode}"].fit_fractions()).sum(),
            )
            makePlots(
                total_pdfs[f"{year}_{mode}"],
                signal_pdfs[f"{year}_{mode}"],
                globaldecinfo,
                datasets[f"{year}_{mode}"],
                mode,
                year,
                plotprefix,
            )


def main():
    all_years = [2016, 2017, 2018]
    all_modes = [
        "SingleTag_D0ToKsPiPiLL",
        "SingleTag_D0ToKsPiPiDD",
        "DoubleTag_D0ToKsPiPiLL",
        "DoubleTag_D0ToKsPiPiDD",
    ]
    # perform single sample fits
    for year in all_years:
        for mode in all_modes:
            makeFit([year], [mode], "single/")
    # for i in range(0, 50):
    # for mode in all_modes:
    #  makeFit(all_years, [mode], f'combined/comb_{mode}_allyears_', i)
    # perform simultaneous fit to all samples
    makeFit(all_years, all_modes, f"combined/comb_allmodes_allyears_")


if __name__ == "__main__":
    sys.exit(main())
