#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

# Standard library stuff
import re
import warnings
from time import time
from contextlib import contextmanager

# GooFit package
from goofit import *

import numpy as np

# Protect the matplotlib call for systems with no graphics
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


try:
    import uncertainties as unc
    from uncertainties import unumpy
except ImportError:
    unc = None
    warnings.warn("uncertainties not found, no unc on plot", RuntimeWarning)

# Optional import of numba, if available
# Make code much faster if found!
try:
    from numba import jit
except ImportError:
    def jit(*args, **kargs):
        def copyf(function):
            return function
        return copyf
    warnings.warn("Numba not found, will be 100x slower. Try `pip install numba` to install.", RuntimeWarning)

# Simple timer in a context manager
@contextmanager
def timed(msg):
    start = time()
    try:
        yield
    finally:
        end = time()
        print(msg, end - start, "s")

@jit(nopython=True)
def integralExpCon(lo, hi):
    return(np.exp(-lo) - np.exp(-hi))

@jit(nopython=True)
def integralExpLin(lo, hi):
    return((lo + 1) * np.exp(-lo) - (hi + 1) * np.exp(-hi))

@jit(nopython=True)
def integralExpSqu(lo, hi):
    return((lo * lo + 2 * lo + 2) * np.exp(-lo) - (hi * hi + 2 * hi + 2) * np.exp(-hi))

@jit(nopython=True)
def generateEvents(numbins, lowerlimit, upperlimit, conCoef, linCoef, squCoef, eventsToGenerate):

    totalRSintegral = integralExpCon(0, 100)
    step            = (upperlimit - lowerlimit) / numbins

    rsEvtVec = np.empty(numbins)
    wsEvtVec = np.empty(numbins)

    for i in range(numbins):

        binStart = i * step
        binStart += lowerlimit
        binFinal = binStart + step

        rsIntegral = integralExpCon(binStart, binFinal)
        wsIntegral = conCoef * integralExpCon(binStart, binFinal)
        wsIntegral += linCoef * integralExpLin(binStart, binFinal)
        wsIntegral += squCoef * integralExpSqu(binStart, binFinal)

        expectedRSevts = eventsToGenerate * rsIntegral / totalRSintegral
        expectedWSevts = eventsToGenerate * wsIntegral / totalRSintegral

        rsEvtVec[i] = np.random.poisson(expectedRSevts)
        wsEvtVec[i] = np.random.poisson(expectedWSevts)

    return rsEvtVec, wsEvtVec


def produce_hist(rsEvts, wsEvts):

    rsEvts[rsEvts == 0] = 1 # Cheating to avoid div by zero.
    ratio = wsEvts / rsEvts
    wsEvts[wsEvts == 0] = 1 # Avoid zero errors

    error = np.sqrt(wsEvts / rsEvts**2 +  wsEvts**2 / rsEvts**3)
    return ratio, error


def fitRatio(decayTime, weights, rsEvts, wsEvts, plotname=None):

    ratioData = BinnedDataSet(decayTime)

    ratio, error = produce_hist(rsEvts, wsEvts)

    for i in range(len(ratio)):
        ratioData.setBinContent(i, ratio[i])
        ratioData.setBinError(i, error[i])

    poly = PolynomialPdf("poly", decayTime, weights)
    poly.setFitControl(BinnedErrorFit())
    poly.setData(ratioData)
    fitter = FitManager(poly)
    fitter.verbosity = 0
    mn_func_min = fitter.fit()
    mn_cov = mn_func_min.UserCovariance()
    cov = mn_cov.to_matrix()
    print(cov)

    xvals = np.array(decayTime.bin_centers)
    p = [w.value for w in reversed(weights)]
    print("GooFit fit:", end=' ')
    plot_ratio(xvals, ratio, error, p, plotname, cov)

    # This would work just fine, too
    # values = poly.evaluateAtPoints(decayTime)
    # plot.plot(xvals, values)

    return fitter


def fitRatioCPU(decayTime, rsEvts, wsEvts, plotname=None):
    ratio, error = produce_hist(rsEvts, wsEvts)
    xvals = np.array(decayTime.bin_centers)

    p, cov = np.polyfit(xvals, ratio, 2, w=1/error, cov=True)

    print("CPU Fit:", end=' ')
    plot_ratio(xvals, ratio, error, p, plotname, cov)

def plot_ratio(xvals, ratio, error, p, plotname, cov=None):
    result = "{:.2e} t^2 + {:.2e} t + {:.2e}".format(*p)
    result = result.replace('e-0','e-')
    print("Fit result:", result)

    if plotname:
        result = re.sub(r'e(\S+?)(\s|$)', r' \\times {10}^{\1} ', result)
        fig, ax = plt.subplots()
        ax.errorbar(xvals, ratio, yerr=error, fmt='k.')
        ax.plot(xvals, np.polyval(p, xvals))
        ax.text(.1, .85, "${}$".format(result),
                transform=ax.transAxes, fontsize=16)

        if cov is not None and unc is not None:

            p_unc = unc.correlated_values(p, cov)
            y = np.polyval(p_unc, xvals)
            ax.fill_between(xvals,
                     unumpy.nominal_values(y)-unumpy.std_devs(y),
                     unumpy.nominal_values(y)+unumpy.std_devs(y),
                     alpha=.5, facecolor='b')

        fig.savefig(plotname)


def main():
    decayTime = Observable("decayTime", 0, 10)

    numbins = 100

    decayTime.value = 100
    decayTime.numbins = numbins

    rSubD = 0.03
    rBarD = 0.03
    delta = 0
    wpPhi = 0
    x_mix = 0.0016
    y_mix = 0.0055
    magPQ = 1.0
    magQP = 1.0 / magPQ

    eventsToGenerate = 10000000

    dZeroLinearCoef = magPQ * np.sqrt(rSubD) * (y_mix * np.cos(delta + wpPhi) - x_mix * np.sin(delta + wpPhi))
    d0barLinearCoef = magQP * np.sqrt(rBarD) * (y_mix * np.cos(delta - wpPhi) - x_mix * np.sin(delta - wpPhi))

    dZeroSecondCoef = 0.25 * magPQ * magPQ * (x_mix * x_mix + y_mix * y_mix)
    d0barSecondCoef = 0.25 * magQP * magQP * (x_mix * x_mix + y_mix * y_mix)

    dZeroEvtsRS, dZeroEvtsWS = generateEvents(
            decayTime.numbins, decayTime.lowerlimit, decayTime.upperlimit,
            rSubD, dZeroLinearCoef, dZeroSecondCoef, eventsToGenerate)

    d0barEvtsRS, d0barEvtsWS = generateEvents(
            decayTime.numbins, decayTime.lowerlimit, decayTime.upperlimit,
            rBarD, d0barLinearCoef, d0barSecondCoef, eventsToGenerate)

    constaCoef = Variable("constaCoef", 0.03, 0.01, -1, 1)
    linearCoef = Variable("linearCoef", 0, 0.01, -1, 1)
    secondCoef = Variable("secondCoef", 0, 0.01, -1, 1)

    weights = (constaCoef, linearCoef, secondCoef)

    retval1 = fitRatio(decayTime, weights, dZeroEvtsRS, dZeroEvtsWS, "chisquare_dzeroEvtRatio_goo_python.png")
    retval2 = fitRatio(decayTime, weights, d0barEvtsRS, d0barEvtsWS, "chisquare_dzbarEvtRatio_goo_python.png")

    with timed("Total CPU (2x fits)"):
        fitRatioCPU(decayTime, dZeroEvtsRS, dZeroEvtsWS, "chisquare_dzeroEvtRatio_cpu_python.png")
        fitRatioCPU(decayTime, d0barEvtsRS, d0barEvtsWS, "chisquare_dzbarEvtRatio_cpu_python.png")

    print("Exit codes (should be 0):", int(retval1), "and", int(retval2))

    return (int(retval1) + int(retval2))



if __name__ == "__main__":
    # Print info about your system and GooFit
    print_goofit_info()

    # Run the main program
    main()
