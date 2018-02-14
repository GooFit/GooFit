#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

# Standard library stuff
import sys
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

# Optional import of numba, if available
# Make code much faster if found!
try:
    from numba import jit
except ImportError:
    def jit(*args, **kargs):
        def copyf(function):
            return function
        return copyf
    warnings.warn("Numba not found, will be 100x slower", RuntimeWarning)


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
    wsEvtVect = np.empty(numbins)

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


def produce_hist(rsEvts, wsEvents):
    ratio = wsEvts

    if rsEvts[0rsEvts == 0] = 1 # Cheating to avoid div by zero.
    ratio /= rsEvts
    wsEvts[wsEvts == 0] = 1 # Avoid zero errors

    error = np.sqrt(wsEvts / rsEvts**2 +  wsEvts**2 / rsEvts**3)
    return ratio, error


def fitRatio(decayTime, weights, rsEvts, wsEvts, plotName = ""):

    ratioData = BinnedDataSet(decayTime)

    ratio, error = produce_hist(rsEvts, wsEvts)


    poly = PolynomialPdf("poly", decayTime, weights)
    print("Setting binned error fit")
    poly.setFitControl(BinnedErrorFit())
    print("Setting data")
    poly.setData(ratioData)
    print("Setting FitManager")
    fitter = FitManager(poly)
    print("Fitting")
    fitter.fit()


    if plotname:
        values = poly.evaluateAtPoints(decayTime)
        xvals = decayTime.bin_centers
        fig, ax = plt.subplots()
        ax.errorbar(xvals, ratio, error)
        ax.plot(xvals, values)
        fig.savefig(plotname)

    return datapdf

def cpvFitFcn(numbins, lowerlimit, upperlimit, npar, gin, fun, fp, iflag):
    print(cpvFitFcn)
    conCoef = fp[0]
    linCoef = fp[1]
    squCoef = fp[2]

    chisq = 0
    step  = (upperlimit - lowerlimit) / numbins

    i=0
    while i < ratios.size():
        currDTime = lowerlimit + (i + 0.5) * step
        pdfval    = conCoef + linCoef * currDTime + squCoef * currDTime * currDTime
        chisq += pow((pdfval - ratios[i]) / errors[i], 2)
        i+=1

    fun = chisq

@jit(nopython=True)
def fitRatioCPU(numbins, lowerlimit, upperlimit, rsEvts, wsEvts):

    ratios.empty_like(wsEvts)
    errors.empty_like(wsEvts)

    for i in range(len(wsEvts)):
        ratio = wsEvts[i]

        if 0 == rsEvts[i]:
            rsEvts[i] = 1 # Cheating to avoid div by zero.

        ratio /= rsEvts[i]

        if 0 == wsEvts[i]:
            wsEvts[i] = 1 # Avoid zero errors

        error = wsEvts[i] / rsEvts[i]**2
        error += wsEvts[i]**2 / rsEvts[i]**3
        error = np.sqrt(error)

        ratios[i] = ratio
        errors[i] = error

    return ratios, errors


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
            decayTime.lowerlimit, decayTime.upperlimit, decayTime.numbins,
            rSubD, dZeroLinearCoef, dZeroSecondCoef, eventsToGenerate)

    d0barEvtsRS, d0barEvtsWS = generateEvents(
            decayTime.lowerlimit, decayTime.upperlimit, decayTime.numbins,
            rBarD, d0barLinearCoef, d0barSecondCoef, eventsToGenerate)

    constaCoef = Variable("constaCoef", 0.03, 0.01, -1, 1)
    linearCoef = Variable("linearCoef", 0, 0.01, -1, 1)
    secondCoef = Variable("secondCoef", 0, 0.01, -1, 1)

    weights = (constaCoef, linearCoef, secondCoef)

    (retval1, fit1) = fitRatio(decayTime, weights, dZeroEvtsRS, dZeroEvtsWS, "dzeroEvtRatio.png")
    (retval2, fit2) = fitRatio(decayTime, weights, d0barEvtsRS, d0barEvtsWS, "dzbarEvtRatio.png")

    with timed("Total CPU (2x fits)"):
        fitRatioCPU(decayTime, dZeroEvtsRS, dZeroEvtsWS)
        fitRatioCPU(decayTime, d0barEvtsRS, d0barEvtsWS)

    print(fit1)
    print(fit2)
    print("Exit codes (should be 0):", retval1, "and", retval2)

    return (retval1 + retval2)



if __name__ == "__main__":
    # Print info about your system and GooFit
    print_goofit_info()

    # Run the main program
    main()
