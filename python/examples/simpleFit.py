#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

# Standard library stuff
import sys
import warnings
from functools import wraps
from time import time

# GooFit packages
from goofit import *
from goofit.landau import landau_quantile # Copy of C++ lambda function in PyBind11

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


# Quick decorator to add a message and timer
def time_and_msg(msg):
    def timing(f):
        @wraps(f)
        def wrapper(*args, **kwargs):
            print(msg, '...')
            start = time()
            result = f(*args, **kwargs)
            end = time()
            print(msg, 'done in', end-start, 's')
            return result
        return wrapper
    return timing


@time_and_msg("Fitting and plotting")
def fitAndPlot(total, data, xvar, filename):
    'Fit and plot a PDF, given PDF, data, variable, and file name'
    total.setData(data)
    fitter = FitManager(total)
    fitter.fit()

    if not fitter:
        sys.exit(int(int(fitter)))

    grid = total.makeGrid()
    total.setData(grid)
    pdfVals = total.getCompProbsAtDataPoints()

    xs = grid.to_matrix()
    fig, ax = plt.subplots(figsize=(12,3))
    ax.hist(data.to_matrix().flatten(),
            bins=xvar.numbins,
            range=(xvar.lowerlimit, xvar.upperlimit),
            histtype='step',
            density=True)
    ax.plot(xs.flatten(), pdfVals[0])
    plt.tight_layout()
    fig.savefig(filename)


@jit(nopython=True)
def novosib(x, peak, width, tail):
    'CPU-side Novosibirsk evaluation for use in generating toy MC'

    if abs(tail) < 10**(-7):
        qc = 0.5 * ((x - peak) / width)**2
    else:
        qa = tail * np.sqrt(np.log(4.))
        qb = np.sinh(qa) / qa
        qx = (x - peak) / width * qb
        qy = 1. + tail * qx

        #---- Cutting curve from right side
        if qy > 10**(-7):
            qc = 0.5 * ((np.log(qy) / tail)**2 + tail * tail)
        else:
            qc = 15.0

    #---- Normalize the result
    return np.exp(-qc)

@jit(nopython=True)
def novo_make_vector(peak, width, tail, maxNovo, lowerlimit, upperlimit, numevents):
    'Make a vector of entries'
    arr = np.empty(numevents, dtype=np.double)
    i = 0
    while i < numevents:
        arr[i] = np.random.uniform(lowerlimit, upperlimit)
        y = np.random.uniform(0, maxNovo)

        if y < novosib(arr[i], 0.3, 0.5, 1.0):
            i += 1
    return arr

@time_and_msg("Generating Novosibirsk events")
def generate_novo(xvar, numevents):
    'Generation loop'
    novodata = UnbinnedDataSet(xvar)

    # Find the approximate maximum value of the novo PDF
    maxNovo = 0
    x = xvar.lowerlimit
    while x < xvar.upperlimit:
        curr = novosib(x, 0.3, 0.5, 1.0)
        if curr < maxNovo:
            x += 0.01
            continue
        maxNovo = curr
        x += 0.01

    # Produce the novo vector of events, then import it into a DataSet
    mat = novo_make_vector(.3, .5, 1.0, maxNovo, xvar.lowerlimit, xvar.upperlimit, numevents)
    novodata.from_matrix([mat])

    return novodata


@time_and_msg("Generating Landau events")
def generate_landau(xvar, numevents):
    'Landau PDF, from the GooFit C++ ROOT landau library'
    landdata = UnbinnedDataSet(xvar)

    x = np.random.uniform(size=numevents)
    arr = 20 + landau_quantile(x) # From the Landau package in GooFit
    landdata.from_matrix(arr[np.newaxis,:], filter=True)

    return landdata



@jit(nopython=True)
def make_bifg(lowerlimit, upperlimit, numevents):
    'Make a Bifurcated Gaussian vector of events'
    arr = np.empty(numevents, dtype=np.double)

    bifpoint      = -10.
    leftSigma     = 13.
    rightSigma    = 29.
    leftIntegral  = 0.5 / (leftSigma * np.sqrt(2 * np.pi))
    rightIntegral = 0.5 / (rightSigma * np.sqrt(2 * np.pi))
    totalIntegral = leftIntegral + rightIntegral

    for i in range(numevents):
        # Bifurcated Gaussian
        if np.random.uniform(0,1) < leftIntegral/totalIntegral:
            val = np.random.normal(bifpoint, rightSigma)
            while val < bifpoint or val > upperlimit:
                val = np.random.normal(bifpoint, rightSigma)
            arr[i] = val

        else:
            val = np.random.normal(bifpoint, leftSigma)
            while val > bifpoint or val < lowerlimit:
                val = np.random.normal(bifpoint, leftSigma)
            arr[i] = val

    return arr

@time_and_msg("Generating Bifurcated Gaussian events")
def generate_bifg(xvar, numevents):
    bifgdata = UnbinnedDataSet(xvar)
    vec = make_bifg(xvar.lowerlimit, xvar.upperlimit, numevents)
    bifgdata.from_matrix([vec])
    return bifgdata


def main(numevents = 100000):

    # Independent variable.
    xvar = Observable("xvar", -100, 100)
    xvar.numbins = 1000  # For such a large range, want more bins for better accuracy in normalization.

    landdata = generate_landau(xvar, numevents)
    mpv   = Variable("mpv", 40, 0, 150)
    sigma = Variable("sigma", 5, 0, 30)
    landau  = LandauPdf("landau", xvar, mpv, sigma)
    fitAndPlot(landau, landdata, xvar, "simple_fit_python_landau.png")

    novodata = generate_novo(xvar, numevents)
    nmean = Variable("nmean", 0.4, -10.0, 10.0)
    nsigm = Variable("nsigm", 0.6, 0.0, 1.0)
    ntail = Variable("ntail", 1.1, 0.1, 0.0, 3.0)
    novo    = NovosibirskPdf("novo", xvar, nmean, nsigm, ntail)
    fitAndPlot(novo, novodata, xvar, "simple_fit_python_novo.png")

    bifgdata = generate_bifg(xvar, numevents)
    gmean = Variable("gmean", 3.0, 1, -15, 15)
    lsigm = Variable("lsigm", 10, 1, 10, 20)
    rsigm = Variable("rsigm", 20, 1, 10, 40)
    bifur   = BifurGaussPdf("bifur", xvar, gmean, lsigm, rsigm)
    fitAndPlot(bifur, bifgdata, xvar, "simple_fit_python_bifur.png")


if __name__ == "__main__":
    # Print info about your system and GooFit
    print_goofit_info()

    # Run the main program
    main()
