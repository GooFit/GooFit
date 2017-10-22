#!/usr/bin/env python
# -*- coding: utf-8 -*-

from goofit import *
import math
import numpy as np
import pylandau

# CPU-side Novosibirsk evaluation for use in generating toy MC.
def novosib(x,peak,width,tail):
    print("novosib")
    qa = 0
    qb = 0
    qc = 0
    qx = 0
    qy = 0

    if abs(tail) < math.pow(10,-7):
        qc = 0.5 * math.pow(((x - peak) / width), 2)
    else:
        qa = tail * math.sqrt(math.log(4.))
        qb = math.sinh(qa) / qa
        qx = (x - peak) / width * qb
        qy = 1. + tail * qx

        #---- Cutting curve from right side

        if qy > math.pow(10,-7):
            qc = 0.5 * (math.pow((math.log(qy) / tail), 2) + tail * tail)
        else:
            qc = 15.0

    #---- Normalize the result

    return math.exp(-qc)



def fitAndPlot(total,data,dataHist,xvar,fname):
    print("firAndPlot")
    total.setData(data)
    fitter = FitManager(total)
    fitter.fit()

    if not fitter:
        exit(fitter)

    grid = UnbinnedDataSet(xvar)
    step = (xvar.upperlimit - xvar.lowerlimit) / xvar.numbins
    i = 0
    while i < xvar.numbins:
        xvar.value(xvar.lowerlimit + (i + 0.5) * step)
        grid.addEvent()
        i+=1


    total.setData(grid)
    pdfVals = total.getCompProbsAtDataPoints()

    totalPdf = 0

    i = 0
    while i < grid.getNumEvents():
        grid.loadEvent(i)
        pdfHist.Fill(xvar.getValue(), pdfVals[0][i])
        totalPdf += pdfVals[0][i]
        i+=1

    i = 0
    while i < xvar.numbins:
        val = pdfHist.GetBinContent(i + 1)
        val /= totalPdf
        val *= data.getNumEvents()
        pdfHist.SetBinContent(i + 1, val)
        i+=1


def main():
    print("main")
    numevents = 100000

    # Independent variable.
    xvar = Variable("xvar", -100, 100)
    xvar.numbins = 1000  # For such a large range, want more bins for better accuracy in normalisation.

    # Data sets for the three fits.
    landdata = UnbinnedDataSet(xvar)
    bifgdata = UnbinnedDataSet(xvar)
    novodata = UnbinnedDataSet(xvar)



    maxNovo = 0
    x = xvar.lowerlimit

    while x < xvar.upperlimit:
        print("main")
        curr = novosib(x, 0.3, 0.5, 1.0)
        if curr < maxNovo:
            x += 0.01
            continue
        maxNovo = curr
        x += 0.01

    leftSigma     = 13
    rightSigma    = 29
    leftIntegral  = 0.5 / (leftSigma * math.sqrt(2 * math.pi))
    rightIntegral = 0.5 / (rightSigma * math.sqrt(2 * math.pi))
    totalIntegral = leftIntegral + rightIntegral
    bifpoint      = -10

    # Generating three sets of toy MC.
    while landdata.getNumEvents() < numevents:
        # Landau
        x = np.random.normal(20,1)
        y_landau = pylandau.landau(x)

        landdata.addEvent()


    while bifgdata.getNumEvents() < numevents:
        # Bifurcated Gaussian
        val
        if(random.uniform() < (leftIntegral / totalIntegral)):
            val = random.gausaus(bifpoint, rightSigma)
            while val < bifpoint or val > xvar.upperlimit:
                val = random.gaus(bifpoint, rightSigma)
            xvar.value(val)

        else:
            val = random.gaus(bifpoint, leftSigma)
            while val > bifpoint or val < xvar.lowerlimit:
                val = random.gaus(bifpoint, leftSigma)
            xvar.value(val)

        bifgdata.addEvent()


    while novodata.getNumEvents() < numevents:
        # And Novosibirsk.
        while True:
            print("main2")
            xvar.value(random.uniform(xvar.lowerlimit, xvar.upperlimit))
            y = random.uniform(0, maxNovo)

            if y < novosib(xvar.getValue(), 0.3, 0.5, 1.0):
                break


        novodata.addEvent()

    mpv   = Variable("mpv", 40, 0, 150)
    sigma = Variable("sigma", 5, 0, 30)
    landau  = LandauPdf("landau", xvar, mpv, sigma)
    fitAndPlot(landau, landdata, landHist, xvar, "landau.png")

    nmean = Variable("nmean", 0.4, -10.0, 10.0)
    nsigm = Variable("nsigm", 0.6, 0.0, 1.0)
    ntail = Variable("ntail", 1.1, 0.1, 0.0, 3.0)
    novo    = NovosibirskPdf("novo", xvar, nmean, nsigm, ntail)
    fitAndPlot(novo, novodata, novoHist, xvar, "novo.png")

    gmean = Variable("gmean", 3.0, 1, -15, 15)
    lsigm = Variable("lsigm", 10, 1, 10, 20)
    rsigm = Variable("rsigm", 20, 1, 10, 40)
    bifur   = BifurGaussPdf("bifur", xvar, gmean, lsigm, rsigm)
    fitAndPlot(bifur, bifgdata, bifgHist, xvar, "bifur.png")


main()
