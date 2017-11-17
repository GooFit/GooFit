#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

from goofit import *
import numpy as np

print_goofit_info()

decayTime  = Observable("decayTime",0,10)

def integralExpCon(lo, hi):
    return(np.exp(-lo) - np.exp(-hi))

def integralExpLin(lo, hi):
    return((lo + 1) * np.exp(-lo) - (hi + 1) * np.exp(-hi))

def integralExpSqu(lo, hi):
    return((lo * lo + 2 * lo + 2) * np.exp(-lo) - (hi * hi + 2 * hi + 2) * np.exp(-hi))

def generateEvents(decayTime, conCoef, linCoef, squCoef, eventsToGenerate):

    totalRSintegral = integralExpCon(0, 100)
    step            = (decayTime.upperlimit - decayTime.lowerlimit) / decayTime.numbins


    rsEvtVec = np.array([])
    wsEvtVect = np.array([])

    for i in range(decayTime.numbins):

        binStart = i * step
        binStart += decayTime.lowerlimit
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


def fitRatio(decayTime, weights, rsEvts, wsEvts, plotName = ""):

    print("fitRatio")
    ratioData = BinnedDataSet(decayTime)

    i=0
    while i < wsEvts.size:

        ratio = wsEvts[i]

        if 0 == rsEvts[i]:
            rsEvts[i] = 1 #Cheating to avoid div by zero.

        ratio /= rsEvts[i]

        if 0 == wsEvts[i]:
            wsEvts[i] = 1 #Avoid zero errors

        error = wsEvts[i] / pow(rsEvts[i], 2)
        error += pow(wsEvts[i], 2) / pow(rsEvts[i], 3)
        error = np.sqrt(error)

        ratioData.setBinContent(i, ratio)
        ratioData.setBinError(i, error)
        ratioHist.SetBinContent(i + 1, ratio)
        ratioHist.SetBinError(i + 1, error)
        i+=1

    '''
    constaCoef = Variable("constaCoef", 0.03, 0.01, -1, 1)
    constaCoef.value = 0.03
    constaCoef.error=0.01
    linearCoef = Variable("linearCoef", 0, 0.01, -1, 1)
    linearCoef.value=0.00
    linearCoef.error=0.01
    secondCoef = Variable("secondCoef", 0, 0.01, -1, 1)
    secondCoef.value=0.00
    secondCoef.error=0.01


    weights = (constaCoef,linearCoef,secondCoef)
    '''

    poly = PolynomialPdf("poly", decayTime, weights)
    print("Setting binned error fit")
    poly.setFitControl(BinnedErrorFit())
    print("Setting data")
    poly.setData(ratioData)
    print("Setting FitManager")
    fitter = FitManager(poly)
    print("Fitting")
    fitter.fit()

    values = poly.evaluateAtPoints(decayTime)
    #pdfHist = TH1D("pdfHist", "", decayTime.numbins, decayTime.lowerlimit, decayTime.upperlimit)


    return datapdf

def cpvFitFcn(npar, gin, fun, fp, iflag):
    print(cpvFitFcn)
    conCoef = fp[0]
    linCoef = fp[1]
    squCoef = fp[2]

    chisq = 0
    step  = (decayTime.upperlimit - decayTime.lowerlimit) / decayTime.numbins

    i=0
    while i < ratios.size():
        currDTime = decayTime.lowerlimit + (i + 0.5) * step
        pdfval    = conCoef + linCoef * currDTime + squCoef * currDTime * currDTime
        chisq += pow((pdfval - ratios[i]) / errors[i], 2)
        i+=1

    fun = chisq


def fitRatioCPU(decayTime, rsEvts, wsEvts):
    print("fitRatioCPU")
    ratioHist = TH1D("ratioHist", "", decayTime.numbins, decayTime.lowerlimit, decayTime.upperlimit)

    ratios.resize(wsEvts.size())
    errors.resize(wsEvts.size())

    i=0
    while i < wsEvts.size():
        ratio = wsEvts[i]

        if 0 == rsEvts[i]:
            rsEvts[i] = 1 # Cheating to avoid div by zero.

        ratio /= rsEvts[i]

        if 0 == wsEvts[i]:
            wsEvts[i] = 1 # Avoid zero errors

        error = wsEvts[i] / pow(rsEvts[i], 2)
        error += pow(wsEvts[i], 2) / pow(rsEvts[i], 3)
        error = sqrt(error)

        ratios[i] = ratio
        errors[i] = error
        ratioHist.SetBinContent(i + 1, ratio)
        ratioHist.SetBinError(i + 1, error)
        i+=1

def main():
    print("main")
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

    dZeroEvtsRS, dZeroEvtsWS = generateEvents(decayTime, rSubD, dZeroLinearCoef, dZeroSecondCoef, eventsToGenerate)
    d0barEvtsRS, d0barEvtsWS = generateEvents(decayTime, rBarD, d0barLinearCoef, d0barSecondCoef, eventsToGenerate)

    constaCoef = Variable("constaCoef",0.03, 0.01, -1, 1)
    linearCoef = Variable("linearCoef",0, 0.01, -1, 1)
    secondCoef = Variable("secondCoef",0, 0.01, -1, 1)

    weights = (constaCoef, linearCoef, secondCoef)

    (retval1, fit1) = fitRatio(decayTime, weights, dZeroEvtsRS, dZeroEvtsWS, "dzeroEvtRatio.png")
    (retval2, fit2) = fitRatio(decayTime, weights, d0barEvtsRS, d0barEvtsWS, "dzbarEvtRatio.png")

    #CLI::Timer timer_cpu{"Total CPU (2x fits)"}
    fitRatioCPU(decayTime, dZeroEvtsRS, dZeroEvtsWS)
    fitRatioCPU(decayTime, d0barEvtsRS, d0barEvtsWS)
    cpu_string = timer_cpu.to_string();

    print (fit1, "\n", fit2, "\n", cpu_string)

    print("Exit codes (should be 0): {} and {}\n", retval1, retval2)

    return (retval1 + retval2)


    '''
    gpuTime = 0
    cpuTime = 0

    retval = fitRatio(dZeroEvtsRS, dZeroEvtsWS, "dzeroEvtRatio.png")
    if retval != 0:
        return retval

    timersub(stopTime, startTime, totalTime)
    gpuTime += totalTime.tv_sec + totalTime.tv_usec / 1000000.0
    retval = fitRatio(d0barEvtsRS, d0barEvtsWS, "dzbarEvtRatio.png")
    if retval != 0:
        return retval

    timersub(stopTime, startTime, totalTime)
    gpuTime += totalTime.tv_sec + totalTime.tv_usec / 1000000.0

    fitRatioCPU(dZeroEvtsRS, dZeroEvtsWS)
    timersub(stopTime, startTime, totalTime)
    cpuTime += totalTime.tv_sec + totalTime.tv_usec / 1000000.0
    fitRatioCPU(d0barEvtsRS, d0barEvtsWS)
    timersub(stopTime, startTime, totalTime)
    cpuTime += totalTime.tv_sec + totalTime.tv_usec / 1000000.0


    print( "GPU time [seconds] : ",gpuTime,"\nCPU time [seconds] : ",cpuTime)
    return 0;
    '''

main()
