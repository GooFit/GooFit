from goofit import *
import sys
import os
from math import sqrt
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import argparse
from distutils.util import strtobool
from root_pandas import to_root

print_goofit_info()

def makeSignalPdf(m12, m13, dtime, sigma, mistag, eventNumber, resolution):
    # Now the amplitudes, masses and widths are fixed.
    #tau = Variable("tau", 0.4101, 0.001, 0.300, 0.500)
    xmixing = Variable("xmixing", 0.0016, 0.001, 0, 0)
    ymixing = Variable("ymixing", 0.0055, 0.001, 0, 0)
    tau = Variable("tau", 0.4101)
    #xmixing = Variable("xmixing", 0)
    #ymixing = Variable("ymixing", 0)
    d0toks0pipi = DecayInfo3t(tau, xmixing, ymixing)
    
    d0toks0pipi.motherMass = 1.86484
    d0toks0pipi.daug1Mass = 0.497614
    d0toks0pipi.daug2Mass = 0.139570
    d0toks0pipi.daug3Mass = 0.139570
    d0toks0pipi.meson_radius = 1.4

    # Gounaris Sakurai (name, amplitude real, amplitude imaginary, mass, width, spin, cyclic index)
    rho_770 = Resonances.GS("rho_770", 
                            Variable("rho_770_amp_real", 1),
                            Variable("rho_770_amp_imag", 0),
                            Variable("rho_770_mass", 0.77526),
                            Variable("rho_770_width", 0.1478),
                            1,
                            PAIR_23)

    # Relativistic Breit Wigner (name, amplitude real, amplitude imaginary, mass, width, spin, cyclic index)
    f2_1270 = Resonances.RBW("f2_1270",
                             Variable("f2_1270_amp_real", 1.152477444, 0.0001, 0, 0),
                             Variable("f2_1270_amp_imag", -0.8465788457, 0.0001, 0, 0),
                             Variable("f2_1270_mass", 1.2755),
                             Variable("f2_1270_width", 0.1867),
                             2,
                             PAIR_23)

    omega_782 = Resonances.RBW("omega_782",
                               Variable("omega_782_amp_real", -0.01980906521, 0.0001, 0, 0),
                               Variable("omega_782_amp_imag", 0.03336226814, 0.0001, 0, 0),
                               Variable("omega_782_mass", 0.78265),
                               Variable("omega_782_width", 0.00849),
                               1,
                               PAIR_23)

    rho_1450 = Resonances.RBW("rho_1450",
                              Variable("rho_1450_amp_real", -0.5974129043, 0.0001, 0, 0),
                              Variable("rho_1450_amp_imag", 2.786682225, 0.0001, 0, 0),
                              Variable("rho_1450_mass", 1.465),
                              Variable("rho_1450_width", 0.400),
                              1,
                              PAIR_23)

    # KS0 pi- resonances 
    Kstarm_892 = Resonances.RBW("Kstarm_892",
                                Variable("Kstarm_892_amp_real", -1.253826039, 0.0001, 0, 0),
                                Variable("Kstarm_892_amp_imag", 1.177421022, 0.0001, 0, 0),
                                Variable("Kstarm_892_mass", 0.89176),
                                Variable("Kstarm_892_width", 0.0503),
                                1,
                                PAIR_13)

    K2starm_1430 = Resonances.RBW("K2starm_1430",
                                  Variable("K2starm_1430_amp_real", 0.9120203982, 0.0001, 0, 0),
                                  Variable("K2starm_1430_amp_imag", -0.8838092517, 0.0001, 0, 0),
                                  Variable("K2starm_1430_mass", 1.4256),
                                  Variable("K2starm_1430_width", 0.0985),
                                  2,
                                  PAIR_13)

    Kstarm_1680 = Resonances.RBW("Kstarm_1680",
                                 Variable("Kstarm_1680_amp_real", -1.564143032, 0.0001, 0, 0),
                                 Variable("Kstarm_1680_amp_imag", -2.917114426, 0.0001, 0, 0),
                                 Variable("Kstarm_1680_mass", 1.718),
                                 Variable("Kstarm_1680_width", 0.322),
                                 1,
                                 PAIR_13)

    Kstarm_1410 = Resonances.RBW("Kstarm_1410",
                                 Variable("Kstarm_1410_amp_real", -0.04736452905, 0.0001, 0, 0),
                                 Variable("Kstarm_1410_amp_imag", 0.2861059269, 0.0001, 0, 0),
                                 Variable("Kstarm_1410_mass", 1.414),
                                 Variable("Kstarm_1410_width", 0.232),
                                 1,
                                 PAIR_13)

    # KS0 pi+ resonances 
    Kstarp_892 = Resonances.RBW("Kstarp_892",
                                Variable("Kstarp_892_amp_real", 0.1214919538, 0.0001, 0, 0),
                                Variable("Kstarp_892_amp_imag", -0.1101621766, 0.0001, 0, 0),
                                Variable("Kstarp_892_mass", 0.89176),
                                Variable("Kstarp_892_width", 0.503),
                                1,
                                PAIR_12)

    K2starp_1430 = Resonances.RBW("K2starp_1430",
                                  Variable("K2starp_1430_amp_real", 0.0006981260298, 0.0001, 0, 0),
                                  Variable("K2starp_1430_amp_imag", -0.09999756307, 0.0001, 0, 0),
                                  Variable("K2starp_1430_mass", 1.4256),
                                  Variable("K2starp_1430_width", 0.0985),
                                  2,
                                  PAIR_12)

    Kstarp_1410 = Resonances.RBW("Kstarp_1410",
                                  Variable("Kstarp_1410_amp_real", -0.1822307452, 0.0001, 0, 0),
                                  Variable("Kstarp_1410_amp_imag", 0.1043645318, 0.0001, 0, 0),
                                  Variable("Kstarp_1410_mass", 1.414),
                                  Variable("Kstarp_1410_width", 0.232),
                                  1,
                                  PAIR_12)

    # LASS: Kpi S wave
    K0star_1430_mass = Variable("K0starp_1430_mass", 1.441)
    K0star_1430_width = Variable("K0starp_1430_width", 0.193)
    a = Variable("a", 0.113)
    r = Variable("r", -33.8)
    R = Variable("R", 1) # Fixed
    phiR = Variable("phiR", -1.91462619)
    F = Variable("F", 0.96)
    phiF = Variable("phiF", 0.00174533)
    K0starp_1430 = Resonances.LASS("K0starp_1430",
                                   Variable("K0starp_1430_amp_real", -0.1047927629, 0.0001, 0, 0),
                                   Variable("K0starp_1430_amp_imag", 0.0334436367, 0.0001, 0, 0),
                                   K0star_1430_mass,
                                   K0star_1430_width,
                                   a,
                                   r,
                                   R,
                                   phiR,
                                   F,
                                   phiF,
                                   0,
                                   PAIR_12)

    K0starm_1430 = Resonances.LASS("K0starm_1430",
                                   Variable("K0starm_1430_amp_real", -0.3854492709, 0.0001, 0, 0),
                                   Variable("K0starm_1430_amp_imag", 2.328310301, 0.0001, 0, 0),
                                   K0star_1430_mass,
                                   K0star_1430_width,
                                   a,
                                   r,
                                   R,
                                   phiR,
                                   F,
                                   phiF,
                                   0,
                                   PAIR_13)

    d0toks0pipi.resonances = (rho_770, f2_1270, omega_782, rho_1450, Kstarm_892, K2starm_1430, Kstarm_1680, Kstarm_1410, 
                              Kstarp_892, K2starp_1430, Kstarp_1410,
                              K0starp_1430, K0starm_1430)
    # Constant efficiency
    constantOne = Variable('One', 1)
    constantZero = Variable('Zero', 0)
    eff = PolynomialPdf('constantEff', [m12, m13], [constantOne], [constantZero,constantZero], 0)
    
    # Amp3Body_TD: name, dtime, dtime_err, m12, m13, eventNumber, decay, mixing resolution, eff
    d = Amp3Body_TD("signalPDF", dtime, sigma, m12, m13, eventNumber, d0toks0pipi, resolution, eff, mistag)
    return d

def inDalitz(m12, m13):
    mD0, mKS, mPi = 1.86486, 0.497614, 0.13957

    if m12 < pow(mKS + mPi, 2) or m12 > pow(mD0 - mPi, 2):
        return False

    # Calculate energies of particles 1 and 3 in m12 rest frame
    E1 = 0.5 * (m12 - mPi*mPi + mKS*mKS) / sqrt(m12)
    E3 = 0.5 * (mD0*mD0 - m12 - mPi*mPi) / sqrt(m12)    

    minimum = pow(E1 + E3, 2) - pow(sqrt(E1*E1 - mKS*mKS) + sqrt(E3*E3 - mPi*mPi), 2)
    maximum = pow(E1 + E3, 2) - pow(sqrt(E1*E1 - mKS*mKS) - sqrt(E3*E3 - mPi*mPi), 2)
    if m13 < minimum or m13 > maximum:
        return False

    return True

def Getm23(m12, m13):
    mD0, mKS, mPi = 1.86486, 0.497614, 0.13957
    m23 = mD0**2 + mKS**2 + 2*mPi**2 - m12 - m13
    return m23

def rebin2D(arr, new_shape):
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)

def makePlots(signal, dataset, dtime, sigma, mistag, m12, m13, eventNumber):
    m12.setNumBins(m12.getNumBins()*5)
    m13.setNumBins(m13.getNumBins()*5)

    data_tmp = UnbinnedDataSet(dtime, sigma, mistag, m12, m13, eventNumber)
    data_dtime = UnbinnedDataSet(dtime, sigma, mistag, m12, m13, eventNumber)

    count = 0
    sigma.setValue(0.1)
    mistag.setValue(0)
    m12val, m13val, m23val, dtval = [], [], [], []

    for i in range(m12.getNumBins()):
        m12.setValue(m12.getLowerLimit() + m12.getBinSize() * (i + 0.5))
        for j in range(m13.getNumBins()):
            m13.setValue(m13.getLowerLimit() + m13.getBinSize() * (j + 0.5))
            if inDalitz(m12.getValue(), m13.getValue()):
                for k in range(10):
                    dtime.setValue(dtime.getLowerLimit() + (dtime.getUpperLimit() - dtime.getLowerLimit()) * (k + 0.5) / 10)
                    m12val.append(m12.getValue())
                    m13val.append(m13.getValue())
                    m23val.append(Getm23(m12.getValue(), m13.getValue()))
                    eventNumber.setValue(count)
                    data_tmp.addEvent()
                    count += 1

    data_tmp.setValueForAllEvents(sigma)
    data_tmp.setValueForAllEvents(mistag)
    print(data_tmp.getNumEvents())

    print(data_tmp.getObservables())
    print(signal.getObservables())

    signal.setDataSize(data_tmp.getNumEvents(), 6)

    prod = ProdPdf('prod', [signal])

    prod.setData(data_tmp)
    prod.normalize()
    vals = prod.getCompProbsAtDataPoints()
    array = np.array(vals)
    print(array.shape)

    m12dat = np.array(dataset.to_numpy())[0]
    m13dat = np.array(dataset.to_numpy())[1]
    dtdat  = np.array(dataset.to_numpy())[2]
    m23dat = Getm23(m12dat, m13dat)

    _, ax = plt.subplots(1, 2, figsize=(10, 5), num=1)
    hist_dat, _, _, _ = ax[0].hist2d(m12dat, m13dat, bins=175, range=[(0,3.5),(0,3.5)], cmap='hot_r')
    hist_fit, _, _, _ = ax[1].hist2d(m12val, m13val, weights=array[0], bins=175, range=[(0,3.5),(0,3.5)], density=True, cmap='hot_r')

    for a in ax:
        a.set_xlabel('m12')
        a.set_ylabel('m13')
        a.axis('equal')
    ax[0].set_title('Toy Data')
    ax[1].set_title('Fitted Function')
    plt.tight_layout()
    plt.savefig('Tddp_dalitz.pdf')

    pulls = (hist_fit/hist_fit.sum() - hist_dat/hist_dat.sum()) / (np.sqrt(hist_dat)/hist_dat.sum())
    print(pulls)
    plt.figure(2)
    image = plt.imshow(pulls, origin='lower', extent=(0,3.5,0,3.5), vmin=-4, vmax=4, cmap='RdYlBu')
    plt.xlabel('m12')
    plt.ylabel('m13')
    plt.title('Pulls')
    plt.colorbar(image)
    plt.savefig('Tddp_pulls.pdf')

    # Plot Dalitz variables 
    hist, _ = np.histogram(m12val, weights=array[0], bins=175*5, range=(0,3.5), density=True)
    hist = np.reshape(hist, (-1, 5)).mean(axis=1)
    bin_edges = np.linspace(0,3.5,176)
    bin_centre = (bin_edges[:-1] + bin_edges[1:])/2
    
    hist_dat, _ = np.histogram(m12dat, bins=175*5, range=(0,3.5), density=True)
    hist_tmp, _ = np.histogram(m12dat, bins=175*5, range=(0,3.5))
    hist_dat = np.reshape(hist_dat, (-1, 5)).mean(axis=1)
    hist_tmp = np.reshape(hist_tmp, (-1, 5)).mean(axis=1)

    yerr = np.sqrt(hist_tmp)/(hist_tmp.sum()*bin_edges[1]) # sqrt(n) / integral
    pulls = (hist - hist_dat)/yerr

    _, ax = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[3, 1]}, num=4)
    ax[0].plot(bin_centre, hist, zorder=2, color='b')
    ax[0].errorbar(bin_centre, hist_dat, yerr=yerr, capsize=2, fmt='o', color='k', ms=4, zorder=1, elinewidth=1)
    ax[0].set_xlim(0,3.5)
    ax[0].set_ylabel('Candidates (normalised)')
    ax[1].bar(bin_centre, pulls, color='k', width=bin_edges[1])
    ax[1].set_xlim(0,3.5)
    ax[1].set_xlabel('m12')
    ax[1].set_ylabel('Pulls')
    plt.tight_layout()
    plt.savefig('m12_test.pdf')

    hist, _ = np.histogram(m13val, weights=array[0], bins=175*5, range=(0,3.5), density=True)
    hist_dat, _ = np.histogram(m13dat, bins=175*5, range=(0,3.5), density=True)
    hist_tmp, _ = np.histogram(m13dat, bins=175*5, range=(0,3.5))
    hist = np.reshape(hist, (-1, 5)).mean(axis=1)
    hist_dat = np.reshape(hist_dat, (-1, 5)).mean(axis=1)
    hist_tmp = np.reshape(hist_tmp, (-1, 5)).mean(axis=1)
    yerr = np.sqrt(hist_tmp)/(hist_tmp.sum()*bin_edges[1]) # sqrt(n) / integral
    pulls = (hist - hist_dat)/yerr

    _, ax = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[3, 1]}, num=5)
    ax[0].plot(bin_centre, hist, zorder=2, color='b')
    ax[0].errorbar(bin_centre, hist_dat, yerr=yerr, capsize=2, fmt='o', color='k', ms=4, zorder=1, elinewidth=1)
    ax[0].set_xlim(0,3.5)
    ax[0].set_ylabel('Candidates (normalised)')
    ax[1].bar(bin_centre, pulls, color='k', width=bin_edges[1])
    ax[1].set_xlim(0,3.5)
    ax[1].set_xlabel('m13')
    ax[1].set_ylabel('Pulls')
    plt.tight_layout()
    plt.savefig('m13_test.pdf')

    hist, _ = np.histogram(m23val, weights=array[0], bins=175*5, range=(0,3.5), density=True)
    bin_centre = (bin_edges[:-1] + bin_edges[1:])/2
    hist_dat, _ = np.histogram(m23dat, bins=175*5, range=(0,3.5), density=True)
    hist_tmp, _ = np.histogram(m23dat, bins=175*5, range=(0,3.5))
    hist = np.reshape(hist, (-1, 5)).mean(axis=1)
    hist_dat = np.reshape(hist_dat, (-1, 5)).mean(axis=1)
    hist_tmp = np.reshape(hist_tmp, (-1, 5)).mean(axis=1)
    yerr = np.sqrt(hist_tmp)/(hist_tmp.sum()*bin_edges[1]) # sqrt(n) / integral
    pulls = np.nan_to_num((hist - hist_dat)/yerr)
    print(hist)
    print(hist_dat)
    print(pulls)

    _, ax = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[3, 1]}, num=6)
    ax[0].plot(bin_centre, hist, zorder=2, color='b')
    ax[0].errorbar(bin_centre, hist_dat, yerr=yerr, capsize=2, fmt='o', color='k', ms=4, zorder=1, elinewidth=1)
    ax[0].set_xlim(0,2)
    ax[0].set_ylabel('Candidates (normalised)')
    ax[1].bar(bin_centre, pulls, color='k', width=bin_edges[1])
    ax[1].set_xlim(0,2)
    ax[1].set_xlabel('m23')
    ax[1].set_ylabel('Pulls')
    plt.tight_layout()
    plt.savefig('m23_test.pdf')

    # For decay time plots
    count = 0
    sigma.setValue(0.1)
    mistag.setValue(0)

    for i in range(80):
        m12.setValue(m12.getLowerLimit() + (m12.getUpperLimit() - m12.getLowerLimit()) * (i + 0.5) / 80)
        for j in range(80):
            m13.setValue(m13.getLowerLimit() + (m13.getUpperLimit() - m13.getLowerLimit()) * (j + 0.5) / 80)
            if inDalitz(m12.getValue(), m13.getValue()):
                for k in range(dtime.getNumBins()):
                    dtime.setValue(dtime.getLowerLimit() + dtime.getBinSize() * (k + 0.5))
                    dtval.append(dtime.getValue())
                    eventNumber.setValue(count)
                    data_dtime.addEvent()
                    count += 1

    data_dtime.setValueForAllEvents(sigma)
    data_dtime.setValueForAllEvents(mistag)
    print(data_dtime.getNumEvents())

    print(data_dtime.getObservables())
    print(signal.getObservables())

    signal.setDataSize(data_dtime.getNumEvents(), 6)

    prod.setData(data_dtime)
    prod.normalize()
    vals2 = prod.getCompProbsAtDataPoints()
    array = np.array(vals2)
    print(array.shape)

    # Plot Decay time 
    hist, bin_edges = np.histogram(dtval, weights=array[0], bins=175, range=(0,6), density=True)
    bin_centre = (bin_edges[:-1] + bin_edges[1:])/2
    hist_dat, _ = np.histogram(dtdat, bins=175, range=(0,6), density=True)
    hist_tmp, _ = np.histogram(dtdat, bins=175, range=(0,6))
    yerr = np.sqrt(hist_tmp)/(hist_tmp.sum()*bin_edges[1]) # sqrt(n) / integral
    pulls = (hist - hist_dat)/yerr

    _, ax = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[3, 1]}, num=3)
    ax[0].plot(bin_centre, hist, zorder=2, color='b')
    ax[0].errorbar(bin_centre, hist_dat, yerr=yerr, capsize=2, fmt='o', color='k', ms=4, zorder=1, elinewidth=1)
    ax[0].set_xlim(0,6)
    ax[0].set_ylabel('Candidates (normalised)')
    ax[1].bar(bin_centre, pulls, color='k', width=bin_edges[1])
    ax[1].set_xlim(0,6)
    ax[1].set_xlabel('Decay Time (ps)')
    ax[1].set_ylabel('Pulls')
    plt.tight_layout()
    plt.savefig('dtime_test.pdf')

def fillToyDataSet(nEvents, dataset, signal, m12, m13, dtime, sigma, mistag, eventNumber):
    m12.setNumBins(175*2)
    m13.setNumBins(175*2)
    dtime.setNumBins(175)
    data_tmp = UnbinnedDataSet(m12, m13, dtime, sigma, mistag, eventNumber)
    count = 0
    sigma.setValue(0.1)
    mistag.setValue(0)

    m12val, m13val, dtval = [], [], []

    for i in range(m12.getNumBins()):
        m12.setValue(m12.getLowerLimit() + m12.getBinSize() * (i + 0.5))
        for j in range(m13.getNumBins()):
            m13.setValue(m13.getLowerLimit() + m13.getBinSize() * (j + 0.5))
            if inDalitz(m12.getValue(), m13.getValue()):
                for k in range(dtime.getNumBins()):
                    dtime.setValue(dtime.getLowerLimit() + dtime.getBinSize() * (k + 0.5))
                    m12val.append(m12.getValue())
                    m13val.append(m13.getValue())
                    dtval.append(dtime.getValue())
                    eventNumber.setValue(count)
                    data_tmp.addEvent()
                    count += 1

    data_tmp.setValueForAllEvents(sigma)
    data_tmp.setValueForAllEvents(mistag)
    print(data_tmp.getNumEvents())

    signal.setDataSize(data_tmp.getNumEvents(), 6)

    prod = ProdPdf('prod', [signal])

    prod.setData(data_tmp)
    prod.normalize()
    vals = prod.getCompProbsAtDataPoints()

    cumsum = np.cumsum(vals[0])
    cumsum = cumsum/cumsum[-1]

    sigma.setValue(0.1)
    mistag.setValue(0)
    for i in range(nEvents):
        r = np.random.rand()
        index = np.argmax(cumsum > r) # index - 1 < r < index 
    
        rand1 = np.random.rand() - 0.5
    
        currm12 = data_tmp.getValue(m12, index) + m12.getBinSize()*rand1
        currm13 = data_tmp.getValue(m13, index) + m13.getBinSize()*rand1
        currdtime = data_tmp.getValue(dtime, index) + dtime.getBinSize()*rand1

        m12.setValue(currm12)
        m13.setValue(currm13)
        dtime.setValue(currdtime)
        eventNumber.setValue(i)
        dataset.addEvent()

    dataset.setValueForAllEvents(sigma)
    dataset.setValueForAllEvents(mistag)

def main():
    parser = argparse.ArgumentParser(description='Toy fits.')
    parser.add_argument('first_toy', type=strtobool, help='Run one or 100 toy fits.', choices=[True, False])
    args = parser.parse_args()
    first_toy_flag = args.first_toy

    m12 = Observable('m12', 0, 3.5)
    m13 = Observable('m13', 0, 3.5)
    dtime = Observable('dtime', 0, 6) # D0 decay time in pico seconds
    sigma = Observable('sigma', 0.00001, 0.101)
    mistag = Observable('mistag', 0, 1) # Try mistag


    eventNumber = EventNumber('eventNumber')
    res = TruthResolution()

    signal = makeSignalPdf(m12, m13, dtime, sigma, mistag, eventNumber, res)

    print(signal.getObservables())
    print(dtime.getNumBins())

    if first_toy_flag:
        dataset = UnbinnedDataSet(m12, m13, dtime, sigma, mistag, eventNumber)

        fillToyDataSet(1000000, dataset, signal, m12, m13, dtime, sigma, mistag, eventNumber)

        signal.setData(dataset)
        signal.setDataSize(dataset.getNumEvents(), 6)

        fitman = FitManager(signal)

        print('Running fit...')
        func_min = fitman.fit()

        makePlots(signal, dataset, dtime, sigma, mistag, m12, m13, eventNumber)
    
    else:
        if os.path.exists('pulls.root'): os.remove('pulls.root')
        # Save initial parameters before running fit 100 times.
        params = {}
        for param in signal.getParameters():
            if not param.IsFixed():
                print('Parameter: {} value: {}'.format(param.getName(), str(param.getValue())))
                params[param.getName()] = param.getValue()
        dictionary = {}
        for i in range(100):
            for param in signal.getParameters():
                if not param.IsFixed():
                    param.setValue(params[param.getName()])
            dataset = UnbinnedDataSet(m12, m13, dtime, sigma, mistag, eventNumber)

            fillToyDataSet(1000000, dataset, signal, m12, m13, dtime, sigma, mistag, eventNumber)

            signal.setData(dataset)
            signal.setDataSize(dataset.getNumEvents(), 6)

            fitman = FitManager(signal)

            print('Running fit...')
            func_min = fitman.fit()
            if ( (not func_min.HasPosDefCovar()) or func_min.HesseFailed()): continue
            for fit_param in signal.getParameters():
                if not fit_param.IsFixed():
                    print('After fit Parameter: {} value: {}'.format(str(fit_param.getName()), str(fit_param.getValue())))
                    dictionary[fit_param.getName()] = (params[fit_param.getName()] - fit_param.getValue())/fit_param.getError()

            df = pd.DataFrame(dictionary, index=[i])
            df.to_root('pulls.root', 'mytree', mode='a')

if __name__ == '__main__':
    sys.exit(main())
