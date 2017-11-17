#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

from goofit import *
import numpy as np

print_goofit_info()

def main():
    mode = 0, data = 0
    plot

    dm = Variable("dm", 0.1395, 0.1665)
    dm.setNumBins(2700)

    mc_dataset = ()
    data_dataset = ()

    if mode == 0:
        mc_dataset.reset(UnbinnedDataSet(dm))
        data_dataset.reset(UnbinnedDataSet(dm))
    else:
        mc_dataset.reset(BinnedDataSet(dm))
        data_dataset.reset(BinnedDataSet(dm))

    mean1 = Variable("kpi_mc_mean1", 0.145402, 0.00001, 0.143, 0.148)
    mean2 = Variable("kpi_mc_mean2", 0.145465, 0.00001, 0.145, 0.1465)
    mean3 = Variable("kpi_mc_mean3", 0.145404, 0.00001, 0.144, 0.147)

    sigma1 = Variable("kpi_mc_sigma1", 0.00010, 0.00001, 0.000001, 0.002)
    sigma2 = Variable("kpi_mc_sigma2", 0.00075, 0.00001, 0.000001, 0.005)
    sigma3 = Variable("kpi_mc_sigma3", 0.00020, 0.00001, 0.000005, 0.001)

    pimass = Variable("kpi_mc_pimass", 0.13957)
    aslope = Variable("kpi_mc_aslope", -20.0, 1, -100.0, 10.0)
    apower = Variable("kpi_mc_apower", 1.3, 0.1, 0.1, 10.0)
    gfrac1 = Variable("kpi_mc_gfrac1", 0.65, 0.01, 0.0, 0.9)
    gfrac2 = Variable("kpi_mc_gfrac2", 0.02, 0.001, 0.0, 0.12)
    afrac = Variable("kpi_mc_afrac", 0.005, 0.003, 0.0, 0.1)

    gauss1 = GaussianPdf("gauss1", dm, mean1, sigma1)
    gauss2 = GaussianPdf("gauss2", dm, mean2, sigma2)
    gauss3 = GaussianPdf("gauss3", dm, mean3, sigma3)
    argus = ArgusPdf("argus", dm, pimass, aslope, false, apower)

    resolution = AddPdf("resolution", (gfrac1, gfrac2, afrac), (gauss1, gauss2, argus, gauss3))

    resolution.setData(mc_dataset.get())

    mcpdf = FitManagr(resolution)
    mcpdf.fit()

    # Locking the MC variables
    mean1.setFixed(True)
    mean2.setFixed(True)
    mean3.setFixed(True)
    sigma1.setFixed(True)
    sigma2.setFixed(True)
    sigma3.setFixed(True)
    pimass.setFixed(True)
    aslope.setFixed(True)
    gfrac1.setFixed(True)
    gfrac2.setFixed(True)
    afrac.setFixed(True)
    apower.setFixed(True)

    dummyzero = Variable("kpi_rd_dummyzero", 0)
    delta = Variable("kpi_rd_delta", 0.000002, -0.00005, 0.00005)
    epsilon = Variable("kpi_rd_epsilon", 0.05, -0.1, 0.2)

    resolution1 = ScaledGaussianPdf("resolution1", dm, dummyzero, sigma1, delta, epsilon)
    resolution2 = ScaledGaussianPdf("resolution2", dm, dummyzero, sigma2, delta, epsilon)
    resolution3 = ScaledGaussianPdf("resolution3", dm, dummyzero, sigma3, delta, epsilon)

    width_bw = Variable("kpi_rd_width_bw", 0.0001, 0.00001, 0.0005)
    rbw1 = KinLimitBWPdf("rbw1", dm, mean1, width_bw)
    rbw2 = KinLimitBWPdf("rbw2", dm, mean2, width_bw)
    rbw3 = KinLimitBWPdf("rbw3", dm, mean3, width_bw)

    signal1 = ConvolutionPdf("signal1", dm, rbw1, resolution1)
    signal2 = ConvolutionPdf("signal2", dm, rbw2, resolution2)
    signal3 = ConvolutionPdf("signal3", dm, rbw3, resolution3)

    signal1.setIntegrationConstants(0.1395, 0.1665, 0.0000027)
    signal2.setIntegrationConstants(0.1395, 0.1665, 0.0000027)
    signal3.setIntegrationConstants(0.1395, 0.1665, 0.0000027)

    signal = AddPdf("signal", (gfrac1, gfrac2, afrac), (signal1, signal2, argus, signal3))

    slope = Variable("kpi_rd_slope", -1.0, 0.1, -35.0, 25.0)
    bpower = nullptr
    bkg = ArgusPdf("bkg", dm, pimass, slope, false, bpower)
    bkg_frac = Variable("kpi_rd_bkg_frac", 0.03, 0.0, 0.3)
    total = AddPdf("total", (bkg_frac), (bkg, signal))
    total.setData(data_dataset.get())

    chi_control = ()
    if(2 == mode) (
        chi_control.reset(BinnedChisqFit)
        total.setFitControl(chi_control.get())
    )

    FitManager datapdf(total)
    datapdf.fit()
    return datapdf
)
main()
