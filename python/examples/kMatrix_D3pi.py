# -*- coding: utf-8 -*-
import sys
from math import cos, pi, sin, sqrt

from goofit import *

import matplotlib.pyplot as plt

print_goofit_info()

def Real(mag,phs): 
    return mag*cos(phs*pi/180.)

def Imag(mag,phs):
    return mag*sin(phs*pi/180.)

def makeSignalPdf(s12, s13, eventNumber):
    # https://arxiv.org/pdf/hep-ex/0312040.pdf
    D3pi = DecayInfo3()
    D3pi.motherMass = 1.86484
    D3pi.daug1Mass = 0.139570
    D3pi.daug2Mass = 0.139570
    D3pi.daug3Mass = 0.139570

    # beta coeeficients
    beta_r = [#ok
        Variable("beta1_r", 1),
        Variable("beta2_r", Real(2.471,82.5), 0.01, 0, 0),
        Variable("beta3_r", Real(1.084,102.8)),
        Variable("beta4_r", 0),
        Variable("beta5_r", 0),
    ]

    beta_i = [#ok
        Variable("beta1_i", 0),
        Variable("beta2_i", Imag(2.471,82.5), 0.01, 0, 0),
        Variable("beta3_i", Imag(1.084,102.8), 0.01, 0, 0),
        Variable("beta4_i", 0),
        Variable("beta5_i", 0),
    ]

    # f_prod coefficients
    
    f_prod_r = [#ok
        Variable("f_prod1_r", Real(2.565,155.4), 0.01, 0, 0),
        Variable("f_prod2_r", Real(6.312,-160.0), 0.01, 0, 0),
        Variable("f_prod3_r", 0),
        Variable("f_prod4_r", 0),
        Variable("f_prod5_r", 0),
    ]

    f_prod_i = [#ok
        Variable("f_prod1_i", Imag(2.565,155.4), 0.01, 0, 0),
        Variable("f_prod2_i", Imag(6.312,-160.0), 0.01, 0, 0),
        Variable("f_prod3_i", 0),
        Variable("f_prod4_i", 0),
        Variable("f_prod5_i", 0),
    ]

    

    # kMatrix
    fscat = [#ok
        Variable("f11_scatt", 0.26681),
        Variable("f12_scatt", 0.16583),
        Variable("f13_scatt", -0.19840),
        Variable("f14_scatt", 0.32808),
        Variable("f15_scatt", 0.31193),
    ]
    
    poles = [ #ok
        Variable("g00", 0.24844),
        Variable("g01", -0.52523),
        Variable("g02", 0),
        Variable("g03", -0.38878),
        Variable("g04", -0.36397),
        Variable("m0", 0.65100),
        Variable("g10", 0.91779),
        Variable("g11", 0.55427),
        Variable("g12", 0),
        Variable("g13", 0.38705),
        Variable("g14", 0.29448),
        Variable("m1", 1.20720),
        Variable("g20", 0.37024),
        Variable("g21", 0.23591),
        Variable("g22", 0.62605),
        Variable("g23", 0.18409),
        Variable("g24", 0.18923),
        Variable("m2", 1.56122),
        Variable("g30", 0.34501),
        Variable("g31", 0.39642),
        Variable("g32", 0.97644),
        Variable("g33", 0.19746),
        Variable("g34", 0.00357),
        Variable("m3", 1.21257),
        Variable("g40", 0.15770),
        Variable("g41", -0.17915),
        Variable("g42", -0.90100),
        Variable("g43", -0.00931),
        Variable("g44", 0.20689),
        Variable("m4", 1.81746),
    ]
   
    print("Define K Matrix Resonance.")
    print("beta_r =", beta_r)
    print("beta_i =", beta_i)
    print("f_prod_r =", f_prod_r)
    print("f_prod_i =", f_prod_i)

    a_r = Variable("a_r", sqrt(56) )
    a_i = Variable("a_i", 0)
    sA0 = Variable("sA0", 0)
    sA  = Variable("sA", 0)
    s0_prod = Variable("s0_prod", -1.0)
    s0_scatt = Variable("s0_scatt", -3.30564)

    kMatrix_12 = Resonances.kMatrix(
        "kMatrix_12",
        a_r,
        a_i,
        sA0, #ok
        sA,#ok
        s0_prod,#ok
        s0_scatt,#ok
        beta_r,
        beta_i,
        f_prod_r,
        f_prod_i,
        fscat,
        poles,
        0,  # L (spin)
        PAIR_12,
    )

    kMatrix_13 = Resonances.kMatrix(
        "kMatrix_13",
        a_r,
        a_i,
        sA0, #ok
        sA,#ok
        s0_prod,#ok
        s0_scatt,#ok
        beta_r,
        beta_i,
        f_prod_r,
        f_prod_i,
        fscat,
        poles,
        0,  # L (spin)
        PAIR_13,
    )


    rho770_re = Variable("rho770_re",sqrt(56./30.82)*Real(1.858,-139.4),0.01,0,0)
    rho770_im = Variable("rho770_im",sqrt(56./30.82)*Imag(1.858,-139.4),0.01,0,0)
    rho770_mass = Variable("rho770_mass",0.77526)
    rho770_width = Variable("rho770_width",0.1491)

    rho770_12 = Resonances.RBW(
        "rho770_12",
        rho770_re,
        rho770_im,
        rho770_mass,
        rho770_width,
        1,
        PAIR_12,
        True,
        False
        )
    
    rho770_13 = Resonances.RBW(
        "rho770_13",
        rho770_re,
        rho770_im,
        rho770_mass,
        rho770_width,
        1,
        PAIR_13,
        True,
        False
        )

    f2_re = Variable("f2_re",sqrt(56./11.74)*Real(1.147,-47.5),0.01,0,0)
    f2_im = Variable("f2_im",sqrt(56./11.74)*Imag(1.147,-47.5),0.01,0,0)
    f2_mass = Variable("f2_mass",1.2755)
    f2_width = Variable("f2_width",0.1867)

    f2_12 = Resonances.RBW(
        "f2_12",
        f2_re,
        f2_im,
        f2_mass,
        f2_width,
        2,
        PAIR_12,
        True,
        False
        )

    f2_13 = Resonances.RBW(
        "f2_13",
        f2_re,
        f2_im,
        f2_mass,
        f2_width,
        2,
        PAIR_13,
        True,
        False
        )   

    print("Add to list of resonances.")
    D3pi.resonances = [rho770_12,rho770_13,f2_12,f2_13,kMatrix_12,kMatrix_13]
    #D3pi.resonances = [kMatrix_12,kMatrix_13]
    #D3pi.resonances = [rho770_12,rho770_13]
    #D3pi.resonances = [f2_12,f2_13]

    print(kMatrix_12.get_amp_real())
    print(kMatrix_12.get_amp_img())
    print(kMatrix_12.getParameters())
    print(kMatrix_12.getObservables())

    print("Define efficiency and Amp3Body Pdf.")
    # Constant efficiency
    constantOne = Variable("One", 1)
    constantZero = Variable("Zero", 0)
    eff = PolynomialPdf(
        "constantEff", [s12, s13], [constantOne], [constantZero, constantZero], 0
    )
    d = Amp3Body("signalPDF", s12, s13, eventNumber, D3pi, eff)

    print(d.getName())
    print(d.getParameters())
    print(d.getObservables())

    return d


def make_toy_data(dp, s12, s13, eventNumber, nEvents):
    data = UnbinnedDataSet(s12, s13, eventNumber)

    print("ProdPdf")
    prod = ProdPdf("prod", [dp])

    print(s12.getNumBins())
    print(s13.getNumBins())

    dplotter = DalitzPlotter(prod,dp)

    plt.figure(0)
    arr = dplotter.make2D()
    extent = dplotter.getExtent()
    plt.imshow(arr, extent=extent, origin="lower")

    plt.savefig("dalitzD3pi_DP_plot.png")
    
    dplotter.fillDataSetMC(data,nEvents)

    data_matrix = data.to_matrix()

    plt.figure(1)
    plt.hist(data_matrix[0],bins=120)
    plt.savefig("dalitzD3pi_s12_plot.png")

    plt.figure(2)
    plt.hist(data_matrix[1],bins=120)
    plt.savefig("dalitzD3pi_s13_plot.png")

    return data

    
def runToyFit(dp, data):
    
    prod = ProdPdf("prod",[dp])
    prod.setData(data)
    dp.setDataSize(data.getNumEvents())
    fitman = FitManager(prod)

    fitman.fit()
    
def main():
    nBins = 100
    piMass = 0.139570
    DpMass = 1.86966
    s12 = Observable("s12", 4*piMass**2, (DpMass-piMass)**2)
    s13 = Observable("s13",  4*piMass**2, (DpMass-piMass)**2)
    s12.setNumBins(nBins)
    s13.setNumBins(nBins)
    eventNumber = EventNumber("eventNumber")

    signal = makeSignalPdf(s12, s13, eventNumber)
    data = make_toy_data(signal, s12, s13, eventNumber,100000)

    runToyFit(signal, data)

if __name__ == "__main__":
    sys.exit(main())