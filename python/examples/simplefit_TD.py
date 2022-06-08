import sys
import argparse
from goofit import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
from math import cos, sin, sqrt, pi, atan2
from root_pandas import read_root
from scipy.optimize import curve_fit
from scipy.stats import norm


mode = 'SingleTag_D0ToKsPiPiDD'
year = '2016'

# Observables
m12 = Observable('m12', 0, 3.5)
m13 = Observable('m13', 0, 3.5)
dtime = Observable('dtime', -1, 6) # D0 decay time in pico seconds
sigma = Observable('sigma', 0.00001, 0.101)
mistag = Observable('mistag', 0, 1) # Try mistag
charmtag = Observable('charmtag', -1, 1)
sigprob = Observable('sigprob', 0, 1)
eventNumber = EventNumber('eventNumber')
category = Observable('category', 0, 20)

tau = Variable("tau", 0.4101, 0.001, 0, 0)
xcp = Variable("xcp", 0.0051, 0.001, 0, 0)
ycp = Variable("ycp", 0.0063, 0.001, 0, 0)
deltax = Variable("deltax", 0)
deltay = Variable("deltay", 0)

fixAmps = True
fixMass = True
fixKmatrix = True

rho_770_amp_real = Variable("rho_770_amp_real", 1)
rho_770_amp_imag = Variable("rho_770_amp_imag", 0)
rho_770_mass = Variable("rho_770_mass", 0.77526, 0.0001, 0, 0)
rho_770_width = Variable("rho_770_width", 0.1478, 0.0001, 0.1300, 0.1500)
f2_1270_amp_real = Variable("f2_1270_amp_real", 1.43*cos(-36.3*pi/180)) if fixAmps else Variable("f2_1270_amp_real", 1.43*cos(-36.3*pi/180), 0.0001, 0, 0)
f2_1270_amp_imag = Variable("f2_1270_amp_imag", 1.43*sin(-36.3*pi/180)) if fixAmps else Variable("f2_1270_amp_imag", 1.43*sin(-36.3*pi/180), 0.0001, 0, 0)
#f2_1270_mass = Variable("f2_1270_mass", 1.2755, 0.0001, 1.2755-20*0.0008, 1.2755+20*0.0008)
f2_1270_width = Variable("f2_1270_width", 0.1867, 0.0001, 0.1867-3*0.0025, 0.1867+3*0.0022)
f2_1270_mass = Variable("f2_1270_mass", 1.2755, 0.0001, 1.250, 1.300)
omega_782_amp_real = Variable("omega_782_amp_real", 0.0388*cos(120.7*pi/180)) if fixAmps else Variable("omega_782_amp_real", 0.0388*cos(120.7*pi/180), 0.0001, 0, 0)
omega_782_amp_imag = Variable("omega_782_amp_imag", 0.0388*sin(120.7*pi/180)) if fixAmps else Variable("omega_782_amp_imag", 0.0388*sin(120.7*pi/180), 0.0001, 0, 0)
#omega_782_width = Variable("omega_782_width", 0.00849, 0.0001, 0.00849-20*0.00008, 0.00849+20*0.00008)
omega_782_mass = Variable("omega_782_mass", 0.78265, 0.0001, 0, 0)
omega_782_width = Variable("omega_782_width", 0.00849, 0.0001, 0.008, 0.1)
rho_1450_amp_real = Variable("rho_1450_amp_real", 2.85*cos(102.1*pi/180)) if fixAmps else Variable("rho_1450_amp_real", 2.85*cos(102.1*pi/180), 0.0001, 0, 0)
rho_1450_amp_imag = Variable("rho_1450_amp_imag", 2.85*sin(102.1*pi/180)) if fixAmps else Variable("rho_1450_amp_imag", 2.85*sin(102.1*pi/180), 0.0001, 0, 0)
#rho_1450_mass = Variable("rho_1450_mass", 1.465, 0.0001, 1.465-20*0.00025, 1.465+20*0.00025)
rho_1450_mass = Variable("rho_1450_mass", 1.465, 0.0001, 1.450, 1.510)
#rho_1450_width = Variable("rho_1450_width", 0.400, 0.0001, 0.400-20*0.060, 0.400+20*0.060)
rho_1450_width = Variable("rho_1450_width", 0.400, 0.0001, 0, 0)

Kstarm_892_amp_real = Variable("Kstarm_892_amp_real", 1.720*cos(136.8*pi/180)) if fixAmps else Variable("Kstarm_892_amp_real", 1.720*cos(136.8*pi/180), 0.0001, 0, 0)
Kstarm_892_amp_imag = Variable("Kstarm_892_amp_imag", 1.720*sin(136.8*pi/180)) if fixAmps else Variable("Kstarm_892_amp_imag", 1.720*sin(136.8*pi/180), 0.0001, 0, 0)
Kstarm_892_mass = Variable("Kstarm_892_mass", 0.89166, 0.0001, 0, 0)
Kstarm_892_width = Variable("Kstarm_892_width", 0.0508, 0.0001, 0.0400, 0.0550)
K2starm_1430_amp_real = Variable("K2starm_1430_amp_real", 1.27*cos(-44.1*pi/180)) if fixAmps else Variable("K2starm_1430_amp_real", 1.27*cos(-44.1*pi/180), 0.0001, 0, 0)
K2starm_1430_amp_imag = Variable("K2starm_1430_amp_imag", 1.27*sin(-44.1*pi/180)) if fixAmps else Variable("K2starm_1430_amp_imag", 1.27*sin(-44.1*pi/180), 0.0001, 0, 0)
#K2starm_1430_mass = Variable("K2starm_1430_mass", 1.4273, 0.0001, 1.4273-3*0.0015, 1.4273+3*0.0015)
K2starm_1430_mass = Variable("K2starm_1430_mass", 1.4273, 0.0001, 1.420, 1.440)
#K2starm_1430_width = Variable("K2starm_1430_width", 0.1000, 0.0001, 0.1000-5*0.0022, 0.1000+5*0.0022)
K2starm_1430_width = Variable("K2starm_1430_width", 0.1000, 0.0001, 0, 0 )
Kstarm_1680_amp_real = Variable("Kstarm_1680_amp_real", 3.31*cos(-118.2*pi/180)) if fixAmps else Variable("Kstarm_1680_amp_real", 3.31*cos(-118.2*pi/180), 0.0001, 0, 0)
Kstarm_1680_amp_imag = Variable("Kstarm_1680_amp_imag", 3.31*sin(-118.2*pi/180)) if fixAmps else Variable("Kstarm_1680_amp_imag", 3.31*sin(-118.2*pi/180), 0.0001, 0, 0)
Kstarm_1680_mass = Variable("Kstarm_1680_mass", 1.718, 0.0001, 1.718-3*0.018, 1.718+3*0.018)
#Kstarm_1680_width = Variable("Kstarm_1680_width", 0.322, 0.0001, 0, 0.322+3*0.110)
#Kstarm_1680_mass = Variable("Kstarm_1680_mass", 1.718, 0.0001, 0, 0)
Kstarm_1680_width = Variable("Kstarm_1680_width", 0.322, 0.0001, 0, 0)
Kstarm_1410_amp_real = Variable("Kstarm_1410_amp_real", 0.29*cos(99.4*pi/180)) if fixAmps else Variable("Kstarm_1410_amp_real", 0.29*cos(99.4*pi/180), 0.0001, 0, 0)
Kstarm_1410_amp_imag = Variable("Kstarm_1410_amp_imag", 0.29*sin(99.4*pi/180)) if fixAmps else Variable("Kstarm_1410_amp_imag", 0.29*sin(99.4*pi/180), 0.0001, 0, 0)
Kstarm_1410_mass = Variable("Kstarm_1410_mass", 1.414, 0.0001, 1.414-2*0.015, 1.414+2*0.015)
#Kstarm_1410_mass = Variable("Kstarm_1410_mass", 1.414, 0.0001, 1.410, 1.440)
Kstarm_1410_width = Variable("Kstarm_1410_width", 0.232, 0.0001, 0.232-3*0.021, 0.232+3*0.021)

Kstarp_892_amp_real = Variable("Kstarp_892_amp_real", 0.164*cos(-42.2*pi/180)) if fixAmps else Variable("Kstarp_892_amp_real", 0.164*cos(-42.2*pi/180), 0.0001, 0, 0)
Kstarp_892_amp_imag = Variable("Kstarp_892_amp_imag", 0.164*sin(-42.2*pi/180)) if fixAmps else Variable("Kstarp_892_amp_imag", 0.164*sin(-42.2*pi/180), 0.0001, 0, 0)
K2starp_1430_amp_real = Variable("K2starp_1430_amp_real", 0.10*cos(-89.6*pi/180)) if fixAmps else Variable("K2starp_1430_amp_real", 0.10*cos(-89.6*pi/180), 0.0001, 0, 0)
K2starp_1430_amp_imag = Variable("K2starp_1430_amp_imag", 0.10*sin(-89.6*pi/180)) if fixAmps else Variable("K2starp_1430_amp_imag", 0.10*sin(-89.6*pi/180), 0.0001, 0, 0)
Kstarp_1410_amp_real = Variable("Kstarp_1410_amp_real", 0.21*cos(150.2*pi/180)) if fixAmps else Variable("Kstarp_1410_amp_real", 0.21*cos(150.2*pi/180), 0.0001, 0, 0)
Kstarp_1410_amp_imag = Variable("Kstarp_1410_amp_imag", 0.21*sin(150.2*pi/180)) if fixAmps else Variable("Kstarp_1410_amp_imag", 0.21*sin(150.2*pi/180), 0.0001, 0, 0)

# LASS Kpi S wave parameters
K0star_1430_mass = Variable("K0star_1430_mass", 1.425, 0.0001, 1.425-3*0.050, 1.425+3*0.050)
K0star_1430_width = Variable("K0star_1430_width", 0.270, 0.0001, 0.270-3*0.080, 0.270+3*0.080)
#K0star_1430_mass = Variable("K0star_1430_mass", 1.425, 0.0001, 1.400, 1.450)
#K0star_1430_width = Variable("K0star_1430_width", 0.270, 0.0001, 0.200, 0.300)

fixLass = True
a = Variable("a", 0.113) if fixLass else Variable("a", 0.113, 0.0001, 0, 0)
r = Variable("r", -33.8) if fixLass else Variable("r", -33.8, 0.0001, 0, 0)
R = Variable("R", 1) # Fixed
phiR = Variable("phiR", -1.91462619) if fixLass else Variable("phiR", -1.91462619, 0.0001, 0, 0)
F = Variable("F", 0.96) if fixLass else Variable("F", 0.96, 0.0001, 0, 0)
phiF = Variable("phiF", 0.00174533) if fixLass else Variable("phiF", 0.00174533, 0.0001, 0, 0)

K0starp_1430_amp_real = Variable("K0starp_1430_amp_real", 0.11*cos(162.3*pi/180)) if fixAmps else Variable("K0starp_1430_amp_real", 0.11*cos(162.3*pi/180), 0.0001, 0, 0)
K0starp_1430_amp_imag = Variable("K0starp_1430_amp_imag", 0.11*sin(162.3*pi/180)) if fixAmps else Variable("K0starp_1430_amp_imag", 0.11*sin(162.3*pi/180), 0.0001, 0, 0)
K0starm_1430_amp_real = Variable("K0starm_1430_amp_real", 2.36*cos(99.4*pi/180)) if fixAmps else Variable("K0starm_1430_amp_real", 2.36*cos(99.4*pi/180), 0.0001, 0, 0)
K0starm_1430_amp_imag = Variable("K0starm_1430_amp_imag", 2.36*sin(99.4*pi/180)) if fixAmps else Variable("K0starm_1430_amp_imag", 2.36*sin(99.4*pi/180), 0.0001, 0, 0)

#kMatrix_amp_real = Variable("kMatrix_amp_real", 8.0*cos(-126.0*pi/180)) if fixAmps else Variable("kMatrix_amp_real", 8.0*cos(-126.0*pi/180), 0.0001, 0, 0)
#kMatrix_amp_imag = Variable("kMatrix_amp_imag", 8.0*sin(-126.0*pi/180)) if fixAmps else Variable("kMatrix_amp_imag", 8.0*sin(-126.0*pi/180), 0.0001, 0, 0)
kMatrix_amp_real = Variable("kMatrix_amp_real", 1) if fixAmps else Variable("kMatrix_amp_real", 1, 0.0001, 0, 0)
kMatrix_amp_imag = Variable("kMatrix_amp_imag", 0) if fixAmps else Variable("kMatrix_amp_imag", 0, 0.0001, 0, 0)

def makeDecayInfo(mesrad, motherrad):
    d0toks0pipi = DecayInfo3t(tau, xcp, ycp, deltax, deltay)

    d0toks0pipi.motherMass = 1.86484
    d0toks0pipi.daug1Mass = 0.497614
    d0toks0pipi.daug2Mass = 0.139570
    d0toks0pipi.daug3Mass = 0.139570
    d0toks0pipi.meson_radius = mesrad
    d0toks0pipi.mother_meson_radius = motherrad


    # Gounaris Sakurai (name, amplitude real, amplitude imaginary, mass, width, spin, cyclic index)
    rho_770 = Resonances.GS("rho_770", 
                            rho_770_amp_real,
                            rho_770_amp_imag,
                            rho_770_mass,
                            rho_770_width,
                            1,
                            PAIR_23)

    # Relativistic Breit Wigner (name, amplitude real, amplitude imaginary, mass, width, spin, cyclic index)
    f2_1270 = Resonances.RBW("f2_1270",
                             f2_1270_amp_real,
                             f2_1270_amp_imag,
                             f2_1270_mass,
                             f2_1270_width,
                             2,
                             PAIR_23)

    omega_782 = Resonances.RBW("omega_782",
                               omega_782_amp_real,
                               omega_782_amp_imag,
                               omega_782_mass,
                               omega_782_width,
                               1,
                               PAIR_23)

    rho_1450 = Resonances.RBW("rho_1450",
                              rho_1450_amp_real,
                              rho_1450_amp_imag,
                              rho_1450_mass,
                              rho_1450_width,
                              1,
                              PAIR_23)

    # KS0 pi- resonances 
    Kstarm_892 = Resonances.RBW("Kstarm_892",
                                Kstarm_892_amp_real,
                                Kstarm_892_amp_imag,
                                Kstarm_892_mass,
                                Kstarm_892_width,
                                1,
                                PAIR_13)

    K2starm_1430 = Resonances.RBW("K2starm_1430",
                                  K2starm_1430_amp_real,
                                  K2starm_1430_amp_imag,
                                  K2starm_1430_mass,
                                  K2starm_1430_width,
                                  2,
                                  PAIR_13)

    Kstarm_1680 = Resonances.RBW("Kstarm_1680",
                                 Kstarm_1680_amp_real,
                                 Kstarm_1680_amp_imag,
                                 Kstarm_1680_mass,
                                 Kstarm_1680_width,
                                 1,
                                 PAIR_13)

    Kstarm_1410 = Resonances.RBW("Kstarm_1410",
                                 Kstarm_1410_amp_real,
                                 Kstarm_1410_amp_imag,
                                 Kstarm_1410_mass,
                                 Kstarm_1410_width,
                                 1,
                                 PAIR_13)

    # KS0 pi+ resonances 
    Kstarp_892 = Resonances.RBW("Kstarp_892",
                                Kstarp_892_amp_real,
                                Kstarp_892_amp_imag,
                                Kstarm_892_mass,
                                Kstarm_892_width,
                                1,
                                PAIR_12)

    K2starp_1430 = Resonances.RBW("K2starp_1430",
                                  K2starp_1430_amp_real,
                                  K2starp_1430_amp_imag,
                                  K2starm_1430_mass,
                                  K2starm_1430_width,
                                  2,
                                  PAIR_12)

    Kstarp_1410 = Resonances.RBW("Kstarp_1410",
                                  Kstarp_1410_amp_real,
                                  Kstarp_1410_amp_imag,
                                  Kstarm_1410_mass,
                                  Kstarm_1410_width,
                                  1,
                                  PAIR_12)

    # LASS: Kpi S wave
    K0starp_1430 = Resonances.LASS("K0starp_1430",
                                   K0starp_1430_amp_real,
                                   K0starp_1430_amp_imag,
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
                                   K0starm_1430_amp_real,
                                   K0starm_1430_amp_imag,
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
    # beta coeeficients
    beta_r = [Variable("beta1_r", 8.5*cos(68.5*pi/180)  , 0.0001, 0, 0),
              Variable("beta2_r", 12.2*cos(24.0*pi/180) , 0.0001, 0, 0),
              Variable("beta3_r", 29.2*cos(-0.1*pi/180) , 0.0001, 0, 0),
              Variable("beta4_r", 10.8*cos(-51.9*pi/180), 0.0001, 0, 0),
              Variable("beta5_r", 0)
             ]

    beta_i = [Variable("beta1_i", 8.5*sin(68.5*pi/180)  , 0.0001, 0, 0),
              Variable("beta2_i", 12.2*sin(24.0*pi/180) , 0.0001, 0, 0),
              Variable("beta3_i", 29.2*sin(-0.1*pi/180) , 0.0001, 0, 0),
              Variable("beta4_i", 10.8*sin(-51.9*pi/180), 0.0001, 0, 0),
              Variable("beta5_i", 0)
             ]
    
    # f_prod coefficients
    f_prod_r = [Variable("f_prod1_r", 8.0*cos(-126.0*pi/180) ),
                Variable("f_prod2_r", 26.3*cos(-152.3*pi/180), 0.0001, 0, 0),
                Variable("f_prod3_r", 33.0*cos(-93.2*pi/180) , 0.0001, 0, 0),
                Variable("f_prod4_r", 26.2*cos(-121.4*pi/180), 0.0001, 0, 0),
                Variable("f_prod5_r", 0)
               ]

    f_prod_i = [Variable("f_prod1_i", 8.0*sin(-126.0*pi/180) ),
                Variable("f_prod2_i", 26.3*sin(-152.3*pi/180), 0.0001, 0, 0),
                Variable("f_prod3_i", 33.0*sin(-93.2*pi/180) , 0.0001, 0, 0),
                Variable("f_prod4_i", 26.2*sin(-121.4*pi/180), 0.0001, 0, 0), 
                Variable("f_prod5_i", 0)                     
               ]

    # kMatrix
    fscat = [Variable("f11_scatt", 0.23399),
             Variable("f12_scatt", 0.15044),
             Variable("f13_scatt", -0.20545),
             Variable("f14_scatt", 0.32825),
             Variable("f15_scatt", 0.35412)]
    poles = [Variable("g00", 0.22889),
             Variable("g01", -0.55377),
             Variable("g02", 0),
             Variable("g03", -0.39899),
             Variable("g04", -0.34639),
             Variable("m0", 0.65100),
             Variable("g10", 0.94128),
             Variable("g11", 0.55095),
             Variable("g12", 0),
             Variable("g13",  0.39065),
             Variable("g14", 0.31503),
             Variable("m1", 1.20360),
             Variable("g20", 0.36856),
             Variable("g21", 0.23888),
             Variable("g22", 0.55639),
             Variable("g23", 0.18340),
             Variable("g24", 0.18681),
             Variable("m2", 1.55817),
             Variable("g30", 0.33650),
             Variable("g31", 0.40907),
             Variable("g32", 0.85679),
             Variable("g33",  0.19906),
             Variable("g34", -0.00984),
             Variable("m3", 1.21000),
             Variable("g40", 0.18171),
             Variable("g41", -0.17558),
             Variable("g42", -0.79658),
             Variable("g43", -0.00355),
             Variable("g44", 0.22358),
             Variable("m4", 1.82206)]
    print('Define K Matrix Resonance.')
    print("beta_r =", beta_r)
    print("beta_i =", beta_i)
    print("f_prod_r =", f_prod_r)
    print("f_prod_i =", f_prod_i)

    kMatrix = Resonances.kMatrix("kMatrix",
                                 kMatrix_amp_real,
                                 kMatrix_amp_imag,
                                 Variable("sA0", -0.15),
                                 Variable("sA", 1),
                                 Variable("s0_prod", -0.07),
                                 Variable("s0_scatt", -3.92637),
                                 beta_r,
                                 beta_i,
                                 f_prod_r,
                                 f_prod_i,
                                 fscat,
                                 poles,
                                 0, # L (spin)
                                 PAIR_23
    )

    d0toks0pipi.resonances = ([rho_770, omega_782, f2_1270, rho_1450, 
                              Kstarm_892, K2starm_1430, Kstarm_1410, Kstarm_1680,
                              Kstarp_892, K2starp_1430, Kstarp_1410,
                              K0starm_1430, K0starp_1430,
                              kMatrix]
                              )

    return d0toks0pipi



def getDatasets(years, modes):
  datasets = {}
  filename = 'root://eoslhcb.cern.ch//eos/lhcb/user/m/mhilton/KSPiPi-ntuples/tuples-BDT/sWeight_{mode}_{year}.root'
  catIndex = 0
  for year in years:
    for mode in modes:
      print(f'using category {catIndex} for {mode} {year}')
      dataset = getData(filename.format(mode=mode, year=year), 'tree_sWeights', catIndex)
      datasets[f'{year}_{mode}'] = dataset
      catIndex = catIndex + 1
  for key, value in datasets.items():
    print(f'Number of events in {key} dataset: ', value.getNumEvents())
  return datasets


def getData(filename, tree, catIndex):
    df = read_root(filename, tree, columns=['m12', 'm13', 'sig', 'B_DTFKS_D0_CTAU',
                                            'D0_ID'])
    df = df.query('B_DTFKS_D0_CTAU/0.3 < 6')

    df["category"] = catIndex
    df.reset_index(drop=True, inplace=True)
    df['index'] = range(df.shape[0])

    df['tau']    = df.eval('B_DTFKS_D0_CTAU/0.3')
    df['index']  = np.arange(df.shape[0])
    df['sigma']  = np.repeat(0.1, df.shape[0])
    df['mistag'] = np.repeat(0.005, df.shape[0])
    df['D0_tag'] = np.where(df.D0_ID == 421, np.ones(df.shape[0]), np.repeat(-1, df.shape[0]))

    print(df)
    array = df[['m12', 'm13', 'tau', 'sig', 'sigma', 'mistag', 'D0_tag', 'category', 'index']].to_numpy()

    dataset = UnbinnedDataSet(m12, m13, dtime, sigprob, sigma, mistag, charmtag, category, eventNumber)

    dataset.from_matrix(array.T)

    return dataset

def makeEfficiency(mode, year):
    # SquareDalitzEfficiency
    path = '../PSEfficiency/'
    eff_params = pd.read_csv(f'{path}Efficiency_params_{mode}_{year}.csv')
    c = eff_params.params

    eff = SquareDalitzEffPdf('eff', m12, m13, 
                             Variable(f'c0_{mode}_{year}', c[0]),
                             Variable(f'c1_{mode}_{year}', c[1]),
                             Variable(f'c2_{mode}_{year}', c[2]),
                             Variable(f'c3_{mode}_{year}', c[3]),
                             Variable(f'c4_{mode}_{year}', c[4]),
                             Variable(f'c5_{mode}_{year}', c[5]),
                             Variable(f'c6_{mode}_{year}', c[6]),
                             Variable(f'c7_{mode}_{year}', c[7]) 
                             )
    eff.setParameterConstantness(True)

    return eff

def fillBGhist(filename, bghist_dalitz, bghist_dtime):
    # Background Pdf. Get data use one file for now.
    columns =['sWeight_bkg', 'm12', 'm13', 'B_DTFKS_D0_CTAU', 'B_DTFKS_D0_M']
    dataframe = read_root(filename, 'tree_sWeights')
    hist, _, _ = np.histogram2d(dataframe.m12, dataframe.m13, weights=dataframe.sWeight_bkg, 
                                density=True, bins=m12.getNumBins(), 
                                range=[(m12.getLowerLimit(),m12.getUpperLimit()),(m12.getLowerLimit(),m12.getUpperLimit())])

    histnan = np.copy(hist) #copy
    histnan[histnan == 0] = np.nan
    hist = np.swapaxes(hist,0,1)
    hist = hist.flatten()
    for i in range(bghist_dalitz.getNumBins()):
        if inDalitz(bghist_dalitz.getBinCenter(0,i), bghist_dalitz.getBinCenter(1,i)):
            bghist_dalitz.setBinContent(i,hist[i])
        elif hist[i] < 0:
            bghist_dalitz.setBinContent(i,0)
        else:
            bghist_dalitz.setBinContent(i,0)


    hdtime, _ = np.histogram(dataframe.B_DTFKS_D0_CTAU/0.3, weights=dataframe.sWeight_bkg, 
                             bins=dtime.getNumBins(), range=(dtime.getLowerLimit(),dtime.getUpperLimit()))

    # Decay time background pdf
    for i in range(bghist_dtime.getNumBins()):
        bghist_dtime.setBinContent(i, hdtime[i])


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

def GetD0Fraction(dataset, charmtag, sigprob):
    nD0 = 0
    sumsigprob = 0
    for i in range(dataset.getNumEvents()):
        sumsigprob += dataset.getValue(sigprob, i)
        if dataset.getValue(charmtag, i) == 1: nD0 += dataset.getValue(sigprob, i)

    D0Fraction = nD0/sumsigprob
    return D0Fraction

def getSigfrac(dataset, sigprob):
    sumsigprob = 0
    for i in range(dataset.getNumEvents()):
        sumsigprob += dataset.getValue(sigprob, i)
    sigfrac = sumsigprob / dataset.getNumEvents()

    print('Total sigprob: ', sumsigprob)

    return sigfrac, sumsigprob

def Getm23(m12, m13):
    mD0, mKS, mPi = 1.86486, 0.497614, 0.13957
    m23 = mD0**2 + mKS**2 + 2*mPi**2 - m12 - m13
    return m23


def makePlots(total, signal, dataset, bkg_flag, string):

    data_tmp   = UnbinnedDataSet(m12, m13, dtime, sigprob, sigma, mistag, charmtag, eventNumber)
    data_dtime = UnbinnedDataSet(m12, m13, dtime, sigprob, sigma, mistag, charmtag, eventNumber)

    d0fraction = GetD0Fraction(dataset, charmtag, sigprob)
    sigfrac, totalsigprob = getSigfrac(dataset, sigprob)
    totalbkprob = dataset.getNumEvents() - totalsigprob
    count = 0
    sigma.setValue(0.1)
    mistag.setValue(0.005)
    charmtag.setValue(1)
    if bkg_flag: sigprob.setValue(sigfrac)
    else: sigprob.setValue(1)

    m12val, m13val, m23val, dtval = [], [], [], []

    for i in range(875):
        m12.setValue(m12.getLowerLimit() + (m12.getUpperLimit() - m12.getLowerLimit()) * (i + 0.5) / 875)
        for j in range(875):
            m13.setValue(m13.getLowerLimit() + (m13.getUpperLimit() - m13.getLowerLimit()) * (j + 0.5) / 875)
            if inDalitz(m12.getValue(), m13.getValue()):
                for k in range(10):
                    dtime.setValue(dtime.getLowerLimit() + (dtime.getUpperLimit() - dtime.getLowerLimit()) * (k + 0.5) / 10)
                    m12val.append(m12.getValue())
                    m13val.append(m13.getValue())
                    m23val.append(Getm23(m12.getValue(), m13.getValue()))
                    if np.random.uniform() < d0fraction: charmtag.setValue(1)
                    else: charmtag.setValue(-1)
                    eventNumber.setValue(count)
                    data_tmp.addEvent()
                    count += 1

    data_tmp.setValueForAllEvents(sigma)
    data_tmp.setValueForAllEvents(mistag)
    data_tmp.setValueForAllEvents(charmtag)
    data_tmp.setValueForAllEvents(sigprob)

    signal.setDataSize(data_tmp.getNumEvents(), 8)

    total.setData(data_tmp)
    total.normalize()
    vals = total.getCompProbsAtDataPoints()
    array = np.array(vals)
    totalPdf = np.sum(array[0])
    totalPdf_sig = np.sum(array[1])
    if bkg_flag: totalPdf_bkg = np.sum(array[2])
    #totalPdf = totalPdf_sig + totalPdf_bkg

    m12dat = np.array(dataset.to_numpy())[0]
    m13dat = np.array(dataset.to_numpy())[1]
    dtdat  = np.array(dataset.to_numpy())[2]
    m23dat = Getm23(m12dat, m13dat)

    _, ax = plt.subplots(1, 2, figsize=(10, 5))
    hist_dat, binsx, binsy, _ = ax[0].hist2d(m12dat, m13dat, bins=175, range=[(m12.getLowerLimit(),m12.getUpperLimit()),(m13.getLowerLimit(),m13.getUpperLimit())], cmap='hot_r')
    bin_centrex = (binsx[:-1] + binsx[1:])/2
    bin_centrey = (binsy[:-1] + binsy[1:])/2
    if bkg_flag:
        hist_sig, _, _ = np.histogram2d(m12val, m13val, weights=array[1], bins=175, range=[(m12.getLowerLimit(),m12.getUpperLimit()),(m12.getLowerLimit(),m12.getUpperLimit())])
        hist_bkg, _, _ = np.histogram2d(m12val, m13val, weights=array[2], bins=175, range=[(m12.getLowerLimit(),m12.getUpperLimit()),(m12.getLowerLimit(),m12.getUpperLimit())])
        hist_sig = hist_sig*totalsigprob/totalPdf_sig
        hist_bkg = hist_bkg*totalbkprob/totalPdf_bkg
        hist_fit = np.add(hist_sig, hist_bkg)
        #hist_fit, _, _ = np.histogram2d(m12val, m13val, weights=array[0], bins=175, range=[(m12.getLowerLimit(),m12.getUpperLimit()),(m12.getLowerLimit(),m12.getUpperLimit())])
        #hist_fit = hist_fit*dataset.getNumEvents()/totalPdf
    else: 
        hist_fit, _, _ = np.histogram2d(m12val, m13val, weights=array[0], bins=175, range=[(m12.getLowerLimit(),m12.getUpperLimit()),(m12.getLowerLimit(),m12.getUpperLimit())])
        hist_fit = hist_fit*dataset.getNumEvents()/totalPdf
    ax[1].imshow(np.swapaxes(hist_fit,0,1), origin='lower', extent=(m12.getLowerLimit(),m12.getUpperLimit(),m12.getLowerLimit(),m12.getUpperLimit()), cmap='hot_r')

    for a in ax:
        a.set_xlabel('m12')
        a.set_ylabel('m13')
        a.axis('equal')
    ax[0].set_title('Data')
    ax[1].set_title('Fitted Function')
    plt.tight_layout()
    plt.savefig('Tddp_dalitz'+string+'.pdf')

    pulls = np.array((hist_dat - hist_fit)/np.sqrt(hist_dat))
    np.nan_to_num(pulls, copy=False, nan=0, posinf=0, neginf=0)
    nbins = np.count_nonzero(np.asarray(hist_dat > 20)) # number of non empty bins
    print('Number of non-empty bins: {}'.format(nbins))
    print('Number of bins contributing: {}'.format(np.count_nonzero(pulls)))
    chi2 = np.sum(np.square(pulls))
    print('Chi2: {}'.format(chi2))
    nfree = 0
    for param in signal.getParameters():
        if not param.IsFixed(): nfree += 1
    print('Number of free parameters: {}'.format(nfree))
    ndof = nbins - nfree
    chi2ndof = chi2/ndof
    print('Chi2/ndof: {:.2f} for ndof: {}'.format(chi2ndof, ndof))

    plt.figure()
    image = plt.imshow(pulls, origin='lower', extent=(m12.getLowerLimit(),m12.getUpperLimit(),m12.getLowerLimit(),m12.getUpperLimit()), vmin=-4, vmax=4, cmap='RdYlBu')
    plt.xlabel('m12')
    plt.ylabel('m13')
    plt.title('Pulls')
    plt.colorbar(image)
    plt.savefig('Tddp_pulls'+string+'.pdf')

    # Plot Dalitz variables 
    if bkg_flag:
        histsig, _ = np.histogram(m12val, weights=array[1], bins=175, range=(m12.getLowerLimit(),m12.getUpperLimit()))
        histbkg, _ = np.histogram(m12val, weights=array[2], bins=175, range=(m12.getLowerLimit(),m12.getUpperLimit()))
        histsig = histsig*totalsigprob/totalPdf_sig
        histbkg = histbkg*totalbkprob/totalPdf_bkg
        hist = np.add(histsig, histbkg)
        #hist, _ = np.histogram(m12val, weights=array[0], bins=175, range=(m12.getLowerLimit(),m12.getUpperLimit()))
        #hist = hist*dataset.getNumEvents()/totalPdf
    else: 
        hist, _ = np.histogram(m12val, weights=array[0], bins=175, range=(m12.getLowerLimit(),m12.getUpperLimit()))
        hist = hist*dataset.getNumEvents()/totalPdf
    hist_dat, bin_edges = np.histogram(m12dat, bins=175, range=(m12.getLowerLimit(),m12.getUpperLimit()))

    hist_tmp, _ = np.histogram(m12dat, bins=175, range=(m12.getLowerLimit(),m12.getUpperLimit()))
    bin_centre = (bin_edges[:-1] + bin_edges[1:])/2
    binwidth = bin_edges[1]-bin_edges[0]

    pulls = (hist_dat - hist)/np.sqrt(hist_dat)

    _, ax = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[3, 1]})
    ax[0].plot(bin_centre, hist, zorder=2, color='b', label='Fit')
    if bkg_flag:
        ax[0].plot(bin_centre, histsig, zorder=2, color='r', label='Signal')
        ax[0].plot(bin_centre, histbkg, zorder=2, color='g', label='Background')    
    ax[0].errorbar(bin_centre, hist_dat, yerr=np.sqrt(hist_dat), capsize=2, fmt='o', color='k', ms=4, zorder=1, elinewidth=1, label='Data')
    ax[0].set_xlim(m12.getLowerLimit(),m12.getUpperLimit())
    ax[0].set_ylabel('Candidates (normalised)')
    ax[1].bar(bin_centre, pulls, color='k', width=binwidth)
    pulls_nozero = np.trim_zeros(pulls)
    npull, binspull, _ = ax[1].hist(pulls_nozero, orientation='horizontal', density=True, bins=40, alpha=0.5, color='k', range=(-4,4))
    centers = (0.5*(binspull[1:]+binspull[:-1]))
    pars, cov = curve_fit(lambda pulls_nozero, mu, sig : norm.pdf(pulls_nozero, loc=mu, scale=sig), centers, npull, p0=[0,1])  
    par_ers = np.sqrt(np.diag(cov))
    gaus = norm.pdf(centers, *pars)
    ax[1].plot(gaus, centers, color='r')
    ax[0].text(0.4,0.8,'mu: {} +/- {} \nsigma: {} +/- {}'.format(round(pars[0],3),round(par_ers[0],3), round(pars[1],3), round(par_ers[1],3)), bbox=dict(facecolor='w', alpha=0.5), transform=ax[0].transAxes)
    ax[0].legend(loc='upper right')
    ax[1].set_ylim(-4,4)
    ax[1].set_xlim(m12.getLowerLimit(),m12.getUpperLimit())
    ax[1].set_xlabel(r'$m^2_{12}$ [GeV]')
    ax[1].set_ylabel('Pulls')
    plt.tight_layout()
    plt.savefig('m12_test'+string+'.pdf')

    # Log plot
    _, ax = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[3, 1]})
    ax[0].set_yscale('log')
    ax[0].plot(bin_centre, hist, zorder=2, color='b', label='Fit')
    if bkg_flag:
        ax[0].plot(bin_centre, histsig, zorder=2, color='r', label='Signal')
        ax[0].plot(bin_centre, histbkg, zorder=2, color='g', label='Background')    
    ax[0].errorbar(bin_centre, hist_dat, yerr=np.sqrt(hist_dat), capsize=2, fmt='o', color='k', ms=4, zorder=1, elinewidth=1, label='Data')
    ax[0].set_xlim(m12.getLowerLimit(),m12.getUpperLimit())
    ax[0].set_ylabel('Candidates (normalised)')
    ax[1].bar(bin_centre, pulls, color='k', width=binwidth)
    ax[1].hist(pulls_nozero, orientation='horizontal', density=True, bins=40, alpha=0.5, color='k', range=(-4,4))
    ax[1].plot(gaus, centers, color='r')
    ax[0].text(0.4,0.8,'mu: {} +/- {} \nsigma: {} +/- {}'.format(round(pars[0],3),round(par_ers[0],3), round(pars[1],3), round(par_ers[1],3)), bbox=dict(facecolor='w', alpha=0.5), transform=ax[0].transAxes)
    ax[0].legend(loc='best')
    ax[1].set_ylim(-4,4)
    ax[1].set_xlim(m12.getLowerLimit(),m12.getUpperLimit())
    ax[1].set_xlabel(r'$m^2_{12}$ [GeV]')
    ax[1].set_ylabel('Pulls')
    plt.tight_layout()

    plt.savefig('m12_log'+string+'.pdf')


    if bkg_flag:
        histsig, _ = np.histogram(m13val, weights=array[1], bins=175, range=(m13.getLowerLimit(),m13.getUpperLimit()))
        histbkg, _ = np.histogram(m13val, weights=array[2], bins=175, range=(m13.getLowerLimit(),m13.getUpperLimit()))
        histsig = histsig*totalsigprob/totalPdf_sig
        histbkg = histbkg*totalbkprob/totalPdf_bkg
        hist = np.add(histsig, histbkg)
        #hist = np.add(sigfrac*histsig, (1-sigfrac)*histbkg)
        #histsig = np.reshape(histsig, (175, 5)).mean(axis=1)
        #histbkg = np.reshape(histbkg, (175, 5)).mean(axis=1)
    else: 
        hist, _ = np.histogram(m13val, weights=array[0], bins=175, range=(m13.getLowerLimit(),m13.getUpperLimit()))
        hist = hist*dataset.getNumEvents()/totalPdf

    hist_dat, bin_edges = np.histogram(m13dat, bins=175, range=(m13.getLowerLimit(),m13.getUpperLimit()))
    #hist_dat = hist_dat*dataset.getNumEvents()/totalPdf
    hist_tmp, _ = np.histogram(m13dat, bins=175, range=(m13.getLowerLimit(),m13.getUpperLimit()))
    bin_centre = (bin_edges[:-1] + bin_edges[1:])/2
    binwidth = bin_edges[1]-bin_edges[0]

    yerr = np.sqrt(hist_tmp)/(hist_tmp.sum()*binwidth) # sqrt(n) / integral
    #pulls = np.nan_to_num((hist - hist_dat)/yerr, posinf=0)
    pulls = (hist_dat - hist)/np.sqrt(hist_dat)

    _, ax = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[3, 1]})
    ax[0].plot(bin_centre, hist, zorder=2, color='b', label='Fit')
    if bkg_flag:
        ax[0].plot(bin_centre, histsig, zorder=2, color='r', label='Signal')
        ax[0].plot(bin_centre, histbkg, zorder=2, color='g', label='Background')    
    ax[0].errorbar(bin_centre, hist_dat, yerr=np.sqrt(hist_dat), capsize=2, fmt='o', color='k', ms=4, zorder=1, elinewidth=1, label='Data')
    ax[0].set_xlim(m13.getLowerLimit(),m13.getUpperLimit())
    ax[0].set_ylabel('Candidates (normalised)')
    ax[1].bar(bin_centre, pulls, color='k', width=binwidth)
    pulls_nozero = np.trim_zeros(pulls)
    npull, binspull, _ = ax[1].hist(pulls_nozero, orientation='horizontal', density=True, bins=40, alpha=0.5, color='k', range=(-4,4))
    centers = (0.5*(binspull[1:]+binspull[:-1]))
    pars, cov = curve_fit(lambda pulls_nozero, mu, sig : norm.pdf(pulls_nozero, loc=mu, scale=sig), centers, npull, p0=[0,1])  
    par_ers = np.sqrt(np.diag(cov))
    gaus = norm.pdf(centers, *pars)
    ax[1].plot(gaus, centers, color='r')
    ax[0].text(0.3,0.8,'mu: {} +/- {} \nsigma: {} +/- {}'.format(round(pars[0],3),round(par_ers[0],3), round(pars[1],3), round(par_ers[1],3)), bbox=dict(facecolor='w', alpha=0.5), transform=ax[0].transAxes)
    ax[0].legend(loc='best')
    ax[1].set_ylim(-4,4)
    ax[1].set_xlim(m13.getLowerLimit(),m13.getUpperLimit())
    ax[1].set_xlabel(r'$m^2_{13}$ [GeV]')
    ax[1].set_ylabel('Pulls')
    plt.tight_layout()
    plt.savefig('m13_test'+string+'.pdf')

    # Log plot
    _, ax = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[3, 1]})
    ax[0].set_yscale('log')
    ax[0].plot(bin_centre, hist, zorder=2, color='b', label='Fit')
    if bkg_flag:
        ax[0].plot(bin_centre, histsig, zorder=2, color='r', label='Signal')
        ax[0].plot(bin_centre, histbkg, zorder=2, color='g', label='Background')    
    ax[0].errorbar(bin_centre, hist_dat, yerr=np.sqrt(hist_dat), capsize=2, fmt='o', color='k', ms=4, zorder=1, elinewidth=1, label='Data')
    ax[0].set_xlim(m13.getLowerLimit(),m13.getUpperLimit())
    ax[0].set_ylabel('Candidates (normalised)')
    ax[1].bar(bin_centre, pulls, color='k', width=binwidth)
    ax[1].hist(pulls_nozero, orientation='horizontal', density=True, bins=40, alpha=0.5, color='k', range=(-4,4))
    ax[1].plot(gaus, centers, color='r')
    ax[0].text(0.3,0.8,'mu: {} +/- {} \nsigma: {} +/- {}'.format(round(pars[0],3),round(par_ers[0],3), round(pars[1],3), round(par_ers[1],3)), bbox=dict(facecolor='w', alpha=0.5), transform=ax[0].transAxes)
    ax[0].legend(loc='best')
    ax[1].set_ylim(-4,4)
    ax[1].set_xlim(m13.getLowerLimit(),m13.getUpperLimit())
    ax[1].set_xlabel(r'$m^2_{13}$ [GeV]')
    ax[1].set_ylabel('Pulls')
    plt.tight_layout()

    plt.savefig('m13_log'+string+'.pdf')

    if bkg_flag:
        histsig, _ = np.histogram(m23val, weights=array[1], bins=175, range=(0,3.5))
        histbkg, _ = np.histogram(m23val, weights=array[2], bins=175, range=(0,3.5))
        histsig = histsig*totalsigprob/totalPdf_sig
        histbkg = histbkg*totalbkprob/totalPdf_bkg
        hist = np.add(histsig, histbkg)
        #hist = np.add(sigfrac*histsig, (1-sigfrac)*histbkg)
        #histsig = np.reshape(histsig, (175, 5)).mean(axis=1)
        #histbkg = np.reshape(histbkg, (175, 5)).mean(axis=1)
    else: 
        hist, _ = np.histogram(m23val, weights=array[0], bins=175, range=(0,3.5))
        hist = hist*dataset.getNumEvents()/totalPdf

    hist_dat, bin_edges = np.histogram(m23dat, bins=175, range=(0,3.5))
    #hist_dat = hist_dat*dataset.getNumEvents()/totalPdf
    hist_tmp, _ = np.histogram(m23dat, bins=175, range=(0,3.5))
    bin_centre = (bin_edges[:-1] + bin_edges[1:])/2
    binwidth = bin_edges[1]-bin_edges[0]

    yerr = np.sqrt(hist_tmp)/(hist_tmp.sum()*binwidth) # sqrt(n) / integral
    #pulls = np.nan_to_num((hist - hist_dat)/yerr, posinf=0)
    pulls = (hist_dat - hist)/np.sqrt(hist_dat)

    _, ax = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[3, 1]})
    ax[0].plot(bin_centre, hist, zorder=2, color='b', label='Fit')
    if bkg_flag:
        ax[0].plot(bin_centre, histsig, zorder=2, color='r', label='Signal')
        ax[0].plot(bin_centre, histbkg, zorder=2, color='g', label='Background') 
    ax[0].errorbar(bin_centre, hist_dat, yerr=np.sqrt(hist_dat), capsize=2, fmt='o', color='k', ms=4, zorder=1, elinewidth=1, label='Data')
    ax[0].set_xlim(0,2)
    ax[0].set_ylabel('Candidates (normalised)')
    ax[1].bar(bin_centre, pulls, color='k', width=binwidth)
    pulls_nozero = np.trim_zeros(pulls)
    npull, binspull, _ = ax[1].hist(pulls_nozero, orientation='horizontal', density=True, bins=40, alpha=0.5, color='k', range=(-4,4))
    centers = (0.5*(binspull[1:]+binspull[:-1]))
    pars, cov = curve_fit(lambda pulls_nozero, mu, sig : norm.pdf(pulls_nozero, loc=mu, scale=sig), centers, npull, p0=[0,1])  
    par_ers = np.sqrt(np.diag(cov))
    gaus = norm.pdf(centers, *pars)
    ax[1].plot(gaus, centers, color='r')
    ax[0].text(0.4,0.8,'mu: {} +/- {} \nsigma: {} +/- {}'.format(round(pars[0],3),round(par_ers[0],3), round(pars[1],3), round(par_ers[1],3)), bbox=dict(facecolor='w', alpha=0.5), transform=ax[0].transAxes)
    ax[0].legend(loc='best')
    ax[1].set_ylim(-4,4)
    ax[1].set_xlim(0,2)
    ax[1].set_xlabel(r'$m^2_{23}$ [GeV]')
    ax[1].set_ylabel('Pulls')
    plt.tight_layout()
    plt.savefig('m23_test'+string+'.pdf')

    # Log plot
    _, ax = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[3, 1]})
    ax[0].set_yscale('log')
    ax[0].plot(bin_centre, hist, zorder=2, color='b', label='Fit')
    if bkg_flag:
        ax[0].plot(bin_centre, histsig, zorder=2, color='r', label='Signal')
        ax[0].plot(bin_centre, histbkg, zorder=2, color='g', label='Background')    
    ax[0].errorbar(bin_centre, hist_dat, yerr=np.sqrt(hist_dat), capsize=2, fmt='o', color='k', ms=4, zorder=1, elinewidth=1, label='Data')
    ax[0].set_xlim(0,2)
    ax[0].set_ylabel('Candidates (normalised)')
    ax[0].set_xlim(0,2)
    ax[1].bar(bin_centre, pulls, color='k', width=binwidth)
    ax[1].hist(pulls_nozero, orientation='horizontal', density=True, bins=40, alpha=0.5, color='k', range=(-4,4))
    ax[1].plot(gaus, centers, color='r')
    ax[0].text(0.4,0.8,'mu: {} +/- {} \nsigma: {} +/- {}'.format(round(pars[0],3),round(par_ers[0],3), round(pars[1],3), round(par_ers[1],3)), bbox=dict(facecolor='w', alpha=0.5), transform=ax[0].transAxes)
    ax[0].legend(loc='best')
    ax[1].set_ylim(-4,4)
    ax[1].set_xlim(0,2)
    ax[1].set_xlabel(r'$m^2_{23}$ [GeV]')
    ax[1].set_ylabel('Pulls')
    plt.tight_layout()

    plt.savefig('m23_log'+string+'.pdf')

    # For decay time plots
    count = 0
    sigma.setValue(0.1)
    mistag.setValue(0.005)
    charmtag.setValue(1)
    if bkg_flag: sigprob.setValue(sigfrac)
    else: sigprob.setValue(1)

    for i in range(80):
        m12.setValue(m12.getLowerLimit() + (m12.getUpperLimit() - m12.getLowerLimit()) * (i + 0.5) / 80)
        for j in range(80):
            m13.setValue(m13.getLowerLimit() + (m13.getUpperLimit() - m13.getLowerLimit()) * (j + 0.5) / 80)
            if inDalitz(m12.getValue(), m13.getValue()):
                for k in range(175*5):
                    dtime.setValue(dtime.getLowerLimit() + (dtime.getUpperLimit() - dtime.getLowerLimit()) * (k + 0.5) / 875)
                    dtval.append(dtime.getValue())
                    if np.random.uniform() < d0fraction: charmtag.setValue(1)
                    else: charmtag.setValue(-1)
                    eventNumber.setValue(count)
                    data_dtime.addEvent()
                    count += 1

    data_dtime.setValueForAllEvents(sigma)
    data_dtime.setValueForAllEvents(mistag)
    data_dtime.setValueForAllEvents(charmtag)
    data_dtime.setValueForAllEvents(sigprob)
    print(data_dtime.getNumEvents())

    print(data_dtime.getObservables())
    print(total.getObservables())

    signal.setDataSize(data_dtime.getNumEvents(), 8)

    total.setData(data_dtime)
    vals2 = total.getCompProbsAtDataPoints()
    array = np.array(vals2)

    totalPdf = np.sum(array[0])
    totalPdf_sig = np.sum(array[1])
    if bkg_flag: totalPdf_bkg = np.sum(array[2])

    # Plot Decay time
    if bkg_flag:
        histsig, _ = np.histogram(dtval, weights=array[1], bins=175, range=(dtime.getLowerLimit(),dtime.getUpperLimit()))
        histbkg, _ = np.histogram(dtval, weights=array[2], bins=175, range=(dtime.getLowerLimit(),dtime.getUpperLimit()))
        histsig = histsig*totalsigprob/totalPdf_sig
        histbkg = histbkg*totalbkprob/totalPdf_bkg
        hist = np.add(histsig, histbkg)
        #hist = np.add(sigfrac*histsig, (1-sigfrac)*histbkg)
        #histsig = np.reshape(histsig, (175, 3)).mean(axis=1)
        #histbkg = np.reshape(histbkg, (175, 3)).mean(axis=1)

    else:
        hist, _ = np.histogram(dtval, weights=array[0], bins=175, range=(dtime.getLowerLimit(),dtime.getUpperLimit()))
        hist = hist*dataset.getNumEvents()/totalPdf

    hist_dat, bin_edges = np.histogram(dtdat, bins=175, range=(dtime.getLowerLimit(),dtime.getUpperLimit()))
    #hist_dat = hist_dat*dataset.getNumEvents()/totalPdf
    bin_centre = (bin_edges[:-1] + bin_edges[1:])/2
    hist_tmp, _ = np.histogram(dtdat, bins=175, range=(dtime.getLowerLimit(),dtime.getUpperLimit()))
    binwidth = bin_edges[1]-bin_edges[0]
    print('Bin width: ', binwidth)
    yerr = np.sqrt(hist_tmp)/(hist_tmp.sum()*binwidth) # sqrt(n) / integral
    #pulls = np.nan_to_num((hist - hist_dat)/yerr, posinf=0, neginf=0)
    pulls = (hist_dat - hist)/np.sqrt(hist_dat)

    _, ax = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[3, 1]})
    ax[0].plot(bin_centre, hist, zorder=2, color='b', label='Fit')
    if bkg_flag:    
        ax[0].plot(bin_centre, histsig, zorder=2, color='r', label='Signal')
        ax[0].plot(bin_centre, histbkg, zorder=2, color='g', label='Background')
    ax[0].errorbar(bin_centre, hist_dat, yerr=np.sqrt(hist_dat), capsize=2, fmt='o', color='k', ms=4, zorder=1, elinewidth=1, label='Data')
    ax[0].set_xlim(dtime.getLowerLimit(),dtime.getUpperLimit())
    ax[0].set_ylabel('Candidates (normalised)')
    ax[0].set_ylim(bottom=0)
    ax[1].bar(bin_centre, pulls, color='k', width=binwidth)
    pulls_nozero = np.trim_zeros(pulls)
    npull, binspull, _ = ax[1].hist(pulls_nozero, orientation='horizontal', density=True, bins=40, alpha=0.5, color='k', bottom=-1, range=(-4,4))
    centers = (0.5*(binspull[1:]+binspull[:-1]))
    pars, cov = curve_fit(lambda pulls_nozero, mu, sig : norm.pdf(pulls_nozero, loc=mu, scale=sig), centers, npull, p0=[0,1])
    par_ers = np.sqrt(np.diag(cov))
    gaus = norm.pdf(centers, *pars)
    ax[1].plot(gaus-1, centers, color='r')
    ax[0].text(0.3,0.8,'mu: {} +/- {} \nsigma: {} +/- {}'.format(round(pars[0],3),round(par_ers[0],3), round(pars[1],3), round(par_ers[1],3)), bbox=dict(facecolor='w', alpha=0.5), transform=ax[0].transAxes)
    ax[0].legend(loc='best')
    ax[1].set_ylim(-4,4)
    ax[1].set_xlim(dtime.getLowerLimit(),dtime.getUpperLimit())
    ax[1].set_xlabel('Decay Time (ps)')
    ax[1].set_ylabel('Pulls')
    plt.tight_layout()
    plt.savefig('dtime_test'+string+'.pdf')

    # Log plot
    _, ax = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[3, 1]})
    #ax[0].semilogy(nonposy='clip')
    ax[0].set_yscale('log', nonposy='clip')
    ax[0].plot(bin_centre, hist, zorder=2, color='b', label='Fit')
    if bkg_flag:
        ax[0].plot(bin_centre, histsig, zorder=2, color='r', label='Signal')
        ax[0].plot(bin_centre, histbkg, zorder=2, color='g', label='Background')
    ax[0].errorbar(bin_centre, hist_dat, yerr=np.sqrt(hist_dat), capsize=2, fmt='o', color='k', ms=4, zorder=1, elinewidth=1, label='Data')
    ax[0].set_xlim(dtime.getLowerLimit(),dtime.getUpperLimit())
    ax[0].set_ylabel('Candidates (normalised)')
    #ax[0].set_ylim(bottom=1e-8)
    ax[1].bar(bin_centre, pulls, color='k', width=binwidth)
    ax[1].hist(pulls_nozero, orientation='horizontal', density=True, bins=40, alpha=0.5, color='k', bottom=-1, range=(-4,4))
    ax[1].plot(gaus-1, centers, color='r')
    ax[0].text(0.3,0.8,'mu: {} +/- {} \nsigma: {} +/- {}'.format(round(pars[0],3),round(par_ers[0],3), round(pars[1],3), round(par_ers[1],3)), bbox=dict(facecolor='w', alpha=0.5), transform=ax[0].transAxes)
    ax[0].legend(loc='best')
    ax[1].set_ylim(-4,4)
    ax[1].set_xlim(dtime.getLowerLimit(),dtime.getUpperLimit())
    ax[1].set_xlabel('Decay time (ps)')
    ax[1].set_ylabel('Pulls')
    plt.tight_layout()

    plt.savefig('dtime_log'+string+'.pdf')

def Cholesky(Corr_matrix):
    npars = Corr_matrix.shape[0]
    print(npars)
    L = np.zeros(npars*npars)
    Corr_matrix = Corr_matrix.reshape(npars*npars)

    for i in range(npars):
        for j in range(i+1):
            s = 0
            for k in range(j):
                s += L[i*npars+k]*L[j*npars+k]
            if i==j: L[i*npars+j] = sqrt(Corr_matrix[i*npars+i] - s)
            else: L[i*npars+j] = 1/L[j*npars+j] * (Corr_matrix[i*npars+j] - s)
            
    print(L)
    return L

def resamplePSefficiency(efficiency):
    corr_matrix = np.load('../PSEfficiency/correlation_matrix_{}_{}.npy'.format(mode, year))
    npars = corr_matrix.shape[0]
    Cholesky_matrix = Cholesky(corr_matrix)

    path = '../PSEfficiency/'
    eff_params = pd.read_csv(path+'Efficiency_params_'+mode+'_'+year+'.csv')
    c = eff_params.params
    errors = eff_params.errors

    # Generate a set of uncorrelated parameters
    a = np.random.normal(size=npars)
    a_shift = [(a[0]*Cholesky_matrix[npars*n+0] + a[1]*Cholesky_matrix[npars*n+1] + a[2]*Cholesky_matrix[npars*n+2] + a[3]*Cholesky_matrix[npars*n+3] + a[4]*Cholesky_matrix[npars*n+4] + a[5]*Cholesky_matrix[npars*n+5] + a[6]*Cholesky_matrix[npars*n+6] + a[7]*Cholesky_matrix[npars*n+7]) for n in range(npars)]

    print(c)
    new_coeffs = a_shift*errors + c
    print(new_coeffs)

    for i, param in enumerate(efficiency.getParameters()):
        print('Resampling parameter {} original value: {}, new value: {}'.format(param.getName(), c[i], new_coeffs[i]))
        param.setValue(new_coeffs[i])

def resample_LASS():
    cov_matrix = np.load('cov_matrix_LASS_SingleTag_D0ToKsPiPiLL_2016.npy')
    npars = cov_matrix.shape[0]

    LASS_params = [K0starm_1430_amp_real, K0starm_1430_amp_imag, K0starp_1430_amp_real, K0starp_1430_amp_imag,
                   a, r, phiR, F, phiF]
    params_csv = pd.read_csv('Fit_params_LASS_Step2_SingleTag_D0ToKsPiPiLL_2016.csv')
    LASS_params_csv = [params_csv.K0starm_1430_amp_real, params_csv.K0starm_1430_amp_imag,
                       params_csv.K0starp_1430_amp_real, params_csv.K0starp_1430_amp_imag, 
                       params_csv.a, params_csv.r, params_csv.phiR, params_csv.F, params_csv.phiF]
    LASS_params_errors = [params_csv.K0starm_1430_amp_real_err, params_csv.K0starm_1430_amp_imag_err,
                          params_csv.K0starp_1430_amp_real_err, params_csv.K0starp_1430_amp_imag_err,
                          params_csv.a_err, params_csv.r_err, params_csv.phiR_err, params_csv.F_err, 
                          params_csv.phiF_err]

    params = np.array(LASS_params_csv).reshape(npars)
    errors = np.array(LASS_params_errors).reshape(npars)

    corr_matrix = np.zeros((npars, npars))
    for i in range(npars):
        for j in range(npars):
            corr_matrix[i][j] = cov_matrix[i][j]/(errors[i]*errors[j])

    Cholesky_matrix = Cholesky(corr_matrix)


    # Generate a set of uncorrelated parameters
    rand = np.random.normal(size=npars)
    shift = np.zeros(npars)

    # Multiply by Cholesky decomposed matrix to get correlated parameters
    for i in range(npars):
        for j in range(npars):
            shift[i] += rand[j]*Cholesky_matrix[npars*i+j]

    print(shift)

    new_coeffs = shift*errors + params
    print(new_coeffs)

    for i, param in enumerate(LASS_params):
        if 'K0star' in param.getName(): continue # Do not resample amplitudes
        print('Resampling parameter {} original value: {}, new value: {}'.format(param.getName(), params[i], new_coeffs[i]))
        param.setValue(new_coeffs[i])



def main():
    print_goofit_info()
    years = [2016, 2017, 2018]
    modes = ['SingleTag_D0ToKsPiPiLL', 'SingleTag_D0ToKsPiPiDD', 'DoubleTag_D0ToKsPiPiLL', 'DoubleTag_D0ToKsPiPiDD']

    mes_radius    = 1.5
    mother_radius = 5.0
    BW_syst = False
    if BW_syst:
        mes_radius    = np.random.uniform(low=1.0, high=1.8, size=None)
        mother_radius = np.random.uniform(low=4.2, high=5.0, size=None)

    decayInfo = makeDecayInfo(mes_radius, mother_radius)


    m12.setNumBins(175*4)
    m13.setNumBins(175*4)
    dtime.setNumBins(175)
    sigma.setNumBins(1)
    sigprob.setNumBins(10)
    charmtag.setNumBins(2)
    mistag.setNumBins(2)

    sigma.setValue(0.1)

    datasets = getDatasets(years, modes)

    arrays = [a.to_matrix() for a in datasets.values()]
    tot_matrix = np.concatenate(arrays, axis=1)
    #tot_matrix[3] = np.arange(tot_matrix.shape[1])
    print(tot_matrix)

    data_combined = UnbinnedDataSet(m12, m13, dtime, sigprob, sigma, mistag, charmtag, category, eventNumber)
    data_combined.from_matrix(tot_matrix) # check this works
    print('Number of events in combined dataset: {}'.format(data_combined.getNumEvents()))


    signal_pdfs = {}
    bkg_pdfs = {}
    total_pdfs = {}
    bkg1_pdfs = {}
    bkg2_pdfs = {}
    efficency_pdfs = {}

    dalitz_smoothing = Variable('dalitz_smoothing', 0.5)
    dtime_smoothing = Variable('dtime_smoothing', 0.5)
    bghist_dalitz = BinnedDataSet(m12, m13)
    bghist_dtime = BinnedDataSet(dtime)

    for year in years:
      for mode in modes:
        res_params = pd.read_csv('../DTResolution/DTResolution_fit_parameters.csv')
        select_year = res_params['year'] == int(year)
        select_mode = res_params['mode'] == mode
        select_res_params = res_params.loc[select_year & select_mode]
        select_res_params.reset_index(inplace=True)

        # Resolution parameters
        coreFrac  = Variable(f'coreFrac_{mode}_{year}', select_res_params['f'][0])
        tailFrac  = Variable(f'tailFrac_{mode}_{year}', select_res_params['g'][0])
        mean      = Variable(f'mean_{mode}_{year}', select_res_params['mean'][0]/(sigma.getValue()))
        coreScale = Variable(f'coreScale_{mode}_{year}', select_res_params['sigma'][0]/(sigma.getValue())) # Check parameters
        tailScale = Variable(f'tailScale_{mode}_{year}', select_res_params['sigma1'][0]/(sigma.getValue()))
        outScale  = Variable(f'outScale_{mode}_{year}', select_res_params['sigma2'][0]/(sigma.getValue()))

        acc_params = pd.read_csv('../DTAcceptance/DTAcceptance_fit_parameters.csv')
        select_year = acc_params['year'] == int(year)
        select_mode = acc_params['mode'] == mode
        select_acc_params = acc_params.loc[select_year & select_mode]
        select_acc_params.reset_index(inplace=True)
        print(select_acc_params['slope'][0])

        selBias   = Variable(f'selBias_{mode}_{year}', select_acc_params['slope'][0])

        resolution = ThreeGaussResolution(coreFrac, 
                                          tailFrac, 
                                          mean, 
                                          coreScale, 
                                          mean, 
                                          tailScale, 
                                          mean, 
                                          outScale,
                                          selBias)

        yearmode = f'{year}_{mode}'
        efficiency = makeEfficiency(mode, year)
        efficency_pdfs[yearmode] = efficiency
        
        signal_pdf = Amp3Body_TD(f'signalPDF_{yearmode}', dtime, sigma, m12, m13, eventNumber, decayInfo, resolution, efficiency, mistag, charmtag)
        signal_pdfs[yearmode] = signal_pdf

        bghist_dalitz = BinnedDataSet(m12, m13)
        bghist_dtime = BinnedDataSet(dtime)
        filename = f'root://eoslhcb.cern.ch//eos/lhcb/user/m/mhilton/KSPiPi-ntuples/tuples-BDT/sWeight_{mode}_{year}.root'
        fillBGhist(filename, bghist_dalitz, bghist_dtime)
        BGPdf = SmoothHistogramPdf(f'BGPdf_{yearmode}', bghist_dalitz, dalitz_smoothing)
        BGPdf_dtime_hist = SmoothHistogramPdf(f'BGPdf_dtime_hist_{yearmode}', bghist_dtime, dtime_smoothing)
        bkg1_pdfs[yearmode] = BGPdf
        bkg2_pdfs[yearmode] = BGPdf_dtime_hist
        BGpdfList = [BGPdf, BGPdf_dtime_hist]
        BGPdf3D = ProdPdf(f'BGPdf3D_{yearmode}', [BGPdf, BGPdf_dtime_hist])
        bkg_pdfs[yearmode] = BGPdf3D

        total_pdfs[yearmode] = EventWeightedAddPdf(f'total_{yearmode}', [sigprob], [signal_pdf, BGPdf3D])

    stepFunction = BinTransformPdf('stepFunction', [category], [-0.5], 
                                    [1], [len(datasets)])
    finalPDF_comb = MappedPdf("finalPDF_comb", stepFunction, list(total_pdfs.values()))

    #params = pd.read_csv('Timeintegrated_Fits/Fit_params_timeintegrated_{}_{}.csv'.format(mode, year))
    params = pd.read_csv('Fit_params_SingleTag_D0ToKsPiPiLL.csv')
    for column in params.columns:
        for param in finalPDF_comb.getParameters():
            if param.getName() == column:
                param.setValue(params[column])
    params['xcp'] = xcp.getValue()
    params['ycp'] = ycp.getValue()
    params['tau'] = tau.getValue()
    params['deltax'] = deltax.getValue()
    params['deltay'] = deltay.getValue()

    for param in finalPDF_comb.getParameters():
        if 'amp_real' in param.getName() or 'amp_imag' in param.getName():
            param.setFixed(True)
        if 'mass' in param.getName() or 'width' in param.getName():
            param.setFixed(True)
        if 'beta' in param.getName() or 'f_prod' in param.getName():
            param.setFixed(True) # Fix K-matrix parameters

    rho_770_amp_imag.setFixed(True)
    rho_770_amp_real.setFixed(True)


    m12.setNumBins(175*8)
    m13.setNumBins(175*8)

    np.random.seed(123)

    blindseedx = np.random.normal(loc=0, scale=0.005)
    blindseedy = np.random.normal(loc=0, scale=0.005)

    #xcp.setBlind(np.random.normal(loc=0, scale=0.005))
    #ycp.setBlind(np.random.normal(loc=0, scale=0.005))

    xcp.setBlind(blindseedx)
    ycp.setBlind(blindseedy)
    
    finalPDF_comb.setData(data_combined)
    offset = 0 
    for year in years:
      for mode in modes:
        signal_pdfs[f'{year}_{mode}'].setDataSize(datasets[f'{year}_{mode}'].getNumEvents(), 9, offset)
        offset = datasets[f'{year}_{mode}'].getNumEvents() + offset

    #resamplePSefficiency(efficiency)

    LASS_resmaple = False
    if LASS_resmaple: resample_LASS()

    fitman = FitManager(finalPDF_comb)
    fitman.setMaxCalls(320000)
    fitman.setMinos(False)
    print('Running fit...')
    func_min = fitman.fit()
    
    FCNval = func_min.Fval()
    #reset cache counter before plotting
    Amp3Body_TD.resetCacheCounter()

    for year in years:
      for mode in modes:
        makePlots(total_pdfs[f'{year}_{mode}'], signal_pdfs[f'{year}_{mode}'], datasets[f'{year}_{mode}'], 1, f'_{mode}_{year}')





if __name__ == '__main__':
    sys.exit(main())
