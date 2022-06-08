#!/usr/bin/env python
# -*- coding: utf-8 -*-

# In[1]:


from math import sqrt

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

import goofit
from goofit import *

# In[2]:


# In[3]:


print_goofit_info()


# ## Setting constants

# In[4]:


B_MASS = 5.27934  # GeV
k_MASS = 0.493677  # GeV
smin = 4 * k_MASS * k_MASS  # GeV^2
smax = (B_MASS - k_MASS) ** 2  # GeV^2
nBins = 200


# In[5]:


s12 = Observable("s12", smin, smax)
s13 = Observable("s13", smin, smax)
eventNumber = EventNumber("eventNumber")

s12.setNumBins(nBins)
s13.setNumBins(nBins)

Mother_Mass = Variable("Mother_Mass", B_MASS)
Daughter1_Mass = Variable("Daughter1_Mass", k_MASS)
Daughter2_Mass = Variable("Daughter2_Mass", k_MASS)
Daughter3_Mass = Variable("Daughter3_Mass", k_MASS)

B2KKK = DecayInfo3()
B2KKK.motherMass = B_MASS
B2KKK.daug1Mass = k_MASS
B2KKK.daug2Mass = k_MASS
B2KKK.daug3Mass = k_MASS
B2KKK.meson_radius = 1.5
B2KKK.mother_meson_radius = 5.0

constantOne = Variable("constantOne", 1)
constantZero = Variable("constantZero", 0)


# ## Phase space info

# In[6]:


print(s12)
print(s13)


# ## Signal PDF

# In[7]:


def makePolyEff():
    observables = (s12, s13)
    offsets = (
        constantZero,
        constantZero,
    )
    coefficients = (constantOne,)

    return PolynomialPdf("constantEff", observables, coefficients, offsets, 0)


# In[8]:


def makesignal(eff):
    phi_re = Variable("phi_re", 1.0)
    phi_im = Variable("phi_im", 0.0)
    phi_mass = Variable("phi_mass", 1.019460)
    phi_width = Variable("phi_width", 0.004247)

    f2p_re = Variable("f2p_re", 0.0, 0.01, 0, 0)
    f2p_im = Variable("f2p_im", 1.0, 0.01, 0, 0)
    f2p_mass = Variable("f2p_mass", 1.5245)
    f2p_width = Variable("f2p_width", 0.0733158)

    nr_re = Variable("nr_re", 1.0, 0.01, 0, 0)
    nr_im = Variable("nr_im", 0.0, 0.01, 0, 0)

    phi = Resonances.RBW(
        "phi", phi_re, phi_im, phi_mass, phi_width, 1, PAIR_12, True, True
    )
    f2p = Resonances.RBW(
        "f2", f2p_re, f2p_im, f2p_mass, f2p_width, 2, PAIR_12, True, True
    )
    nr = Resonances.NonRes("nr", nr_re, nr_im)

    B2KKK.resonances = (phi, f2p, nr)

    d = Amp3Body("signalPDF", s12, s13, eventNumber, B2KKK, eff)

    return d


# ## Signal with veto

# In[9]:


def veto():
    DMass = 1.86966  # GeV
    veto_12 = VetoInfo(
        Variable("veto_min_12", (DMass - k_MASS) ** 2),
        Variable("veto_max_12", smax),
        PAIR_12,
    )
    veto_13 = VetoInfo(
        Variable("veto_min_13", (DMass - k_MASS) ** 2),
        Variable("veto_max_13", smax),
        PAIR_13,
    )
    vetos = [veto_12, veto_13]
    return DalitzVetoPdf(
        "veto",
        s12,
        s13,
        Mother_Mass,
        Daughter1_Mass,
        Daughter2_Mass,
        Daughter3_Mass,
        vetos,
    )


# In[10]:


def makesignalwithveto(eff):
    phi_re = Variable("phi_re", 1.0)
    phi_im = Variable("phi_im", 0.0)
    phi_mass = Variable("phi_mass", 1.019460)
    phi_width = Variable("phi_width", 0.004247)

    f2p_re = Variable("f2p_re", 0.0, 0.01, 0, 0)
    f2p_im = Variable("f2p_im", 1.0, 0.01, 0, 0)
    f2p_mass = Variable("f2p_mass", 1.5245)
    f2p_width = Variable("f2p_width", 0.0733158)

    nr_re = Variable("nr_re", 1.0, 0.01, 0, 0)
    nr_im = Variable("nr_im", 0.0, 0.01, 0, 0)

    phi = Resonances.RBW(
        "phi", phi_re, phi_im, phi_mass, phi_width, 1, PAIR_12, True, True
    )
    f2p = Resonances.RBW(
        "f2", f2p_re, f2p_im, f2p_mass, f2p_width, 2, PAIR_12, True, True
    )
    nr = Resonances.NonRes("nr", nr_re, nr_im)

    B2KKK.resonances = (phi, f2p, nr)

    observables = (s12, s13)
    offsets = (constantZero, constantZero)
    coefficients = (constantOne,)

    vetoDp = veto()

    effwithveto = ProdPdf("effwithveto", [vetoDp, eff])

    d = Amp3Body("signalPDFwithveto", s12, s13, eventNumber, B2KKK, effwithveto)

    return d


# ## Make toy and plotting

# In[11]:


def maketoy(dp):
    print(B2KKK)
    prod = ProdPdf("totalSignal", [dp])
    dplotted = DalitzPlotter(prod, dp)
    toyData = UnbinnedDataSet(s12, s13, eventNumber)
    dplotted.fillDataSetMC(toyData, 10000)
    return toyData


# In[12]:


def plot(toyData, name):
    # plt.figure(0,figsize=(15,5))
    plt.subplot(131)
    plt.hist2d(toyData[0], toyData[1], bins=[100, 100], norm=LogNorm())
    plt.subplot(132)
    plt.hist(toyData[0], bins=100, log=False)
    plt.subplot(133)
    plt.hist(toyData[1], bins=100, log=False)
    plt.savefig(name)
    plt.show()


# ## Filling toy sample

# In[13]:


eff = makePolyEff()
dp = makesignal(eff)
toyData = maketoy(dp)
plt.figure(0, figsize=(15, 5))
plot(toyData, "B2KKK_Simple_toyData.png")


# ## Initial Fit fractions

# In[14]:


print("ff Laura++:")
print("phi = 0.344")
print("f'2 = 0.333")
print("NR = 0.326")
print("\n")

ffs = np.matrix(dp.fit_fractions(True))
print("ff interferences GooFit (%):")
print(ffs.view())


# ## Fitting signal

# In[14]:


prod = ProdPdf("totalSignal", [dp])
prod.setData(toyData)
dp.setDataSize(toyData.getNumEvents())
fitter = FitManager(prod)
fitter.setVerbosity(1)
fitter.setMaxCalls(200000)
print("Running fit...")
func_min = fitter.fit()
print(dp.normalize())


# In[16]:


ffs = np.matrix(dp.fit_fractions(True))
print("ff interferences GooFit (%):")
print(ffs.view())


# In[15]:


toyAfterFit = maketoy(dp)
plt.figure(1, figsize=(15, 5))
plot(toyAfterFit, "B2KKK_Simple_toyData_after.png")


# ## Fitting signal with veto (Not Running)

# In[ ]:


dpwithveto = makesignalwithveto(eff)
prodwithveto = ProdPdf("totalSignal_withVeto", [dpwithveto])


# In[ ]:


prodwithveto.setData(toyData)
dpwithveto.setDataSize(toyData.getNumEvents())
fitter = FitManager(prodwithveto)
fitter.setVerbosity(1)
fitter.setMaxCalls(200000)
print("Running fit...")
func_min = fitter.fit()


# In[ ]:


ffs = np.matrix(dpwithveto.fit_fractions(True))
print("ff interferences GooFit (%):")
print(ffs.view())


# In[ ]:


toyAfterFitWithVeto = maketoy(dpwithveto)
plt.figure(2, figsize=(15, 5))
plot(toyAfterFitWithVeto)


# In[ ]:
