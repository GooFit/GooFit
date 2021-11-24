from goofit import *
import numpy as np
import ROOT
import pandas as pd
import matplotlib.pyplot as plt
import joblib 
import uproot

#from hep_ml.reweight import  BinsReweighter, GBReweighter,FoldingReweighter
#from hep_ml import reweight

print_goofit_info()

_mD0       = 1.8645
piPlusMass = 0.13957018
KmMass     = .493677

columns = ['m12','m34','c12','c34','phi','dtime']

#plotting labels
m12label = r'm$_{\pi^{\mp}\pi^{\pm}}$ [GeV/c$^{2}$]'
m34label = r'm$_{K^{\mp}\pi^{\pm}}$ [GeV/c$^{2}$]'
c12label = r'$cos(\theta_{1})$ [rad]'
c34label = r'$cos(\theta_{2})$ [rad]'
philabel = r'$\phi$ [rad]' 
dTimeLabel = r'$D^{0}$ Decay Time [ns]'

#read in generated MC data that we will use to test the fit
mc_df = pd.read_csv('amp4body_cpp_generated_events.txt',names=columns,delimiter=" ")
pdf_df = pd.read_csv('amp4body_example_pdf_values.txt',names=['total_pdf','sig_pdf','bkg_pdf'],delimiter=" ")
mc_df['eff_weight'] = np.ones(len(mc_df))
print(f"number of phase space events: {len(mc_df)} number of pdf values: {len(pdf_df)}")
mc_df['sig_pdf_val'] = pdf_df['sig_pdf']
mc_df['total_pdf_val'] = pdf_df['total_pdf']
total_pdf_sum = sum(mc_df['total_pdf_val'])
print(total_pdf_sum)
#print(mc_df.head())
"""
m12_n,bins = np.histogram(mc_df['m12']*mc_df['total_pdf_val'],bins=100)
m12_n = (m12_n/total_pdf_sum)*len(mc_df)
plt.bar(bins[:-1],m12_n,width=np.diff(bins))
"""
#_,bins,_ = plt.hist((mc_df['m12']*mc_df['total_pdf_val'])/total_pdf_sum,density=True,bins=100,histtype='step',label='Signal PDF weighted')
#plt.hist(mc_df['m12'],density=False,bins=bins,label='MC',histtype='step')
#plt.legend(loc='best')
#_,bins,_ = plt.hist(mc_df['m12'],weights=mc_df['total_pdf_val'],density=True,label='Signal PDF',histtype='step',bins=100)
m12_n,bins = np.histogram(mc_df['m12'],bins=100)
plt.hist(mc_df['m12'],density=False,bins=bins,label='MC',histtype='step')
plt.legend(loc='best')
plt.xlabel(m12label)
bin_width = abs(bins[1]-bins[0])/len(bins)
plt.ylabel(fr'Entries/{{{bin_width:.3e}}} GeV/c$^{2}$')
plt.savefig('plots/signal_pdf_component.png')