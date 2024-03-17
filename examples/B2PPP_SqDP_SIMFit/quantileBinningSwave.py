import physt
import matplotlib.pyplot as plt
import uproot
import pandas
import numpy
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline
from math import cos,sin,pi

file = uproot.open("MC/MC_Swave.root:genResults")
s12 = file['s13'].array(library="pd")

bins = 30-1
s12_hist = physt.h1(s12,"quantile",bin_count=bins)

# s12_hist.plot(density=True)

bins_left_edges = s12_hist.bin_left_edges
bins_right_edges= s12_hist.bin_right_edges


Bins = numpy.sqrt(numpy.append(bins_left_edges,bins_right_edges[-1]))

print(Bins, len(Bins))

pwa_coefs = pandas.DataFrame(Bins,columns=['bins'])
pwa_coefs['mags'] = numpy.ones(len(Bins))
pwa_coefs['phases'] = numpy.zeros(len(Bins))

pwa_coefs.to_csv('PWAFile.bin',sep='\t',index=False,header=False, float_format='%.3f')

print(pwa_coefs)

# plt.figure(figsize=(10,5))
# plt.subplot(1,2,1)
# plt.plot(Bins,Mag,'o',x,f_mag(x),'-',s_babar, mag_babar,'x')
# plt.subplot(1,2,2)
# plt.plot(Bins,Phs,'o',x,f_phs(x),'-',s_babar, phs_babar,'x')
# plt.show()
