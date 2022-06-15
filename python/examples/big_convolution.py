## This is (essentially) a copy of GooFit/python/examples/convolution.py
## It seems to work, out of the box, converging and giving credible
## central values and errors for fit parameters

## 210930 --  unclear why my very similar code does not work as well
## have change central value of x0 as a test; otherwise is original code

from __future__ import print_function, division

from goofit import *
import numpy as np
import matplotlib.pyplot as plt

print_goofit_info()

import sys
# sys.stderr = open('/share/lazy/sokoloff/D2KKpi/PGun/myErr.convolution', 'w')
##sys.stdout = open('/share/lazy/sokoloff/D2KKpi/PGun/myOut.convolution', 'w')

def cpu_bw(x, x0, gamma):
    ret = gamma
    ret /= 2 * np.sqrt(np.pi)
    ret /= (x - x0) * (x - x0) + 0.25 * gamma * gamma
    return ret

trueX0     = 1.2
trueGamma  = 2.0
trueSigma  = 0.15

xvar  = Observable("xvar", -10, 10)
gamma = Variable("gamma", 1.1*trueGamma, 0.1, 0.1, 5.5)
sigma = Variable("sigma", 0.8*trueSigma, 0.1, 0.1, 5.00)
x0    = Variable("x0", trueX0+0.05, 0.05, -2, 2)
zero  = Variable("zero", 0)

data  = UnbinnedDataSet(xvar)

xvals = np.arange(-10.,10.,0.1)
fvals = np.ones(len(xvals))
for ii in range(len(xvals)):
    fvals[ii] = cpu_bw(xvals[ii],trueX0,trueGamma)
##    print("ii, xvals[ii],fvals[ii] = ",ii, xvals[ii],fvals[ii] )

plt.figure()
plt.plot(xvals,fvals)
plt.show()

raw_bw        = []
random_offset = []
smeared_x     = []

i = -1
while i < 100000:
    i = i + 1

    raw = np.random.uniform(0, 20) - 10
    

    bwvalue = cpu_bw(raw, trueX0, trueGamma)
    roll = np.random.uniform() * (2.0 / (np.sqrt(np.pi) * gamma.value))

    if roll > bwvalue:
        i -= 1
        continue
   
    offset = np.random.normal(loc=0.0, scale=sigma.value)
    xvar.value = raw + offset

    if xvar.value < -10 or xvar.value > 10:
        i -= 1
        continue

    data.addEvent()
    
    raw_bw.append(raw)
    random_offset.append(offset)
    smeared_x.append(xvar.value)
##    if (i<11):
##`        print("i, bwvalue, offset, xvar.value = ", i, bwvalue, offset, xvar.value)
    


breit = BWPdf("breit", xvar, x0, gamma)
gauss = GaussianPdf("gauss", xvar, zero, sigma)


convolution = ConvolutionPdf("convolution", xvar, breit, gauss)
convolution.setData(data)

plt.clf()
plt.figure()
nC, binsC, patches = plt.hist(raw_bw, bins=100, range=[-10.,10.], density=False, facecolor="b", alpha=0.75)
plt.xlabel('xval')
plt.ylabel('frequency (arbitrary units)')
plt.title('raw_bw distribution')
plt.axis([-10.,10.,0,1.15*max(nC)])
plt.grid(True)
plt.savefig('big_convolution_bw.png')

plt.clf()
plt.figure()
nC, binsC, patches = plt.hist(random_offset, bins=100, range=[-10.,10.], density=False, facecolor="b", alpha=0.75)
plt.xlabel('xval')
plt.ylabel('frequency (arbitrary units)')
plt.title('random_offset')
plt.axis([-10.,10.,0,1.15*max(nC)])
plt.grid(True)
plt.savefig('big_convolution_random_offset.png')

plt.clf()
plt.figure()
nC, binsC, patches = plt.hist(smeared_x, bins=100, range=[-10.,10.], density=False, facecolor="b", alpha=0.75)
plt.xlabel('random_offset')
plt.ylabel('frequency (arbitrary units)')
plt.title('smeared_x')
plt.axis([-10.,10.,0,1.15*max(nC)])
plt.grid(True)
plt.savefig('big_convolution_smeared.png')

convolution.setIntegrationConstants(-10,10,0.02)

fitter = FitManager(convolution)
fitter.fit()