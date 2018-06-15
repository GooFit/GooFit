#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

from goofit import *
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

print_goofit_info()

xvar = Observable("xvar", -5, 5)
data = UnbinnedDataSet(xvar)

# Classic C++ style method:
#i = 0
#while i < 100000:
#    xvar.value = np.random.normal(0.2,1.1)
#
#    if np.random.uniform() < 0.1:
#        xvar.value = np.random.uniform(-5,5)
#
#    if abs(xvar.value) <= 5:
#        data.addEvent()
#        i += 1
# dat = data.to_matrix().flatten()

# Better method in Python:
dat = np.random.normal(0.2,1.1,100000)
dat[:10000] = np.random.uniform(-5,5,10000)
data.from_matrix([dat], filter=True)


xmean = Variable("xmean", 0, 1, -10, 10)
xsigm = Variable("xsigm", 1, 0.5, 1.5)
signal = GaussianPdf("signal", xvar, xmean, xsigm)

constant = Variable("constant", 1.0)
backgr = PolynomialPdf("backgr", xvar, [constant])

sigfrac = Variable("sigFrac", 0.9, 0.75, 1.00)

total = AddPdf("total", [sigfrac], [signal, backgr])
total.fitTo(data)

# Plot data
print(dat, dat.shape)
plt.hist(dat, bins='auto', label='data', density=True)

# Make grid and evaluate on it
grid = total.makeGrid()
total.setData(grid)
main, gauss, flat = total.getCompProbsAtDataPoints()
xvals = grid.to_matrix().flatten()

plt.plot(xvals, main, label='total')
plt.plot(xvals, np.array(gauss)*sigfrac.value, label='signal')
plt.plot(xvals, np.array(flat)*(1-sigfrac.value), label='background')

plt.legend()
plt.savefig('addition_plot.pdf')
#plt.show()

