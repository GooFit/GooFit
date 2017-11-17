#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

from goofit import *
import numpy as np

print_goofit_info()

# Create random data
xarr = np.random.normal(.2, 1.1, size=100000)
yarr = np.random.normal(.5, .3, size=100000)
xyarr = np.array([xarr, yarr])

xvar = Observable("xvar", -5, 5)
yvar = Observable("yvar", -5, 5)

data = UnbinnedDataSet(xvar, yvar)

data.from_numpy(xyarr, filter=True)

#for x,y in zip(xarr, yarr):
#    xvar.value = x
#    yvar.value = y
#    data.addEvent()

xmean = Variable("xmean", 0, 1, -10, 10)
xsigm = Variable("xsigm", 1, 0.5, 1.5)
xgauss = GaussianPdf("xgauss", xvar, xmean, xsigm)

ymean = Variable("ymean", 0, 1, -10, 10)
ysigm = Variable("ysigm", 0.4, 0.1, 0.6)
ygauss = GaussianPdf("ygauss", yvar, ymean, ysigm)

total = ProdPdf("total", [xgauss, ygauss])
total.setData(data)

fitter = FitManager(total)
fitter.fit()

print(xmean)
print(xsigm)
print(ymean)
print(ysigm)

print(data.to_numpy())
grid = total.makeGrid()
print(grid.to_numpy())
total.setData(grid)
