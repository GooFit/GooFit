#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

from goofit import *
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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

grid = total.makeGrid()
total.setData(grid)

val = total.getCompProbsAtDataPoints()

mat = np.array(val[0]).reshape([xvar.numbins, yvar.numbins])

fig, axs = plt.subplots(1, 2, figsize=(10,5))

axs[0].hist2d(xyarr[0], xyarr[1], bins=(100,100), range=((xvar.lowerlimit, xvar.upperlimit), (yvar.lowerlimit, yvar.upperlimit)))
axs[1].imshow(mat, extent=(xvar.lowerlimit, xvar.upperlimit, yvar.lowerlimit, yvar.upperlimit), origin='lower')

for ax in axs:
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.axis("equal")

axs[0].set_title("Original input (histogram)")
axs[1].set_title("Fitted function")

fig.savefig("plot_2d.png")

