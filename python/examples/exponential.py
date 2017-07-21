#!/usr/bin/env python

from goofit import *
import numpy as np

xdata = np.random.exponential(size=100000)
xvar = Variable("xvar", 0, np.max(xdata) + 1)
data = UnbinnedDataSet(xvar)

for v in xdata:
    xvar.value = v
    data.addEvent()

alpha = Variable("alpha", -2, 0.1, -10, 10)
exppdf = ExpPdf("exppdf", xvar, alpha)

exppdf.fitTo(data)

#exppdf.setData(data)
#fitter = FitManager(exppdf)
#fitter.fit()
