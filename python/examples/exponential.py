#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

from goofit import *
import numpy as np

print_goofit_info()

xdata = np.random.exponential(size=100000)
xvar = Observable("xvar", 0, np.max(xdata) + 1)

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
