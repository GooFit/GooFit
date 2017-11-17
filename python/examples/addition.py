#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

from goofit import *
import numpy as np

print_goofit_info()

xvar = Observable("xvar", -5, 5)

data = UnbinnedDataSet(xvar)

totalData = 0


i=-1

while i < 100000:
    i+=1
    xvar.value = np.random.normal(0.2,1.1)

    if np.random.uniform() < 0.1:
        xvar.value = np.random.uniform(-5,5)


    if abs(xvar.value) > 5:
        i-=1
        continue

    data.addEvent()


xmean = Variable("xmean", 0, 1, -10, 10)
xsigm = Variable("xsigm", 1, 0.5, 1.5)
signal = GaussianPdf("signal", xvar, xmean, xsigm)


constant = [Variable("constant", 1.0)]
backgr = PolynomialPdf("backgr", xvar, constant)

sigFrac = Variable("sigFrac", 0.9, 0.75, 1.00)

total = AddPdf("total", [sigFrac], [signal,backgr])
total.setData(data)

fitter = FitManager(total)
fitter.fit()
