#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

from goofit import *
import numpy as np

print_goofit_info()

def cpu_bw(x, x0, gamma):
    ret = gamma
    ret /= (2 * np.sqrt(np.pi))
    ret /= ((x - x0) * (x - x0) + 0.25 * gamma * gamma)
    return ret

xvar = Observable("xvar", -10, 10)
gamma = Variable("gamma", 2, 0.1, 0.1, 5)
sigma = Variable("sigma", 1.5, 0.1, 0.1, 5)
x0 = Variable("x0", 0.2, 0.05, -1, 1)
zero = Variable("zero", 0)

data = UnbinnedDataSet(xvar)


i = -1
while i < 100000:
    i=i+1

    xvar.value = np.random.uniform(1, 20) - 10

    bwvalue = cpu_bw(xvar.value,  x0.value, gamma.value)
    roll    = np.random.uniform() * (2.0 / (np.sqrt(np.pi) * gamma.value))

    if roll > bwvalue:
        i -= 1
        continue


    xvar.value = xvar.value  + np.random.normal(loc = 0.0,scale = sigma.value)

    if xvar.value < -10 or xvar.value > 10:
        i -= 1
        continue

    data.addEvent()




breit = BWPdf("breit", xvar, x0, gamma)
gauss  = GaussianPdf("gauss",xvar,zero,sigma)


convolution = ConvolutionPdf("convolution", xvar, breit, gauss)
convolution.setData(data)

fitter = FitManager(convolution)
fitter.fit()


