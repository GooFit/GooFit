#!/usr/bin/env python
"""
This is the Quad1F example from Minuit2, implemented in Python.
"""

from goofit import minuit2

class Quad1F(minuit2.FCNBase):
    def Up(self):
        return 1.0
    def __call__(self, vect):
        return vect[0]**2

fcn = Quad1F();
upar = minuit2.MnUserParameters()
upar.Add("x", 1., 0.1)
migrad = minuit2.MnMigrad(fcn, upar)

minuit2.MnPrint.SetLevel(3)

min = migrad()

print("min =", min)
