# Only needed if run from plain CMake build
# If installed with pip, this is not needed
import sys
sys.path.append('.')

from goofit import *

import numpy as np

def test_exp():
    xdata = np.random.exponential(size=100000)
    xvar = Variable("xvar", 0, np.max(xdata) + 1)
    data = UnbinnedDataSet(xvar)

    for v in xdata:
        xvar.value = v
        data.addEvent()

    alpha = Variable("alpha", -2, 0.1, -10, 10)
    exppdf = ExpPdf("exppdf", xvar, alpha)
    exppdf.setData(data)

    fitter = FitManager(exppdf)
    fitter.fit()

    assert(abs(alpha.value + 1) < .01)
    assert(alpha.error < .05)

def test_exp_simple():
    xdata = np.random.exponential(size=100000)
    xvar = Variable("xvar", 0, np.max(xdata) + 1)
    data = UnbinnedDataSet(xvar)
    data.from_numpy(xdata.reshape(1,-1))

    alpha = Variable("alpha", -2, 0.1, -10, 10)
    exppdf = ExpPdf("exppdf", xvar, alpha)

    exppdf.fitTo(data)

    assert(abs(alpha.value + 1) < .01)
    assert(alpha.error < .05)

def test_exp_eigen():
    xdata = np.random.exponential(size=100000).reshape([1,-1])
    xvar = Variable("xvar", 0, np.max(xdata) + 1)
    data = UnbinnedDataSet(xvar)
    data.from_matrix(xdata, False)
    new_mat = data.to_matrix()
    new_mat2 = data.to_numpy()

    np.testing.assert_array_equal(xdata, new_mat)
    np.testing.assert_array_equal(xdata, new_mat2)

    alpha = Variable("alpha", -2, 0.1, -10, 10)
    exppdf = ExpPdf("exppdf", xvar, alpha)

    exppdf.fitTo(data)

    assert(abs(alpha.value + 1) < .01)
    assert(alpha.error < .05)
