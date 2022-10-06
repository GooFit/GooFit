import numpy as np
import pytest

import goofit as gf


def test_min_conv():
    xvar = gf.Observable("xvar", 1690, 1890)
    xvar.setNumBins(200)
    gf.UnbinnedDataSet(xvar)

    # We want to define a PDF shape corresponding to a single basic Pdfs
    # This doe NOT seem to break the indexing in ConvolutionPdf

    #  mds 211021 -- create a "minimal" version of photosShape using two pure
    #  exponential added using ProdPdf

    expX0_value = 1860.0
    expX0 = gf.Variable("expX0", expX0_value)
    expAlpha_value = 0.4676863526561826
    expAlpha = gf.Variable("expAlpha", expAlpha_value)
    pureExp = gf.ExpPdf("pureExp", xvar, expAlpha, expX0)

    sigma = gf.Variable("sigma", 3.4, 0.1, 0.1, 5)
    zero = gf.Variable("zero", 0.0)
    gauss = gf.GaussianPdf("gauss", xvar, zero, sigma)

    #  here we are creating the simplest ConvolutionPdf
    convolution = gf.ConvolutionPdf("convolution", xvar, pureExp, gauss)

    # a ConvolutionPdf needs an integration range; set this one explicitly
    convolution.setIntegrationConstants(1695, 1885, 0.1)

    grid = convolution.makeGrid()
    grid.to_matrix().flatten()
    convolution.setData(grid)

    # this is the point at which the program crashses;
    # it never gets as far as ' print("convolution computedProbs[0] = " ...
    # convolution.getCompProbsAtDataPoints()
    computedProbs = convolution.getCompProbsAtDataPoints()

    assert len(computedProbs) == 3

    computedProbsArray = np.asarray(computedProbs)
    assert computedProbsArray.shape == (3, 200)

    sums = np.sum(computedProbsArray, axis=1)

    assert sums == pytest.approx(np.array([1, 1, 0.0]))
