# -*- coding: utf-8 -*-
# modified from December 30, 2021
# the goal is to produce a ConvolutionPdf from a model with two components and a resolution function that is "basic".
# the existing version of GooFit should fail to calculate probabilities at data points; with luck
# a to-be-corrected version will succeed, and then ConvolutionPdfs will work "recursively" in all cases.
#
# -- this subminimal version of the code does not have the StepPdf modification noted here.
# uses a kernel gFit-22Sept2021 that was built 210922 (no surprise) using my 22sept2021 conda environment
# plus a local version of GooFit that includes an extra argument for StepPdf to indicate whether the
# function steps "up" (normal behavior) or "down", a feature to be used here.
#
# to make this work, the following commands should be executed in the shell before starting Jupyter
# <ul>
#     <li> conda activate 22sept2021 </li>
#     <li>  export PYTHONPATH=/data/sokoloff/GooFit/CPP_211021:/opt/root-6.18.00/lib </li>
# </ul>
#
# Note that this conda environment uses uproot4.  As a result, the code  uses the packakge **vector** rather than **TVectorLorentzArray** to handle four-vectors
# --


import goofit as gf


def test_comp_advanced():
    xvar = gf.Observable("xvar", 1690, 1890)
    xvar.setNumBins(200)
    gf.UnbinnedDataSet(xvar)

    # We want to define a PDF shape corresponding the fit values found above, except for the peakMass,
    # and then use this is a convolution to find the peakMass and the parameters of the resolution function.
    #
    # To do this, we need to set the parameters of the exponentials to those found above.
    #
    # Also, as we will be fitting the smeared distribution rather than the raw distribution, we
    # need to specify that the data set to be fit is "xvar"

    expX0_A_value = 1860.0
    expX0_A = gf.Variable("expX0_A", expX0_A_value)
    expAlpha_A_value = 0.4676863526561826
    expAlpha_A = gf.Variable("expAlpha_A", expAlpha_A_value)
    Exp_A = gf.ExpPdf("Exp_A", xvar, expAlpha_A, expX0_A)

    expX0_B_value = 1860.0
    expX0_B = gf.Variable("expX0_B", expX0_B_value)
    expAlpha_B_value = 0.02429336729138952
    expAlpha_B = gf.Variable("expAlpha_B", expAlpha_B_value)
    Exp_B = gf.ExpPdf("Exp_B", xvar, expAlpha_B, expX0_B)

    fracA = gf.Variable("fracA", 0.5)
    fracB = gf.Variable("fracB", 0.5)
    AddedPdf = gf.AddPdf("photosShape", [fracA, fracB], [Exp_A, Exp_B])

    grid = AddedPdf.makeGrid()
    AddedPdf.setData(grid)
    grid.to_matrix().flatten()

    computedProbs = AddedPdf.getCompProbsAtDataPoints()

    sigma = gf.Variable("sigma", 3.4, 0.1, 0.1, 5)
    zero = gf.Variable("zero", 0.0)

    gauss = gf.GaussianPdf("gauss", xvar, zero, sigma)

    convolution = gf.ConvolutionPdf("convolution", xvar, AddedPdf, gauss)

    grid = convolution.makeGrid()

    convolution.setIntegrationConstants(1695, 1885, 0.1)

    convolution.setData(grid)

    grid.to_matrix().flatten()

    convolution.getCompProbsAtDataPoints()
    computedProbs = convolution.getCompProbsAtDataPoints()

    assert len(computedProbs) == 3
    assert len(computedProbs[0]) == 200
