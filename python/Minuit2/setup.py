from pybind11.setup_helpers import Pybind11Extension
from setuptools import setup

__version__ = "6.24.00"


minuit2_pybind_sources = [
    "FCNBase.cpp",
    "MinuitParameter.cpp",
    "MnUserParameters.cpp",
    "MnUserParameterState.cpp",
    "MnUserCovariance.cpp",
    "FunctionMinimum.cpp",
    "MnApplication.cpp",
    "MnMigrad.cpp",
    "MnPrint.cpp",
    "main.cpp",
]

math_sources = [
    "GenAlgoOptions.cxx",
    "MinimizerOptions.cxx",
]

minuit2_sources = [
    "AnalyticalGradientCalculator.cxx",
    "BasicMinimumError.cxx",
    "CombinedMinimumBuilder.cxx",
    "DavidonErrorUpdator.cxx",
    "FumiliBuilder.cxx",
    "FumiliErrorUpdator.cxx",
    "FumiliGradientCalculator.cxx",
    "FumiliMinimizer.cxx",
    "FumiliStandardChi2FCN.cxx",
    "FumiliStandardMaximumLikelihoodFCN.cxx",
    "HessianGradientCalculator.cxx",
    "InitialGradientCalculator.cxx",
    "LaEigenValues.cxx",
    "LaInnerProduct.cxx",
    "LaInverse.cxx",
    "LaOuterProduct.cxx",
    "LaSumOfElements.cxx",
    "LaVtMVSimilarity.cxx",
    "MPIProcess.cxx",
    "MinimumBuilder.cxx",
    "Minuit2Minimizer.cxx",
    "MnApplication.cxx",
    "MnContours.cxx",
    "MnCovarianceSqueeze.cxx",
    "MnEigen.cxx",
    "MnFcn.cxx",
    "MnFumiliMinimize.cxx",
    "MnFunctionCross.cxx",
    "MnGlobalCorrelationCoeff.cxx",
    "MnHesse.cxx",
    "MnLineSearch.cxx",
    "MnMachinePrecision.cxx",
    "MnMinos.cxx",
    "MnParabolaFactory.cxx",
    "MnParameterScan.cxx",
    "MnPlot.cxx",
    "MnPosDef.cxx",
    "MnPrint.cxx",
    "MnPrintImpl.cxx",
    "MnScan.cxx",
    "MnSeedGenerator.cxx",
    "MnStrategy.cxx",
    "MnTiny.cxx",
    "MnTraceObject.cxx",
    "MnUserFcn.cxx",
    "MnUserParameterState.cxx",
    "MnUserParameters.cxx",
    "MnUserTransformation.cxx",
    "ModularFunctionMinimizer.cxx",
    "NegativeG2LineSearch.cxx",
    "Numerical2PGradientCalculator.cxx",
    "ParametricFunction.cxx",
    "ScanBuilder.cxx",
    "SimplexBuilder.cxx",
    "SimplexParameters.cxx",
    "SimplexSeedGenerator.cxx",
    "SinParameterTransformation.cxx",
    "SqrtLowParameterTransformation.cxx",
    "SqrtUpParameterTransformation.cxx",
    "VariableMetricBuilder.cxx",
    "VariableMetricEDMEstimator.cxx",
    "mnbins.cxx",
    "mndasum.cxx",
    "mndaxpy.cxx",
    "mnddot.cxx",
    "mndscal.cxx",
    "mndspmv.cxx",
    "mndspr.cxx",
    "mnlsame.cxx",
    "mnteigen.cxx",
    "mntplot.cxx",
    "mnvert.cxx",
    "mnxerbla.cxx",
]


# Add the correct directories
math_sources = ["src/Math/" + m for m in math_sources]
minuit2_sources = ["src/Minuit2/" + m for m in minuit2_sources]
minuit2_pybind_sources = ["python/" + m for m in minuit2_pybind_sources]

ext_modules = [
    Pybind11Extension(
        "minuit2",
        math_sources + minuit2_sources + minuit2_pybind_sources,
        include_dirs=[
            "inc",
        ],
        language="c++",
        std_cxx="11",
        define_macros=[("WARNINGMSG", None), ("MATH_NO_PLUGIN_MANAGER", None)],
    ),
]


setup(
    name="minuit2",
    version=__version__,
    author="Henry Schriener",
    author_email="schreihf@cern.ch",
    url="https://github.com/GooFit/Minuit2",
    description="A Pybind11 Minuit2 binding",
    long_description="",
    ext_modules=ext_modules,
    install_requires=["numpy>=1.10"],
    zip_safe=False,
)
