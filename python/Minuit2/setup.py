from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import setuptools

__version__ = '6.08.06'


class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)

minuit2_pybind_sources = [
    'FCNBase.cpp',
    'MinuitParameter.cpp',
    'MnUserParameters.cpp',
    'MnUserParameterState.cpp',
    'MnUserCovariance.cpp',
    'FunctionMinimum.cpp',
    'MnApplication.cpp',
    'MnMigrad.cpp',
    'MnPrint.cpp',
    'main.cpp',
]

math_sources = [
    'GenAlgoOptions.cxx',
    'MinimizerOptions.cxx',
]

minuit2_sources = [
    'AnalyticalGradientCalculator.cxx',
    'BasicMinimumError.cxx',
    'CombinedMinimumBuilder.cxx',
    'DavidonErrorUpdator.cxx',
    'FumiliBuilder.cxx',
    'FumiliErrorUpdator.cxx',
    'FumiliGradientCalculator.cxx',
    'FumiliMinimizer.cxx',
    'FumiliStandardChi2FCN.cxx',
    'FumiliStandardMaximumLikelihoodFCN.cxx',
    'HessianGradientCalculator.cxx',
    'InitialGradientCalculator.cxx',
    'LaEigenValues.cxx',
    'LaInnerProduct.cxx',
    'LaInverse.cxx',
    'LaOuterProduct.cxx',
    'LaSumOfElements.cxx',
    'LaVtMVSimilarity.cxx',
    'MPIProcess.cxx',
    'MinimumBuilder.cxx',
    'Minuit2Minimizer.cxx',
    'MnApplication.cxx',
    'MnContours.cxx',
    'MnCovarianceSqueeze.cxx',
    'MnEigen.cxx',
    'MnFcn.cxx',
    'MnFumiliMinimize.cxx',
    'MnFunctionCross.cxx',
    'MnGlobalCorrelationCoeff.cxx',
    'MnHesse.cxx',
    'MnLineSearch.cxx',
    'MnMachinePrecision.cxx',
    'MnMinos.cxx',
    'MnParabolaFactory.cxx',
    'MnParameterScan.cxx',
    'MnPlot.cxx',
    'MnPosDef.cxx',
    'MnPrint.cxx',
    'MnScan.cxx',
    'MnSeedGenerator.cxx',
    'MnStrategy.cxx',
    'MnTiny.cxx',
    'MnTraceObject.cxx',
    'MnUserFcn.cxx',
    'MnUserParameterState.cxx',
    'MnUserParameters.cxx',
    'MnUserTransformation.cxx',
    'ModularFunctionMinimizer.cxx',
    'NegativeG2LineSearch.cxx',
    'Numerical2PGradientCalculator.cxx',
    'ParametricFunction.cxx',
    'ScanBuilder.cxx',
    'SimplexBuilder.cxx',
    'SimplexParameters.cxx',
    'SimplexSeedGenerator.cxx',
    'SinParameterTransformation.cxx',
    'SqrtLowParameterTransformation.cxx',
    'SqrtUpParameterTransformation.cxx',
    'VariableMetricBuilder.cxx',
    'VariableMetricEDMEstimator.cxx',
    'mnbins.cxx',
    'mndasum.cxx',
    'mndaxpy.cxx',
    'mnddot.cxx',
    'mndscal.cxx',
    'mndspmv.cxx',
    'mndspr.cxx',
    'mnlsame.cxx',
    'mnteigen.cxx',
    'mntplot.cxx',
    'mnvert.cxx',
    'mnxerbla.cxx',
]


# Add the correct directories
math_sources = ['src/Math/' + m for m in math_sources]
minuit2_sources = ['src/Minuit2/' + m for m in minuit2_sources]
minuit2_pybind_sources = ['python/' + m for m in minuit2_pybind_sources]

ext_modules = [
    Extension(
        'minuit2',
        math_sources + minuit2_sources + minuit2_pybind_sources ,
        include_dirs=[
            # Path to pybind11 headers
            get_pybind_include(),
            get_pybind_include(user=True),
            'inc',
        ],
        language='c++',
        define_macros=[('WARNINGMSG', None),
                       ('MATH_NO_PLUGIN_MANAGER', None)],
    ),
]


# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14] compiler flag.
    The c++14 is preferred over c++11 (when it is available).
    """
    if has_flag(compiler, '-std=c++14'):
        return '-std=c++14'
    elif has_flag(compiler, '-std=c++11'):
        return '-std=c++11'
    else:
        raise RuntimeError('Unsupported compiler -- at least C++11 support '
                           'is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }

    if sys.platform == 'darwin':
        c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args = opts
        build_ext.build_extensions(self)

setup(
    name='minuit2',
    version=__version__,
    author='Henry Schriener',
    author_email='schreihf@cern.ch',
    url='https://github.com/GooFit/Minuit2',
    description='A Pybind11 Minuit2 binding',
    long_description='',
    ext_modules=ext_modules,
    install_requires=['pybind11>=2.2', 'numpy>=1.10'],
    cmdclass={'build_ext': BuildExt},
    zip_safe=False,
)
