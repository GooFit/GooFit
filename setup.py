#!/usr/bin/env python
try:
    from skbuild import setup
except ImportError:
    print("Failed to find scikit-build, please run `pip install scikit-build cmake`)
    raise

setup(
        name='goofit',
        version='2.1.0.dev1',
        description='GooFit fitting package',
        author='Henry Schreiner',
        author_email='hschrein@cern.ch',
        url='https://goofit.github.io',
        platforms = ["POSIX"],
        provides = ["goofit"],
        cmake_args=[
            '-DGOOFIT_PYTHON=ON',
            '-DGOOFIT_EXAMPLES=OFF'],
        license="LGPL 3.0",
        packages=['goofit']
        )

