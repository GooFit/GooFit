#!/usr/bin/env python
try:
    from skbuild import setup
except ImportError:
    print("Failed to find scikit-build, please run `pip install scikit-build cmake`")
    raise

setup(
        name='goofit',
        version='2.1.0.dev2',
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
        packages=['goofit'],
        long_description='''\
GooFit for Python
-----------------

GooFit is a highly parallel fitting framework originally designed for High Energy Physics.

Installation
============

This package can be installed with pip, but uses SciKit-Build, and is build, fully optimized, on your system. Because of this, there are a few caveats when running a pip install. Make sure you have SciKit-Build (``pip install scikit-build``) before you attempt an install. Also, if you don't have a recent version of CMake (3.8 or better recommended), also run ``pip install cmake``. When you build, you should also use pip's ``-v`` flag, so that you can see it build (and observe the
configuration options). Otherwise, you might wait a very long time without output (especially if CUDA was found).

In practice, this looks like this::

    pip install scikit-build cmake
    pip install -v goofit


Building a source package from git
==================================

For developers:

To make a source package, start with a clean (such as new) git GooFit package with all submodules checked out::

    git clone --branch=master --recursive --depth=10 git@github.com:GooFit/GooFit.git
    cd goofit
    python setup.py sdist
    twine upload dist/*

'''
        )

