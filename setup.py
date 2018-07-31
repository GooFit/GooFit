#!/usr/bin/env python
try:
    from skbuild import setup
except ImportError:
    print("Failed to find scikit-build, please run `pip install scikit-build`")
    raise

setup(
        name='goofit',
        version='2.2.0',
        description='GooFit fitting package',
        author='Henry Schreiner',
        author_email='hschrein@cern.ch',
        url='https://goofit.github.io',
        platforms = ["POSIX"],
        provides = ["goofit"],
        install_requires = [
            'numpy>=1.11.1',
        ],
        classifiers = [
            "Development Status :: 5 - Production/Stable",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
            "Natural Language :: English",
            "Operating System :: Unix",
            "Programming Language :: C++",
            "Programming Language :: Python",
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 2.7",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.4",
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Topic :: Scientific/Engineering :: Physics"
        ],
        cmake_args=[
            '-DGOOFIT_PYTHON=ON',
            '-DGOOFIT_CERNROOT=OFF',
            '-DGOOFIT_EXAMPLES=OFF'],
        license="LGPL 3.0",
        packages=['goofit'],
        extras_require={
            'dev': [
                'pytest',
                'matplotlib>=1.5',
                'pandas>=0.15.1',
                'uncertainties>=3.0.2',
                'scipy'
            ]
        },
        long_description='''\
GooFit for Python
-----------------

GooFit is a highly parallel fitting framework originally designed for High Energy Physics.

Installation basics
===================

This package can be installed with pip, but uses SciKit-Build, and is build, fully optimized, on your system. Because of this, there are a few caveats when running a pip install. Make sure you have SciKit-Build (``pip install scikit-build``) before you attempt an install. Also, if you don't have a recent version of CMake (3.8 or better recommended), also run ``pip install cmake``. When you build, you should also use pip's ``-v`` flag, so that you can see it build (and observe the
configuration options). Otherwise, you might wait a very long time without output (especially if CUDA was found).


Installation: pip
=================

Using pip 10+::

    pip install -v goofit

Using pip < 10`::

    pip install scikit-build
    pip install -v goofit


If you want to send commands to CMake through PIP, use (for example)::

    PIP_INSTALL_OPTIONS="-- -DGOOFIT_PACKAGES=OFF" pip install -v goofit


Installation: local
===================

If you want to add PDFs to GooFit, or use GooFit pacakges, you should be working in a local directory using git. In the following example, I'm assuming you've set up SSH keys with GitHub; you can use https instead if you perfer by changing the URL to ``https://github.com/GooFit/GooFit.git``::

    git clone --recursive git@github.com:GooFit/GooFit.git
    cd goofit

Pipenv
~~~~~~

You can set up a quick environment using pipenv::

    pipenv install --dev

Then activate that environment::

    pipenv shell

Local pip
~~~~~~~~~

The normal install here works, though as usual you should include verbose output::

    pip install -v .


You can set the ``PIP_INSTALL_OPTIONS`` variable to pass through build command, for example::

    PIP_INSTALL_OPTIONS="-- -DGOOFIT_PACKAGES=OFF" pip install -v .


Building a source package from git
==================================

For developers only:

To make a source package, start with a clean (such as new) git GooFit package with all submodules checked out::

    git clone --branch=master --recursive --depth=10 git@github.com:GooFit/GooFit.git
    cd goofit
    python setup.py sdist
    python -m twine upload dist/*

To make a binary package, use instead::

    python setup.py bdist_wheel -- -DGOOFIT_OPTI="-march=core2"

'''
        )

# To set this up on Docker for linux, use::
#
#    docker run -it quay.io/pypa/manylinux1_x86_64 -v goofit-py:goofit-py
#    export PATH=/opt/python/cp36-cp36m/bin:$PATH
#    cd goofit-py
#    python -m pip install scikit-build cmake
#    python setup.py bdist_wheel -- -DGOOFIT_OPTI="-march=core2"

