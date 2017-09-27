#!/usr/bin/env python
try:
    from skbuild import setup
except ImportError:
    print("Failed to find scikit-build, please run `pip install scikit-build cmake`")
    raise

setup(
        name='goofit',
        version='2.1.0.beta1',
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
            extras_require={
            'dev': [
                'pytest',
                'pylandau'
            ]
        }
    }

        long_description='''\
GooFit for Python
-----------------

GooFit is a highly parallel fitting framework originally designed for High Energy Physics.

Installation: pipenv
====================

Best method for an install you don't expect to change!
Use pipenv, the officially recommended method for managing installs, virtual environments, and dependencies.

To install::

    pipenv install goofit

To run a shell::

    pipenv shell


Installation: pip
=================

This package can be installed with pip, but uses SciKit-Build, and is build, fully optimized, on your system. Because of this, there are a few caveats when running a pip install. Make sure you have SciKit-Build (``pip install scikit-build``) before you attempt an install. Also, if you don't have a recent version of CMake (3.8 or better recommended), also run ``pip install cmake``. When you build, you should also use pip's ``-v`` flag, so that you can see it build (and observe the
configuration options). Otherwise, you might wait a very long time without output (especially if CUDA was found).

In practice, this looks like this::

    pip install scikit-build cmake
    pip install -v goofit


Installation: local
===================

If you want to add PDFs to GooFit, or use GooFit pacakges, you should be working in a local directory using git. In the following example, I'm assuming you've set up SSH keys with GitHub; you can use https instead if you perfer by changing the URL to ``https://github.com/GooFit/GooFit.git``::

    git clone --recursive git@github.com:GooFit/GooFit.git
    cd goofit
    pipenv install --verbose -e .

And, to use::

    pipenv shell

You can set the ``PIP_INSTALL_OPTIONS`` variable to pass through build command, for example::

    PIP_INSTALL_OPTIONS="-- -DGOOFIT_PACKAGES=OFF" pipenv install --verbose -e .

You can also install development requirements with the ``--dev`` flag.

Building a source package from git
==================================

For developers only:

To make a source package, start with a clean (such as new) git GooFit package with all submodules checked out::

    git clone --branch=master --recursive --depth=10 git@github.com:GooFit/GooFit.git
    cd goofit
    python setup.py sdist
    twine upload dist/*

'''
        )

