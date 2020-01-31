#!/usr/bin/env python
try:
    from skbuild import setup
except ImportError:
    print("Failed to find scikit-build, please run `pip install scikit-build`")
    raise

import os

ITEMS = [
    '-DGOOFIT_PYTHON=ON',
    '-DGOOFIT_TESTS=OFF',
    '-DGOOFIT_CERNROOT=OFF',
    '-DGOOFIT_EXAMPLES=OFF',
    '-DCMAKE_UNITY_BUILD=ON', # Faster build on CMake 3.16+
]

# Add GOOFIT_* from env.
for item in os.environ:
    if item.startswith('GOOFIT_'):
        ITEMS.append('-D{0}={1}'.format(item, os.environ[item]))



setup(
        name='goofit',
        version='2.2.3',
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
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Topic :: Scientific/Engineering :: Physics"
        ],
        cmake_args=ITEMS,
        license="LGPL 3.0",
        packages=['goofit'],
        extras_require={
            'dev': [
                'pytest',
                'matplotlib>=1.5',
                'pandas>=0.15.1',
                'uncertainties>=3.0.2',
                'scipy',
                'plumbum'
            ]
        },
        long_description='''\
GooFit for Python
-----------------

GooFit is a highly parallel fitting framework originally designed for High Energy Physics.

Installation basics
===================

This package can be installed with pip, but uses SciKit-Build, and is build,
fully optimized, on your system. Because of this, there are a few caveats when
running a pip install if you use an old version of pip. When you build, you
should also use pip's ``-v`` flag, so that you can see it build (and observe
the configuration options). Otherwise, you might wait a very long time without
output (especially if CUDA was found).


Installation: pip
=================

Using pip 10+::

    pip install -v goofit

Using pip < 10::

    pip install scikit-build # optionally cmake ninja
    pip install -v goofit


GooFit will automatically look for CUDA, and build in GPU mode if it finds CUDA. You can pick a specific version by passing through a CMake option (see below), or by setting an environment variable, `GOOFIT_DEVICE` before building. You may want to build with OpenMP as a backend to avoid using your GPU, or you might want the CPP version if you are using Anaconda on macOS. Here are the three common backends::

    GOOFIT_DEVICE=CUDA pip install -v goofit
    GOOFIT_DEVICE=OMP pip install -v goofit
    GOOFIT_DEVICE=CPP pip install -v goofit

The lines above use environment variables; GooFit will find any environment variables that start with ``GOOFIT_*`` and set them as CMake defines. If you want to send arbitrary commands to CMake through PIP, you will need to pass each option through, starting with a ``--`` option. Pip will try to reuse the built version if you do not pass options, but will rebuild if you pass options, so this works for a rebuild, unlike the lines above. This is how you would do this to set OMP as the backend::

    pip install -v goofit --install-option="--" --install-option="-DGOOFIT_DEVICE=OMP"
    # OR
    PIP_INSTALL_OPTION="-- -DGOOFIT_DEVICE=OMP" pip install -v goofit


Installation: local
===================

If you want to add PDFs to GooFit, or use GooFit pacakges, you should be working in a local directory using git. In the following example, I'm assuming you've set up SSH keys with GitHub; you can use https instead if you prefer by changing the URL to ``https://github.com/GooFit/GooFit.git``::

    git clone --recursive git@github.com:GooFit/GooFit.git
    cd goofit

Local pip
~~~~~~~~~

The normal install here works, though as usual you should include verbose output and you should be in a virtual environment (standard practice)::

    pip install -v .

'''
        )



# Building a source package from git
# ==================================
# 
# For developers only:
# 
# To make a source package, start with a clean (such as new) git GooFit package with all submodules checked out::
# 
#     git clone --branch=master --recursive --depth=10 git@github.com:GooFit/GooFit.git
#     cd goofit
#     python setup.py sdist
#     python -m twine upload dist/*
# 
# To make a binary package, use instead::
# 
#     GOOFIT_OPTI="" python setup.py bdist_wheel
# 
# To set this up on Docker for linux, use::
#
#    docker run -it quay.io/pypa/manylinux1_x86_64 -v goofit-py:goofit-py
#    export PATH=/opt/python/cp36-cp36m/bin:$PATH
#    cd goofit-py
#    python -m pip install scikit-build cmake
#    python setup.py bdist_wheel -- -DGOOFIT_OPTI="-march=core2"

