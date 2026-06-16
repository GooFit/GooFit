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


GooFit will automatically look for CUDA, and build in GPU mode if it finds CUDA. You can pick a specific version by passing through a CMake option (see below), or by setting an environment variable, ``GOOFIT_DEVICE`` before building. You may want to build with OpenMP as a backend to avoid using your GPU, or you might want the CPP version if you are using Anaconda on macOS. Here are the three common backends::

    pip install -v goofit --config-settings=cmake.define.GOOFIT_DEVICE=CUDA
    pip install -v goofit --config-settings=cmake.define.GOOFIT_DEVICE=OMP
    pip install -v goofit --config-settings=cmake.define.GOOFIT_DEVICE=CPP

The lines above use environment variables; GooFit will find any environment variables that start with ``GOOFIT_*`` and set them as CMake defines. If you want to send arbitrary commands to CMake through PIP, you will need to pass each option through, starting with a ``--`` option. Pip will try to reuse the built version if you do not pass options, but will rebuild if you pass options, so this works for a rebuild, unlike the lines above. This is how you would do this to set OMP as the backend::

    pip install -v goofit --config-settings=cmake.define.GOOFIT_DEVICE=OMP


Installation: local
===================

If you want to add PDFs to GooFit, or use GooFit packages, you should be working in a local directory using git. In the following example, I'm assuming you've set up SSH keys with GitHub; you can use https instead if you prefer by changing the URL to ``https://github.com/GooFit/GooFit.git``::

    git clone --recurse-submodules git@github.com:GooFit/GooFit.git
    cd goofit

Local pip
~~~~~~~~~

The normal install here works, though as usual you should include verbose output and you should be in a virtual environment (standard practice)::

    pip install -v .
