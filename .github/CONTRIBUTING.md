## Help setting up

See [`docs/SYSTEM_INSTALL.md`](../docs/SYSTEM_INSTALL.md) for help building and installing.

## Development

All formatting is done with [pre-commit][]. Install with brew (`brew install
pre-commit`) (macOS) or pipx/pip (any OS).

[pre-commit]: https://pre-commit.com


## Making a release (for developers)

* Make version bump PR
  * Bump the version numbers in setup.py and CMakeLists.txt
  * Set the `GOOFIT_TAG` to `"RELEASE"` in CMakeLists.txt
* Make GitHub release (web interface)
  * Github page -> Releases -> Draft a new release
  * Copy the release notes (should be up to date!) to the release (include link definitions)
  * CI should automatically append source package and upload docs / PyPI packages for you
* Return `GOOFIT_TAG` to `"dev"`

## Building a source package from git

For developers only:

To make a source package, start with a clean (such as new) git GooFit package with all submodules checked out:

    git clone --branch=master --recursive --depth=10 git@github.com:GooFit/GooFit.git
    cd goofit
    pipx run build --sdist
    pipx run twine upload dist/*

To make a binary package, use instead:

    pipx run build --wheel -Ccmake.define.GOOFIT_OPTI=""

To set this up on Docker for linux, use:

   docker run -it quay.io/pypa/manylinux1_x86_64 -v goofit-py:goofit-py
   export PATH=/opt/python/cp39-cp39m/bin:$PATH
   cd goofit-py
   python -m pip install -v . -Ccmake.define.GOOFIT_OPTI="-march=core2"
