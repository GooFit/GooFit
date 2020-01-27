#!/usr/bin/env bash

echo -en 'travis_fold:start:script.build\\r'
echo "Building..."
set -evx

mkdir -p build
cd build
cmake -DGOOFIT_DEVICE=OMP -DCMAKE_BUILD_TYPE=RelWithDebInfo -DGOOFIT_PYTHON=ON -DGOOFIT_SPLASH=OFF ..
cmake --build . -- -j2

set +evx
echo -e 'travis_fold:end:script.build\\r'
echo -en 'travis_fold:start:script.pytest\\r'
echo "Python testing..."
set -evx

export PYTHONPATH=`pwd`
pytest

set +evx
echo -e 'travis_fold:end:script.pytest\\r'
echo -en 'travis_fold:start:script.pyexamples\\r'
echo "Python examples..."
set -evx

./pyexamples/RunAll.sh

set +evx
echo -e 'travis_fold:end:script.pyexamples\\r'
