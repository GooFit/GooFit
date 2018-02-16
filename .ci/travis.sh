#!/usr/bin/env bash
set -evx

cd ${TRAVIS_BUILD_DIR}

mkdir build || true
cd build
cmake -DGOOFIT_DEVICE=OMP -DCMAKE_BUILD_TYPE=RelWithDebInfo -DGOOFIT_PYTHON=ON -DGOOFIT_SPLASH=OFF ..
cmake --build . -- -j2

./examples/RunAll.py

export PYTHONPATH=`pwd`

pytest
./pyexamples/RunAll.sh

set +evx
