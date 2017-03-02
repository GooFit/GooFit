#!/usr/bin/env sh
set -evx
env | sort

mkdir build || true
cd build
cmake -DGOOFIT_DEVICE=OMP -DGOOFIT_SEPARATE_COMP=ON -DGOOFIT_TESTS=ON ..
make -j2
CTEST_OUTPUT_ON_FAILURE=1 make test
