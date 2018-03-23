#!/usr/bin/env bash

echo -en 'travis_fold:start:script.build\\r'
echo "Building..."
set -evx

mkdir -p build
cd build
cmake -DGOOFIT_DEVICE=OMP -DGOOFIT_TESTS=ON -DCMAKE_BUILD_TYPE=Coverage -DGOOFIT_PYTHON=OFF ..
cmake --build . -- -j2
cmake --build . --target GooFit_coverage

set +evx
echo -en 'travis_fold:end:script.build\\r'
echo -en 'travis_fold:start:script.lcov\\r'
echo "Running lcov and upload..."
set -evx

lcov --directory . --capture --output-file coverage.info # capture coverage info
lcov --remove coverage.info '*/tests/*' '*gtest*' '*gmock*' '/usr/*' --output-file coverage.info # filter out system
lcov --list coverage.info #debug info
# Uploading report to CodeCov
bash <(curl -s https://codecov.io/bash) || echo "Codecov did not collect coverage reports"

set +evx
echo -en 'travis_fold:end:script.lcov\\r'
