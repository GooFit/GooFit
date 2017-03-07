#!/usr/bin/env sh
set -evx

cd ${TRAVIS_BUILD_DIR}

mkdir build || true
cd build
cmake -DGOOFIT_DEVICE=OMP -DGOOFIT_SEPARATE_COMP=ON -DGOOFIT_TESTS=ON -DCMAKE_BUILD_TYPE=Coverage ..
cmake --build .
cmake --build . --target CLI_coverage

lcov --directory . --capture --output-file coverage.info # capture coverage info
lcov --remove coverage.info '*/tests/*' '*gtest*' '*gmock*' '/usr/*' --output-file coverage.info # filter out system
lcov --list coverage.info #debug info
# Uploading report to CodeCov
bash <(curl -s https://codecov.io/bash) || echo "Codecov did not collect coverage reports"
set +evx
