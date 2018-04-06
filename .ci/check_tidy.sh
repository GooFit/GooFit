#!/usr/bin/env bash

echo -en 'travis_fold:start:script.build\\r'
echo "Building..."
set -evx

mkdir -p build
cd build
cmake .. -DGOOFIT_TIDY_FIX=ON
cmake --build .

set +evx
echo -e 'travis_fold:end:script.build\\r'
echo -en 'travis_fold:start:script.diff\\r'
echo "Diff..."
set -evx

git diff --exit-code --color

set +evx
echo -e 'travis_fold:end:script.diff\\r'
