#!/usr/bin/env sh
set -evx

mkdir build || true
cd build
cmake .. -DGOOFIT_TIDY_FIX=ON
cmake --build .

git diff --exit-code --color

set +evx
