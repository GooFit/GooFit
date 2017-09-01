#!/usr/bin/env sh
set -evx

clang-format --version

git ls-files -- '*.cpp' '*.hpp' '*.cu' '*.h' '*.cc' | xargs clang-format -i -style=file

git diff --exit-code --color

mkdir build || true
cd build
cmake .. -DGOOFIT_TIDY_FIX=ON
cmake --build .

git diff --exit-code --color

set +evx
