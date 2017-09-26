#!/usr/bin/env sh
set -evx

clang-format-3.9 --version

git ls-files -- '*.cpp' '*.hpp' '*.cu' '*.cc' | xargs clang-format-3.9 -style=file -i

git diff --exit-code --color

mkdir build || true
cd build
cmake .. -DGOOFIT_TIDY_FIX=ON
cmake --build .

git diff --exit-code --color

set +evx
