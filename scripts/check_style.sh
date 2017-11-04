#!/usr/bin/env sh
set -evx

clang-format --version

git ls-files -- '*.cpp' '*.hpp' '*.cu' '*.cc' '*.h' | xargs clang-format -style=file -i

git diff --exit-code --color

set +evx
