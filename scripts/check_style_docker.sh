#!/usr/bin/env sh

# Also good but untagged: CLANG_FORMAT=unibeautify/clang-format
CLANG_FORMAT=saschpe/clang-format:5.0.1

set -evx

docker run -it --rm ${CLANG_FORMAT} --version
docker run -it --rm -v "$(pwd)":/workdir -w /workdir ${CLANG_FORMAT} -style=file -sort-includes -i $(git ls-files -- '*.cpp' '*.hpp' '*.cu' '*.h')

git diff --exit-code --color

set +evx
