#!/usr/bin/env sh

set -e

if [ -x "$(command -v clang-format)" ] ; then

    clang-format --version
    git ls-files -- '*.cpp' '*.hpp' '*.cu' '*.h' | xargs clang-format -style=file -sort-includes -i

elif [ -x "$(command -v docker)" ] ; then

    docker run -it unibeautify/clang-format --version
    docker run -it -v "$(pwd)":/workdir -w /workdir unibeautify/clang-format -style=file -sort-includes -i $(git ls-files -- '*.cpp' '*.hpp' '*.cu' '*.h')

fi

git diff --exit-code --color

set +e

