#!/usr/bin/env bash

set -evx

pushd "${DEPS_DIR}"
ROOT_URL="https://root.cern/download/root_v6.18.04.Linux-ubuntu18-x86_64-gcc7.4.tar.gz"

if [[ ! -f "${DEPS_DIR}/root/bin/root-config" ]] ; then
  echo "Downloading Root"
  mkdir -p root
  travis_retry wget --no-check-certificate --quiet -O - "${ROOT_URL}" | tar --strip-components=1 -xz -C root
fi

source "${DEPS_DIR}/root/bin/thisroot.sh"
popd

set +evx
