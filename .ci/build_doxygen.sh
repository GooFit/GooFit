set -evx
DOXYGEN_URL="http://doxygen.nl/files/doxygen-1.8.17.linux.bin.tar.gz"
cd "${DEPS_DIR}"

if [[ ! -f "${DEPS_DIR}/doxygen/bin/doxygen" ]] ; then
  echo "Downloading Doxygen"
  mkdir -p doxygen
  travis_retry wget --no-check-certificate --quiet -O - "${DOXYGEN_URL}" | tar --strip-components=1 -xz -C doxygen
fi

export PATH="${DEPS_DIR}/doxygen/bin:${PATH}"

cd "${TRAVIS_BUILD_DIR}"
set +evx
