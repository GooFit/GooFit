set -evx
cd "${DEPS_DIR}"
ROOT_URL="https://root.cern.ch/download/root_v6.08.04.Linux-ubuntu14-x86_64-gcc4.8.tar.gz"
 
if [[ ! -f "${DEPS_DIR}/root/bin/root-config" ]] ; then
  echo "Downloading Root"
  mkdir -p root
  travis_retry wget --no-check-certificate --quiet -O - "${ROOT_URL}" | tar --strip-components=1 -xz -C root
fi

source "${DEPS_DIR}/root/bin/thisroot.sh"
set +evx
