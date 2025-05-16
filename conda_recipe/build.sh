#!/bin/bash
rm -rf build
mkdir -p build
cd build

find "${VIENNARNA_LIB}" -name libRNA.a
find "${VIENNARNA_INCLUDE_DIR}"-name pair_mat.h

cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER="${CXX}" \
  -DCMAKE_C_FLAGS="-DHAVE_STRDUP=1" \
  -DVIENNARNA_LIB="${VIENNARNA_LIB}" \
  -DVIENNARNA_INCLUDE_DIR="${VIENNARNA_INCLUDE_DIR}"

cmake --build . --parallel -- VERBOSE=1
cmake --install . --prefix="${PREFIX}"
