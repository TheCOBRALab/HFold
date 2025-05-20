#!/bin/bash
rm -rf build
mkdir -p build
cd build

cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER="${CXX}" \
  -DCMAKE_C_FLAGS="-DHAVE_STRDUP=1" \
  -DVIENNARNA_LIB="${PREFIX}" \
  -DVIENNARNA_INCLUDE_DIR="${PREFIX}" \
  -DENABLE_TESTS=OFF \

cmake --build . --parallel -- VERBOSE=1
cmake --install . --prefix="${PREFIX}"
