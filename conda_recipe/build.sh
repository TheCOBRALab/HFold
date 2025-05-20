#!/usr/bin/env bash
set -euo pipefail

# ──────────────────────────────────────────────────────────────
# 0 · Global compiler / pre-processor flags
# ──────────────────────────────────────────────────────────────

export CXXFLAGS="${CXXFLAGS:-} -std=c++17 -include cstdint"
if [[ "$OSTYPE" == "darwin"* ]]; then
  export CXXFLAGS="${CXXFLAGS} -DNI_MAXHOST=1025"
fi
export CFLAGS="${CFLAGS:-}  -std=c11"
export CPPFLAGS="${CPPFLAGS:-} -D_POSIX_C_SOURCE=200809L -D_GNU_SOURCE"

# where (if any) little helper *.py files should end up
: "${PY_VER:=}"                       # silence "unbound variable" if undefined
PY_SITE_DIR="${PREFIX}/lib/python${PY_VER:-no}/site-packages/ViennaRNA"

# ──────────────────────────────────────────────────────────────
# 1 · Fresh build directory
# ──────────────────────────────────────────────────────────────
rm -rf build
mkdir  build
cd     build

# ──────────────────────────────────────────────────────────────
# 2 · Unpack ViennaRNA and fix two outdated prototypes
# ──────────────────────────────────────────────────────────────
curl -L -O https://github.com/ViennaRNA/ViennaRNA/releases/download/v2.7.0/ViennaRNA-2.7.0.tar.gz
tar -xf  ./ViennaRNA-2.7.0.tar.gz
pushd   ViennaRNA-2.7.0 >/dev/null

# Detect OS and set sed in-place flag
if [[ "$OSTYPE" == "darwin"* ]]; then
  # macOS/BSD sed
  SED_INPLACE=(-i '')
else
  # Linux/GNU sed
  SED_INPLACE=(-i)
fi

sed -E "${SED_INPLACE[@]}" 's@void[[:space:]]+lst_kill[^(]*\(LIST[[:space:]]*\*[[:space:]]*l[[:space:]]*,[[:space:]]*void[^(]*\([^)]*\)\)@void lst_kill (LIST *l, void (*freeNode)(void *))@' \
    src/ViennaRNA/datastructures/lists.h

sed -E "${SED_INPLACE[@]}" 's@void[[:space:]]+lst_mergesort[^(]*\(LIST[[:space:]]*\*[[:space:]]*l[[:space:]]*,[[:space:]]*int[^(]*\([^)]*\)\)@void lst_mergesort (LIST *l, int (*cmp_func)(void *, void *))@' \
    src/ViennaRNA/datastructures/lists.h

# ─  Make dlib’s Unicode code libc++-friendly  ─
sed -E "${SED_INPLACE[@]}" 's@typedef[[:space:]]+uint32_t[[:space:]]+unichar;@using unichar = char32_t;@' \
        src/dlib-19.24/dlib/unicode/unicode.h

# ──────────────────────────────────────────────────────────────
# 3 · Configure, build, *install*  (Python bindings OFF)
#    We still let the helper scripts be installed, just not into /RNA
# ──────────────────────────────────────────────────────────────
./configure --without-perl --without-python --prefix="${PREFIX}"

make  -j"${CPU_COUNT}"
# override pkgpyexecdir so "make install" never touches /RNA
make  -j"${CPU_COUNT}" pkgpyexecdir="${PY_SITE_DIR}" install
popd >/dev/null   # ← back to build/

# ──────────────────────────────────────────────────────────────
# 4 · CMake phase for your own project that links ViennaRNA
# ──────────────────────────────────────────────────────────────
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="${CC}" \
  -DCMAKE_CXX_COMPILER="${CXX}" \
  -DCMAKE_C_FLAGS="-DHAVE_STRDUP=1" \
  -DVIENNARNA_LIB="${PREFIX}" \
  -DVIENNARNA_INCLUDE_DIR="${PREFIX}/include" \
  -DENABLE_TESTS=OFF

cmake --build . --parallel -- VERBOSE=1
cmake --install . --prefix "${PREFIX}"
