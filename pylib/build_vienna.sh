#!/bin/bash
set -e

if [ -f /tmp/vienna-cache/lib/libRNA.a ]; then
  echo "✅ Using cached ViennaRNA from /tmp/vienna-cache"
else
  echo "🔧 Building ViennaRNA from source..."
  curl -L -O https://github.com/ViennaRNA/ViennaRNA/releases/download/v2.7.0/ViennaRNA-2.7.0.tar.gz
  tar -xzf ViennaRNA-2.7.0.tar.gz
  cd ViennaRNA-2.7.0

  # Bug fixes in dlib
  if [[ "$OSTYPE" == "darwin"* ]]; then
      # typedef with the standard UTF-32 code unit.
      sed -i.bak 's/typedef uint32 unichar;/typedef char32_t unichar;/' \
          src/dlib-19.24/dlib/unicode/unicode.h

      # make the new type serialisable
      sed -i.bak '/USE_DEFAULT_INT_SERIALIZATION_FOR(uint64)/a\
      USE_DEFAULT_INT_SERIALIZATION_FOR(char32_t)
      ' src/dlib-19.24/dlib/serialize.h

      # “template” keyword.  Add the empty <> that upstream already has.
      sed -i.bak 's/::template go(/::template go<>(/' \
          src/dlib-19.24/dlib/global_optimization/find_max_global.h
  fi

  ./configure --without-perl --without-python --disable-shared --enable-static --prefix=/tmp/vienna-cache
  make -j$(nproc)

  # override pkgpyexecdir so "make install" never touches /RNA (/RNA needs sudo access)
  make  -j"${CPU_COUNT}" pkgpyexecdir="${PY_SITE_DIR}" install

  cd ..
  rm -rf ViennaRNA-2.7.0
  rm -f ViennaRNA-2.7.0.tar.gz
fi

# Set paths for build
export VIENNARNA_LIB=/tmp/vienna-cache/lib
export VIENNARNA_INCLUDE_DIR=/tmp/vienna-cache/include
