#!/bin/bash
set -e

if [ -f /tmp/vienna-cache/lib/libRNA.a ]; then
  echo "✅ Using cached ViennaRNA from /tmp/vienna-cache"
else
  echo "🔧 Building ViennaRNA from source..."
  curl -L -O https://github.com/ViennaRNA/ViennaRNA/releases/download/v2.7.1/ViennaRNA-2.7.1.tar.gz
  tar -xzf ViennaRNA-2.7.1.tar.gz
  cd ViennaRNA-2.7.1

  ./configure --without-perl --without-python --disable-shared --enable-static --prefix=/tmp/vienna-cache
  make -j$(nproc)

  # override pkgpyexecdir so "make install" never touches /RNA (/RNA needs sudo access)
  make  -j"${CPU_COUNT}" pkgpyexecdir="${PY_SITE_DIR}" install

  cd ..
  rm -rf ViennaRNA-2.7.1
  rm -f ViennaRNA-2.7.1.tar.gz
fi

# Set paths for build
export VIENNARNA_LIB=/tmp/vienna-cache/lib
export VIENNARNA_INCLUDE_DIR=/tmp/vienna-cache/include
