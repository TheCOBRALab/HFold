#!/bin/bash
set -e

if [ -f /tmp/vienna-cache/lib/libRNA.a ]; then
  echo "✅ Using cached ViennaRNA from /tmp/vienna-cache"
else
  echo "🔧 Building ViennaRNA from source..."
  curl -L -O https://github.com/ViennaRNA/ViennaRNA/releases/download/v2.7.0/ViennaRNA-2.7.0.tar.gz
  tar -xzf ViennaRNA-2.7.0.tar.gz
  cd ViennaRNA-2.7.0
  ./configure --without-perl --without-python --disable-shared --enable-static --prefix=/tmp/vienna-cache
  make -j$(nproc)
  make install
fi

# Set paths for build
export VIENNARNA_LIB=/tmp/vienna-cache/lib
export VIENNARNA_INCLUDE_DIR=/tmp/vienna-cache/include
