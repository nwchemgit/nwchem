#!/usr/bin/env bash
set -v
arch=`uname -m`
VERSION=0.2.0-ilp64-alpha
if [ -f  tblite-${VERSION}.tar.gz ]; then
    echo "using existing"  tblite-${VERSION}.tar.gz
else
    rm -rf tblite*
    curl -L https://github.com/dmejiar/tblite/archive/v${VERSION}.tar.gz -o tblite-${VERSION}.tar.gz
fi

tar xzf tblite-${VERSION}.tar.gz
ln -sf tblite-${VERSION} tblite
cd tblite
rm -rf _build

if [[  -z "${BLAS_SIZE}" ]]; then
   BLAS_SIZE=8
fi
if [[ ${BLAS_SIZE} == 8 ]]; then
  ilp64=true
else
  ilp64=false
fi

TEST=$( echo $BLASOPT | cut -f1 -d" " | sed 's/-L//' )

SLP=$LIBRARY_PATH
LIBRARY_PATH=$TEST

meson setup _build -Dlapack=mkl -Dilp64=$ilp64 -Ddftd4:ilp64=$ilp64 -Dmulticharge:ilp64=$ilp64 -Ds-dftd3:ilp64=$ilp64
meson compile -C _build

LIBRARY_PATH=$SLP

mkdir -p ../../lib
cp _build/libtblite.a ../../lib/libnwc_tblite.a
