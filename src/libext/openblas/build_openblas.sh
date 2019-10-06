#!/bin/bash 
VERSION=0.3.7
rm -rf OpenBLAS*
curl -LJO https://github.com/xianyi/OpenBLAS/archive/v${VERSION}.tar.gz
tar xzf OpenBLAS-${VERSION}.tar.gz
ln -sf OpenBLAS-${VERSION} OpenBLAS
cd OpenBLAS-${VERSION}
if [[ ! -z "${USE_64TO32}" ]]; then
  sixty4_int=0
else
  sixty4_int=1
fi
if [[ "${NWCHEM_TARGET}" == "LINUX" ]]; then
  binary=32
  sixty4_int=0
else
  binary=64
fi
make  INTERFACE64="$sixty4_int" BINARY="$binary" USE_THREAD=0 NO_CBLAS=1 NO_LAPACKE=1 DEBUG=0 NUM_THREADS=1 libs netlib
mkdir -p ../../lib
cp libopenbla*.* ../../lib
#make PREFIX=. install
#tail -4 openblas.log
