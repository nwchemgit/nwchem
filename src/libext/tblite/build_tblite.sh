#!/usr/bin/env bash
set -v
arch=`uname -m`
VERSION=0.2.4-ilp64-alpha
if [ ! -z "${USE_INTERNALBLAS}" ]; then
    echo USE_TBLITE not compatible with USE_INTERNALBLAS
    echo Please set BUILD_OPENBLAS or
    echo BLASOPT/LAPACK_LIB
    exit 1
fi
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

LAPACK_TYPE=mkl
if [[ ! -z "$BUILD_OPENBLAS"   ]] ; then
    BLASOPT="-L`pwd`/../lib -lnwc_openblas"
    LAPACK_TYPE=openblas
fi
TEST=$( echo $BLASOPT | cut -f1 -d" " | sed 's/-L//' )
SLP=$LIBRARY_PATH
LIBRARY_PATH=$TEST

meson setup _build -Dlapack=$LAPACK_TYPE  -Dilp64=$ilp64 -Ddftd4:ilp64=$ilp64 -Dmulticharge:ilp64=$ilp64 -Ds-dftd3:ilp64=$ilp64 --prefix=$NWCHEM_TOP/src/libext --libdir $NWCHEM_TOP/src/libext/lib --includedir $NWCHEM_TOP/src/libext/include
meson compile -C _build
meson install -C _build

LIBRARY_PATH=$SLP

cd ../lib
ln -sf  libtblite.a  libnwc_tblite.a
