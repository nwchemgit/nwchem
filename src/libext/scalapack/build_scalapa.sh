#!/bin/bash
if [[ "$FC" = "ftn"  ]] ; then
    MPIF90="ftn"
else
    if ! [ -x "$(command -v mpif90)" ]; then
	echo
	echo mpif90 not installed
	echo mpif90 is required for building Scalapack
	echo
	exit 1
    else
	MPIF90="mpif90"
    fi
fi
if [[  -z "${NWCHEM_TOP}" ]]; then
    dir3=$(dirname `pwd`)
    dir2=$(dirname "$dir3")
    NWCHEM_TOP=$(dirname "$dir2")
fi
if [[ "$SCALAPACK_SIZE" != "4"  ]] ; then
    echo SCALAPACK_SIZE must be equal to 4
    exit 1
fi
if [[ "$BLAS_SIZE" != "4"  ]] ; then
    echo BLAS_SIZE must be equal to 4 for SCALAPACK
    exit 1
fi
if [[ -z "$USE_64TO32"   ]] ; then
    echo USE_64TO32 must be set
    exit 1
fi
if [[ ! -z "$BUILD_OPENBLAS"   ]] ; then
    BLASOPT="-L`pwd`/../lib -lopenblas"
fi
#git clone https://github.com/scibuilder/scalapack.git
#svn co --non-interactive --trust-server-cert https://icl.utk.edu/svn/scalapack-dev/scalapack/trunk/ scalapack
rm -rf scalapack*
VERSION=2.1.0
curl -L https://github.com/Reference-ScaLAPACK/scalapack/archive/v${VERSION}.tar.gz -o scalapack.tgz
#curl -L http://www.netlib.org/scalapack/scalapack-${VERSION}.tgz -o scalapack.tgz
tar xzf scalapack.tgz
if [[ ${VERSION} == 2.0.2 ]] ; then
    patch -p0 < mpistruct.patch
fi
ln -sf scalapack-${VERSION} scalapack
mkdir -p scalapack/build
cd scalapack/build
if  [[ -n ${FC} ]] &&   [[ ${FC} == xlf ]] || [[ ${FC} == xlf_r ]] || [[ ${FC} == xlf90 ]]|| [[ ${FC} == xlf90_r ]]; then
    Fortran_FLAGS=" -qintsize=4 -qextname "
elif [[ -n ${FC} ]] &&   [[ ${FC} == ftn ]]; then
    Fortran_FLAGS="-O2 -g -axCORE-AVX2"
#elif [[ -n ${FC} ]] &&   [[ ${FC} == flang ]]; then
# unset FC=flang since cmake gets lost
#       unset FC
fi
if [[ ! -z "$BUILD_SCALAPACK"   ]] ; then
    Fortran_FLAGS+=-I"$NWCHEM_TOP"/src/libext/include
fi
FC=$MPIF90 FFLAGS="$Fortran_FLAGS" cmake -Wno-dev ../ -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_Fortran_FLAGS="$Fortran_FLAGS" -DTEST_SCALAPACK=OFF  -DBUILD_TESTING=OFF -DBUILD_SHARED_LIBS=OFF  -DBLAS_openblas_LIBRARY="$BLASOPT"  -DBLAS_LIBRARIES="$BLASOPT"  -DLAPACK_openblas_LIBRARY="$BLASOPT"  -DLAPACK_LIBRARIES="$BLASOPT" 
make V=0 -j3 scalapack/fast
mkdir -p ../../../lib
cp lib/libscalapack.a ../../../lib
