#!/usr/bin/env bash
source ../libext_utils/cmake.sh

check_tgz() {
    myexit=0
    [ -f $1 ] && gunzip -t $1 > /dev/null && myexit=1
    echo $myexit
}

VERSION=0.2.6-ilp64-alpha
TGZ=tblite-${VERSION}.tar.gz
if [ ! -z "${USE_INTERNALBLAS}" ]; then
    echo USE_TBLITE not compatible with USE_INTERNALBLAS
    echo Please set BUILD_OPENBLAS or
    echo BLASOPT/LAPACK_LIB
    exit 1
fi
if [ `check_tgz $TGZ` == 1 ]; then
    echo "using existing $TGZ"
else
    rm -rf tblite*
    curl -L https://github.com/dmejiar/tblite/archive/v${VERSION}.tar.gz -o $TGZ
fi

tar -xzf tblite-${VERSION}.tar.gz
ln -sf tblite-${VERSION} tblite


if [[  -z "${CC}" ]]; then
    CC=cc
fi
if [[  -z "${FC}" ]]; then
#FC not defined. Look for gfortran
    if [[ ! -x "$(command -v gfortran)" ]]; then
	echo ' '
	echo 'please define FC to compile tblite'
	echo ' '
	exit 1
    else
	echo 'FC not defined, defaulting FC=gfortran'
	FC=gfortran
    fi
fi

if [[ -z "${CMAKE}" ]]; then
    #look for cmake
    if [[ -z "$(command -v cmake)" ]]; then
	cmake_instdir=../libext_utils
	get_cmake_release $cmake_instdir
	status=$?
	if [ $status -ne 0 ]; then
	    echo cmake required to build tblite
	    echo Please install cmake
	    echo define the CMAKE env. variable
	    exit 1
	fi
    else
	CMAKE=cmake
    fi
fi
CMAKE_VER_MAJ=$(${CMAKE} --version|cut -d " " -f 3|head -1|cut -d. -f1)
CMAKE_VER_MIN=$(${CMAKE} --version|cut -d " " -f 3|head -1|cut -d. -f2)
echo CMAKE_VER is ${CMAKE_VER_MAJ} ${CMAKE_VER_MIN}
if ((CMAKE_VER_MAJ < 3)) || (((CMAKE_VER_MAJ > 2) && (CMAKE_VER_MIN < 11))); then
    get_cmake_release  $cmake_instdir
    status=$?
    if [ $status -ne 0 ]; then
	echo cmake required to build scalapack
	echo Please install cmake
	echo define the CMAKE env. variable
	exit 1
    fi
fi



if [[  -z "${BLAS_SIZE}" ]]; then
   BLAS_SIZE=8
fi
if [[ ${BLAS_SIZE} == 8 ]]; then
  ilp64=ON
else
  ilp64=OFF
fi

if [[ ! -z "$BUILD_OPENBLAS"   ]] ; then
    BLASOPT="-L`pwd`/../lib -lnwc_openblas -lpthread"
fi
# check gfortran version
FFLAGS_IN=" "
if [[ `${FC} -dM -E - < /dev/null 2> /dev/null | grep -c GNU` > 0 ]] ; then
    let GFORTRANVERSIONGT8=$(expr `${FC} -dumpversion | cut -f1 -d.` \> 8)
    if [[ ${GFORTRANVERSIONGT8} != 1 ]]; then
	echo
	echo you have gfortran version $(${FC} -dumpversion | cut -f1 -d.)
	echo gcc version 9 and later needed for tblite
	echo
	exit 1
    fi
fi
# stop if FC is flang from AOCC
if [[ ${FC} == flang ]]; then
    if [[ `${FC} -dM -E - < /dev/null 2> /dev/null | grep -c AOCC` > 0 ]] ; then
	echo
	echo flang from AOCC not compatible with tblite
	echo https://tblite.readthedocs.io/en/latest/installation.html#supported-compilers
	echo
	exit 1
    fi
fi
#nvfortran
if [[ ${FC} == nvfortran ]]; then
	echo
	echo nvfortran not compatible with tblite
	echo https://tblite.readthedocs.io/en/latest/installation.html#supported-compilers
	echo
	exit 1
fi
if [[ -z "$USE_OPENMP" ]]; then
  DOOPENMP=OFF
else
  DOOPENMP=ON
fi

cd tblite
rm -rf _build

FC=$FC CC=$CC $CMAKE -B _build -DLAPACK_LIBRARIES="$BLASOPT" -DWITH_ILP64=$ilp64 -DWITH_OpenMP=$DOOPENMP -DCMAKE_INSTALL_PREFIX="../.." -DWITH_TESTS=OFF -DWITH_API=OFF -DCMAKE_INSTALL_LIBDIR="lib" -DCMAKE_IGNORE_PATH=/usr/local
$CMAKE --build _build --parallel 4
status=$?
if [ $status -ne 0 ]; then
    echo tblite compilation failed
    exit 1
fi

$CMAKE --install _build

cd ..

ln -sf  ../lib/libtblite.a  ../lib/libnwc_tblite.a
