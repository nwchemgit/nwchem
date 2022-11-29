#!/usr/bin/env bash
source ../libext_utils/cmake.sh

check_tgz() {
    myexit=0
    [ -f $1 ] && gunzip -t $1 > /dev/null && myexit=1
    echo $myexit
}

VERSION=ilp64
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
    curl  -sS -L https://github.com/dmejiar/tblite/tarball/${VERSION} -o $TGZ
fi

tar -xzf tblite-${VERSION}.tar.gz
ln -sf *tblite-??????? tblite


if [[ -z "${MPIF90}" ]]; then
  if [[ "$FC" = "ftn"  ]] ; then
    MPIF90="ftn"
    MPICC="cc"
  else
    if ! [ -x "$(command -v mpif90)" ]; then
	echo
	echo mpif90 not installed
	echo mpif90 is required for building Scalapack
	echo
	exit 1
    else
	MPIF90="mpif90"
        MPICC=mpicc
    fi
  fi
fi

if [[  -z "${FC}" ]]; then
    FC=$($MPIF90 -show|cut -d " " -f 1)
fi

if [[  -z "${NWCHEM_TOP}" ]]; then
    dir3=$(dirname `pwd`)
    dir2=$(dirname "$dir3")
    NWCHEM_TOP=$(dirname "$dir2")
fi


if [[ -n ${FC} ]] &&   [[ ${FC} == ftn ]]; then
    if [[ ${PE_ENV} == PGI ]]; then
          FC=pgf90
    fi
    if [[ ${PE_ENV} == INTEL ]]; then
	FC=ifort
    fi
    if [[ ${PE_ENV} == GNU ]]; then
	FC=gfortran
    fi
    if [[ ${PE_ENV} == AOCC ]]; then
	FC=flang
    fi
    if [[ ${PE_ENV} == NVIDIA ]]; then
	FC=nvfortran
    fi
    if [[ ${PE_ENV} == CRAY ]]; then
	FC=crayftn
	CC=clang
	#fix for libunwind.so link problem
        export LD_LIBRARY_PATH=/opt/cray/pe/cce/$CRAY_FTN_VERSION/cce-clang/x86_64/lib:/opt/cray/pe/lib64/cce/:$LD_LIBRARY_PATH
    fi
fi

FC_EXTRA=$(${NWCHEM_TOP}/src/config/strip_compiler.sh ${FC})

#Intel MPI
if [[  -z "$I_MPI_F90"   ]] ; then
    export I_MPI_F90="$FC"
    echo I_MPI_F90 is "$I_MPI_F90"
fi
if [[  -z "$PE_ENV"   ]] ; then
    #check if mpif90 and FC are consistent
    MPIF90_EXTRA=$(${NWCHEM_TOP}/src/config/strip_compiler.sh `${MPIF90} -show`)
    if [[ $MPIF90_EXTRA != $FC_EXTRA ]]; then
        echo which mpif90 is `which mpif90`
        echo mpif90show `${MPIF90} -show`
	echo FC and MPIF90 are not consistent
	echo FC is $FC_EXTRA
	echo MPIF90 is $MPIF90_EXTRA
	exit 1
    fi
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
  if  [[ ${FC} == f95 ]] || [[ ${FC_EXTRA} == gfortran ]] ; then
      Fortran_FLAGS=" -fdefault-integer-8 -w "
    elif  [[ ${FC} == xlf ]] || [[ ${FC} == xlf_r ]] || [[ ${FC} == xlf90 ]]|| [[ ${FC} == xlf90_r ]]; then
      Fortran_FLAGS=" -qintsize=8 -qextname "
    elif  [[ ${FC} == crayftn ]]; then
      Fortran_FLAGS=" -s integer64 -h nopattern"
    elif  [[ ${FC} == frtpx ]] || [[ ${FC} == frt ]]; then
      Fortran_FLAGS=" -fs -CcdLL8 -CcdII8 "
    else
      Fortran_FLAGS=" -i8 "
    fi
else
  ilp64=OFF
  Fortran_FLAGS=""
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
if [[ ${FC} == nvfortran ]] || [[ ${FC} == pgf90 ]]; then
  Fortran_FLAGS+="-Mbackslash -fast -tp host"
fi

if [[ -z "$USE_OPENMP" ]]; then
  DOOPENMP=OFF
else
  DOOPENMP=ON
fi
# 2022 Intel compilers generate buggy code when USE_OPENMP=1
if [[ ${FC} == ifort ]] || [[ ${FC} == ifx ]]; then
    IFORTVER=$(ifort -v 2>&1|cut -d " " -f 3)
    IFORTVER_YEAR=$(echo $IFORTVER | cut -d . -f 1)
    if [[ "$IFORTVER_YEAR" -gt "2021" ]]; then
	DOOPENMP=OFF
    fi
    if [[ ${FC} -eq "ifx" ]]; then
	DOOPENMP=OFF
    fi
fi



cd tblite
rm -rf _build

echo compiling TBlite stack with FC=$FC CC=$CC $CMAKE -B _build -DLAPACK_LIBRARIES="$BLASOPT" -DWITH_ILP64=$ilp64 -DWITH_OpenMP=$DOOPENMP -DCMAKE_INSTALL_PREFIX="../.." -DWITH_TESTS=OFF -DWITH_API=OFF -DWITH_APP=OFF -DCMAKE_INSTALL_LIBDIR="lib" -DCMAKE_IGNORE_PATH="/usr/local" -DCMAKE_Fortran_FLAGS="$Fortran_FLAGS"

FC=$FC CC=$CC $CMAKE -B _build -DLAPACK_LIBRARIES="$BLASOPT" -DWITH_ILP64=$ilp64 -DWITH_OpenMP=$DOOPENMP -DCMAKE_INSTALL_PREFIX="../.." -DWITH_TESTS=OFF -DWITH_API=OFF -DWITH_APP=OFF -DCMAKE_INSTALL_LIBDIR="lib" -DCMAKE_IGNORE_PATH="/usr/local" -DCMAKE_Fortran_FLAGS="$Fortran_FLAGS"
$CMAKE --build _build --parallel 4
status=$?
if [ $status -ne 0 ]; then
    echo tblite compilation failed
    exit 1
fi

$CMAKE --install _build

cd ..

ln -sf  ../lib/libtblite.a  ../lib/libnwc_tblite.a
