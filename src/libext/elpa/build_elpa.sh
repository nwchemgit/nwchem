#!/usr/bin/env bash
#set -v
arch=`uname -m`
#SHORTVERSION=2020.11.001
#ERSION=new_release_2020.11.001
#SHORTVERSION=2021.05.002
#VERSION=new_release_2021_05_002
SHORTVERSION=2021.11.001
VERSION=new_release_2021.11.001
#https://gitlab.mpcdf.mpg.de/elpa/elpa/-/archive/new_release_2020.11.001/elpa-new_release_2020.11.001.tar.gz
export ARFLAGS=rU
if [ -f  elpa-${VERSION}.tar.gz ]; then
    echo "using existing"  elpa-${VERSION}.tar.gz
else
    rm -rf elpa*
echo    curl -L https://gitlab.mpcdf.mpg.de/elpa/elpa/-/archive/${VERSION}/elpa-${VERSION}.tar.gz -o elpa-${VERSION}.tar.gz
    curl -L https://gitlab.mpcdf.mpg.de/elpa/elpa/-/archive/${VERSION}/elpa-${VERSION}.tar.gz -o elpa-${VERSION}.tar.gz
fi
tar xzf elpa-${VERSION}.tar.gz
ln -sf elpa-${VERSION} elpa
cd elpa
UNAME_S=$(uname -s)
if [[ ${UNAME_S} == Darwin ]]; then
    export FORTRAN_CPP=$(find  /usr/local/Cellar/gcc/`brew list --versions gcc|cut -c 5-`/bin -name cpp*)
    if ! [ -x "$(command -v $FORTRAN_CPP)" ]; then
	echo
	echo cpp from gcc homebrew missing
	echo homebrew cpp is required for building Elpa
	echo
	exit 1
    else
	echo FORTRAN_CPP is $FORTRAN_CPP
    fi
fi
if [[ "$FC" = "ftn"  ]] ; then
    MPIF90="ftn"
    MPICC="cc"
else
    if ! [ -x "$(command -v mpif90)" ]; then
	echo
	echo mpif90 not installed
	echo mpif90 is required for building Elpa
	echo
	exit 1
    else
	MPIF90="mpif90"
        MPICC=mpicc
	#fix include path
	FCFLAGS+="-I`${NWCHEM_TOP}/src/tools/guess-mpidefs --mpi_include`"
	CFLAGS+="-I`${NWCHEM_TOP}/src/tools/guess-mpidefs --mpi_include`"
    fi
fi
if [[  -z "${FC}" ]]; then
    FC=$($MPIF90 -show|cut -d " " -f 1)
fi

if [[  -z "${BLAS_SIZE}" ]]; then
   BLAS_SIZE=8
fi
if [[ ${BLAS_SIZE} == 8 ]]; then
  sixty4_int+=" --enable-64bit-integer-math-support "

else
  sixty4_int+=" "
fi

if [[   -z "${CC}" ]]; then
    CC=cc
fi
if [[ ${FC} == flang ]] || [[ ${PE_ENV} == AOCC ]]; then
    GOTCLANG=1
else
    GOTCLANG=$( "$MPICC" -dM -E - </dev/null 2> /dev/null |grep __clang__|head -1|cut -c19)
fi
if [[ ${GOTCLANG} == "1" ]] ; then
#    if [[ ${UNAME_S} == Linux ]]; then
#	export FORTRAN_CPP=/usr/bin/cpp
#    fi
    CFLAGS+=" -Wno-error=implicit-function-declaration "
fi
# check gfortran version for arg check
GFORTRAN_EXTRA=$(echo $FC | cut -c 1-8)
if [[ ${GFORTRAN_EXTRA} == gfortran ]] || [[ ${PE_ENV} == GNU ]] || [[ ${FC} == flang ]] || [[ ${PE_ENV} == AOCC ]]; then
    let GFOVERSIONGT7=$(expr `${FC} -dumpversion | cut -f1 -d.` \> 7)
    if [[ ${GFOVERSIONGT7} == 1 ]]; then
	FCFLAGS+=' -std=legacy '
    fi
  sixty4_int+=" --disable-mpi-module "
fi
if [[ ${FC} == nvfortran ]]  || [[ ${PE_ENV} == NVIDIA ]] ; then
    sixty4_int+=" --disable-mpi-module "
    FCFLAGS+=" -fPIC"
    CFLAGS+=" -fPIC"
fi
if [[ ${FC} == ifort ]] || [[ ${FC} == ifx ]] || [[ ${PE_ENV} == INTEL ]] ; then
    FCFLAGS+=' -fpp'
    FCFLAGS+=" -fPIC"
    CFLAGS+=" -fPIC"
    #force CC=gcc
    export I_MPI_CC=gcc
    export I_MPI_FC=ifort
    export CC=gcc
    MYLINK+=" -fPIC"
fi

if [[  -z "${FORCETARGET}" ]]; then
FORCETARGET="-disable-sse -disable-sse-assembly --disable-avx --disable-avx2  --disable-avx512  "
fi #FORCETARGET
if [[ "${USE_HWOPT}" == "1" ]] && [[ "${USE_HWOPT}" == "y" ]] &&[[ "${USE_HWOPT}" != "Y" ]] && [[ ${UNAME_S} == Linux ]]; then
if [[ ${UNAME_S} == Linux ]]; then
    CPU_FLAGS=$(cat /proc/cpuinfo | grep flags |tail -n 1)
    CPU_FLAGS_2=$(cat /proc/cpuinfo | grep flags |tail -n 1)
elif [[ ${UNAME_S} == Darwin ]]; then
    CPU_FLAGS=$(/usr/sbin/sysctl -n machdep.cpu.features)
    CPU_FLAGS_2=$(/usr/sbin/sysctl -n machdep.cpu.leaf7_features)
fi
  GOTSSE2=$(echo ${CPU_FLAGS}   | tr  'A-Z' 'a-z'| awk ' /sse2/   {print "Y"}')
   GOTAVX=$(echo ${CPU_FLAGS}   | tr  'A-Z' 'a-z'| awk ' /avx/    {print "Y"}')
  GOTAVX2=$(echo ${CPU_FLAGS_2} | tr  'A-Z' 'a-z'| awk ' /avx2/   {print "Y"}')
GOTAVX512=$(echo ${CPU_FLAGS}   | tr  'A-Z' 'a-z'| awk ' /avx512f/{print "Y"}')
GOTCLZERO=$(echo ${CPU_FLAGS}   | tr  'A-Z' 'a-z'| awk ' /clzero/{print "Y"}')
if [[ ${CC} == icc ]] ; then
    CFLAGS+=" -xhost "
elif [[ ${CC} == nvc ]]  || [[ ${PE_ENV} == NVIDIA ]] ; then
    CFLAGS+=" -tp native"
elif [[ ${CC} == gcc ]] || [[ ${GOTCLANG} == "1" ]] || [[ ${CC} == cc ]]; then
    CFLAGS+=" -mtune=native -march=native "
fi    
    if [[ "${GOTAVX}" == "Y" ]]; then
	echo "using AVX instructions"
	FORCETARGET=" --disable-sse-assembly --enable-avx --disable-avx2  --disable-avx512  "
    fi
    if [[ "${GOTAVX2}" == "Y" ]]; then
	echo "using AVX2 instructions"
	FORCETARGET=" --enable-sse-assembly --enable-avx --enable-avx2  --disable-avx512  "
#	CFLAGS+=" -mmmx -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -maes -mavx -mfma -mavx2 "
    fi
    if [[ "${GOTAVX512}" == "Y" ]]; then
	echo "using AVX512 instructions"
	FORCETARGET=" --disable-sse-assembly --enable-avx --enable-avx2  --enable-avx512 "
    fi
#    if [[ ! -z ${USE_KNL} ]]; then
#	if [[ ${FC} == ifort ]] || [[ ${PE_ENV} == INTEL ]]; then
#	echo "using compiling for KNL "
#	sixty4_int+=" --enable-cross-compile --host=x86_64-unknown-linux-gnu "
#	CFLAGS+=" -xMIC-AVX512 "
#	FCFLAGS+=" -xMIC-AVX512 "
##	FORCETARGET=" --disable-sse-assembly --disable-avx --disable-avx2  --enable-avx512 "
#	FORCETARGET=" "
#	fi
#    fi
fi #USE_HWOPT
if [[ `${CC} -dM -E - < /dev/null 2> /dev/null | grep -c GNU` > 0 ]] ; then
    if [[ "$(expr `${CC} -dumpversion | cut -f1 -d.` \< 8)" == 1 ]]; then
	echo
	echo you have gcc version $(${CC} -dumpversion | cut -f1 -d.)
	echo gcc version 8 and later needed for elpa
	echo
	exit 1
    fi
fi
## check gcc version for skylake
#let GCCVERSIONGT5=$(expr `${CC} -dumpversion | cut -f1 -d.` \> 5)
#if [[ "$FORCETARGET" == *"SKYLAKEX"* ]]; then
#    if [[ ${GCCVERSIONGT5} != 1 ]]; then
#	echo
#	echo you have gcc version $(${CC} -dumpversion | cut -f1 -d.)
#	echo gcc version 6 and later needed for skylake
#	echo
#	exit 1
#    fi
#fi

if [ ! -f  configure ]; then
    sh ./autogen.sh
fi    
mkdir -p build
cd build
if [[ !  -z "${BUILD_SCALAPACK}" ]]; then
    MYLINK+=" -L${NWCHEM_TOP}/src/libext/lib -lnwc_scalapack"
fi
if [[ !  -z "${BUILD_OPENBLAS}" ]]; then
    MYLINK+=" -L${NWCHEM_TOP}/src/libext/lib -lnwc_openblas -lpthread"
fi
if [[ !  -z "${SCALAPACK_LIB}" ]]; then
    MYLINK+=" ${SCALAPACK_LIB} "
fi
if [[ !  -z "${LAPACK_LIB}" ]]; then
    MYLINK+=" ${LAPACK_LIB} "
fi
if [[ !  -z "${BLASOPT}" ]]; then
    MYLINK+=" ${BLASOPT} "
fi
export CFLAGS
export FCFLAGS
echo FCFLAGS is $FCFLAGS
echo CFLAGS is $CFLAGS
echo 64ints is $sixty4_int
echo MYLINK is "${MYLINK}"
export SCALAPACK_FCFLAGS="${MYLINK}"
export SCALAPACK_LDFLAGS="${MYLINK}"
export LIBS="${MYLINK}"
export    FC=$MPIF90
export CC=$MPICC
../configure \
    $sixty4_int \
  --disable-option-checking \
 --disable-dependency-tracking \
 --disable-shared --enable-static  \
 ${FORCETARGET} \
--prefix=${NWCHEM_TOP}/src/libext
unset LIBS
unset FCFLAGS
unset CFLAGS
#unset I_MPI_CC
#unset CC
#unset FC
#unset MPIF90
#unset MPICC
unset SCALAPACK_FCFLAGS
unset SCALAPACK_LDFLAGS
if [[  -z "$MPICH_FC"   ]] ; then
    MPICH_FC="$FC"
    echo MPICH_FC is nvfortran
fi    
if [[ "${FORTRAN_CPP}" != "" ]] ; then
    make V=0 -j4
    make V=0 -j4 install
#    make V=0 -j4 install-libLTLIBRARIES install-data
else    
    top_srcdir=`pwd`/..
    make V=0 FC="${top_srcdir}/remove_xcompiler ${top_srcdir}/manual_cpp $MPIF90" -j4
    make V=0 FC="${top_srcdir}/remove_xcompiler ${top_srcdir}/manual_cpp $MPIF90" -j4 install
fi

if [[ "$?" != "0" ]]; then
    echo " "
    echo "Elpa compilation failed"
    echo " "
    exit 1
fi
ln -sf ${NWCHEM_TOP}/src/libext/lib/libelpa.a  ${NWCHEM_TOP}/src/libext/lib/libnwc_elpa.a
ln -sf ${NWCHEM_TOP}/src/libext/include/elpa-${SHORTVERSION}  ${NWCHEM_TOP}/src/libext/include/elpa
