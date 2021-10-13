#!/usr/bin/env bash
#set -v
arch=`uname -m`
SHORTVERSION=2020.11.001
VERSION=new_release_2020.11.001
#SHORTVERSION=2021.05.002
#VERSION=new_release_2021_05_002
#https://gitlab.mpcdf.mpg.de/elpa/elpa/-/archive/new_release_2020.11.001/elpa-new_release_2020.11.001.tar.gz
if [ -f  elpa-${VERSION}.tar.gz ]; then
    echo "using existing"  elpa-${VERSION}.tar.gz
else
    rm -rf elpa*
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
    fi
fi
if [[  -z "${FC}" ]]; then
    FC=$($MPIF90 -show|cut -d " " -f 1)
fi

if [[  -z "${FORCETARGET}" ]]; then
FORCETARGET=" -disable-sse-assembly --disable-avx --disable-avx2  --disable-avx512  "
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
if [[ ${UNAME_S} == Linux ]]; then
    if [[ "${GOTAVX}" == "Y" ]]; then
	echo "using AVX instructions"
	FORCETARGET=" --enable-sse-assembly --enable-avx --disable-avx2  --disable-avx512  "
    fi
    if [[ "${GOTAVX2}" == "Y" ]]; then
	echo "using AVX2 instructions"
	FORCETARGET=" --enable-sse-assembly --enable-avx --enable-avx2  --disable-avx512  "
    fi
#    if [[ "${GOTAVX512}" == "Y" ]]; then
#	echo "using AVX512 instructions"
#	FORCETARGET=" --enable-sse-assembly --enable-avx --enable-avx2  --enable-avx512  "
#    fi
fi #Linux
fi #FORCETARGET
if [[  -z "${BLAS_SIZE}" ]]; then
   BLAS_SIZE=8
fi
if [[ ${BLAS_SIZE} == 8 ]]; then
  sixty4_int=" --enable-64bit-integer-math-support "

else
  sixty4_int=" "
fi

if [[   -z "${CC}" ]]; then
    CC=cc
fi
GOTCLANG=$( "$MPICC" -dM -E - </dev/null 2> /dev/null |grep __clang__|head -1|cut -c19)
if [[ ${GOTCLANG} == "1" ]] ; then
    sixty4_int+='CFLAGS=-Wno-error=implicit-function-declaration '
#    C_FLAGS=" -Wno-error=implicit-function-declaration "
fi
# check gfortran version for arg check
let GFOVERSIONGT7=$(expr `${FC} -dumpversion | cut -f1 -d.` \> 7)
if [[ ${GFOVERSIONGT7} == 1 ]]; then
    sixty4_int+=' FCFLAGS=-std=legacy '
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
echo 64ints is $sixty4_int
../configure \
    FC=$MPIF90 CC=$MPICC \
    $sixty4_int \
 ${FORCETARGET} \
	      --disable-shared --enable-static  \
	      --disable-gpu \
--prefix=${NWCHEM_TOP}/src/libext \
SCALAPACK_FCFLAGS="-L${NWCHEM_TOP}/src/libext/lib -lnwc_scalapack  -lnwc_openblas" \
SCALAPACK_LDFLAGS="-L${NWCHEM_TOP}/src/libext/lib -lnwc_scalapack  -lnwc_openblas" \
	      LIBS="-L${NWCHEM_TOP}/src/libext/lib -lnwc_scalapack  -lnwc_openblas"
make VERBOSE=0 V=0 -j3
if [[ "$?" != "0" ]]; then
    echo " "
    echo "Elpa compilation failed"
    echo " "
    exit 1
fi
make install
ln -sf ${NWCHEM_TOP}/src/libext/lib/libelpa.a  ${NWCHEM_TOP}/src/libext/lib/libnwc_elpa.a
ln -sf ${NWCHEM_TOP}/src/libext/include/elpa-${SHORTVERSION}  ${NWCHEM_TOP}/src/libext/include/elpa
