#!/usr/bin/env bash
#set -v
arch=`uname -m`
VERSION=0.3.21
#COMMIT=974acb39ff86121a5a94be4853f58bd728b56b81
BRANCH=develop
if [ -f  OpenBLAS-${VERSION}.tar.gz ]; then
    echo "using existing"  OpenBLAS-${VERSION}.tar.gz
else
    rm -rf OpenBLAS*
    tries=1 ; until [ "$tries" -ge 6 ] ; do
		  if [ "$tries" -gt 1 ]; then sleep 9; echo attempt no.  $tries ; fi
		  curl -L https://github.com/xianyi/OpenBLAS/archive/v${VERSION}.tar.gz -o OpenBLAS-${VERSION}.tar.gz ;
		  # check tar.gz integrity
		  gzip -t OpenBLAS-${VERSION}.tar.gz >&  /dev/null
		  if [ $? -eq 0 ]; then break ;  fi
		  tries=$((tries+1)) ;  done
fi
gzip -t OpenBLAS-${VERSION}.tar.gz >&  /dev/null
if [ $? -ne 0 ]; then echo  "openBLAS tarball not ready"; rm -f OpenBLAS-${VERSION}.tar.gz; exit 1 ; fi
tar xzf OpenBLAS-${VERSION}.tar.gz
ln -sf OpenBLAS-${VERSION} OpenBLAS
cd OpenBLAS
# patch for apple clang -fopenmp
patch -p0 -s -N < ../makesys.patch
#patch -p0 -s -N < ../icc_avx512.patch
# patch for pgi/nvfortran missing -march=armv8
patch -p0 -s -N < ../arm64_fopt.patch
patch -p1 -s -N < ../9402df5604e69f86f58953e3883f33f98c930baf.patch
if [[  -z "${FORCETARGET}" ]]; then
FORCETARGET=" "
UNAME_S=$(uname -s)
if [[ ${UNAME_S} == Linux ]]; then
    CPU_FLAGS=$(cat /proc/cpuinfo | grep flags |tail -n 1)
    CPU_FLAGS_2=$(cat /proc/cpuinfo | grep flags |tail -n 1)
elif [[ ${UNAME_S} == Darwin ]]; then
    CPU_FLAGS=$(/usr/sbin/sysctl -n machdep.cpu.features)
    if [[ "$arch" == "x86_64" ]]; then
#	CPU_FLAGS_2=$(/usr/sbin/sysctl -n machdep.cpu.leaf7_features)
	if [[ $(/usr/sbin/sysctl -n hw.optional.avx2_0) == 1 ]]; then
	    echo got AVX2
	    CPU_FLAGS_2="AVX2"
	fi
    fi
fi
  GOTSSE2=$(echo ${CPU_FLAGS}   | tr  'A-Z' 'a-z'| awk ' /sse2/   {print "Y"}')
   GOTAVX=$(echo ${CPU_FLAGS}   | tr  'A-Z' 'a-z'| awk ' /avx/    {print "Y"}')
  GOTAVX2=$(echo ${CPU_FLAGS_2} | tr  'A-Z' 'a-z'| awk ' /avx2/   {print "Y"}')
GOTAVX512=$(echo ${CPU_FLAGS}   | tr  'A-Z' 'a-z'| awk ' /avx512f/{print "Y"}')
GOTCLZERO=$(echo ${CPU_FLAGS}   | tr  'A-Z' 'a-z'| awk ' /clzero/{print "Y"}')
if [[ "${GOTAVX2}" == "Y" ]]; then
    echo "forcing Haswell target when AVX2 is available"
    FORCETARGET=" TARGET=HASWELL "
fi
if [[ "${GOTCLZERO}" == "Y" ]]; then
    echo "forcing Zen target when CLZERO is available"
    FORCETARGET=" TARGET=ZEN "
fi
if [[ "${GOTAVX512}" == "Y" ]]; then
    echo "forcing Haswell target on SkyLake"
    FORCETARGET=" TARGET=HASWELL "
fi
fi #FORCETARGET
if [[ "$arch" == "riscv64" ]]; then
	 FORCETARGET=" TARGET=RISCV64_GENERIC "
fi
if [[  -z "${BLAS_SIZE}" ]]; then
   BLAS_SIZE=8
fi
if [[ ${BLAS_SIZE} == 8 ]]; then
  sixty4_int=1
else
  sixty4_int=0
fi
if [[ "${NWCHEM_TARGET}" == "LINUX" ]]; then
  binary=32
  sixty4_int=0
else
  binary=64
fi
if [[ -n "${USE_DYNAMIC_ARCH}" ]] || [[ "${USE_HWOPT}" == "n" ]]; then
    if [[ "$arch" == "x86_64" ]]; then
	FORCETARGET+="DYNAMIC_ARCH=1 DYNAMIC_OLDER=1"
    fi
fi
#cray ftn wrapper
if [[ ${FC} == ftn ]]; then
    FCORG=ftn
    CRAY_ACCEL_TARGET_ORG=$CRAY_ACCEL_TARGET
    unset CRAY_ACCEL_TARGET
    if [[ ${PE_ENV} == PGI ]]; then
          FC=pgf90
#          _CC=pgcc
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
        FORCETARGET+=' CC=gcc '
	CC=gcc
	unset CPATH
    fi
    if [[ ${PE_ENV} == CRAY ]]; then
#	echo ' '
#	echo 'openblas installation not ready for crayftn '
#	echo ' '
	if ! [ -x "$(command -v gfortran)" ]; then
	    echo " please load the gcc module (not prgenv)"
	    echo " by executing"
	    echo "     module load gcc "
	    echo " "
	    exit 1
	fi
	FC=gfortran
	CCORG=${CC}
	CC=clang
	export PATH=/opt/cray/pe/cce/default/cce-clang/x86_64/bin:$PATH
        FORCETARGET+=' FC=gfortran CC=clang '
#	exit 1
#        exit 1
    fi
fi
if [[ -n ${FC} ]] &&  [[ ${FC} == xlf ]] || [[ ${FC} == xlf_r ]] || [[ ${FC} == xlf90 ]]|| [[ ${FC} == xlf90_r ]]; then
    FORCETARGET+=" CC=gcc "
    _FC=xlf
    LAPACK_FPFLAGS_VAL=" -qstrict=ieeefp -O2 -g" 
elif  [[ -n ${FC} ]] && [[ "${FC}" == "flang" ]] || [[ "${FC}" == "amdflang" ]]; then
    FORCETARGET+=' F_COMPILER=FLANG '
    LAPACK_FPFLAGS_VAL=" -O1 -g -Kieee"
elif  [[ -n ${FC} ]] && [[ "${FC}" == "pgf90" ]] || [[ "${FC}" == "nvfortran" ]]; then
    FORCETARGET+=' F_COMPILER=PGI '
  if  [[ "${FC}" == "nvfortran" ]]; then
    FORCETARGET+=' NEWPGI2=1 '
  fi
#    if [[ "$arch" == "aarch64" ]]; then
	LAPACK_FLAGS="-O2  -Mrecursive -Kieee -fPIC"
        if [[ ${BLAS_SIZE} == 8 ]]; then
	    LAPACK_FLAGS+=" -i8"
	fi
#    fi
    LAPACK_FPFLAGS_VAL=" -O1 -g -Kieee"
elif  [[ -n ${FC} ]] && [[ "${FC}" == "ifort" ]] || [[ "${FC}" == "ifx" ]]; then
    FORCETARGET+=' F_COMPILER=INTEL '
    LAPACK_FPFLAGS_VAL=" -fp-model source -O2 -g "
else
	#assuming gfortran
    FORCETARGET+=' F_COMPILER=GFORTRAN '
    LAPACK_FPFLAGS_VAL=" "
        if [[ ${BLAS_SIZE} == 8 ]]; then
           LAPACK_FPFLAGS_VAL+=" -fdefault-integer-8"
       fi
fi
if  [[ -n ${CC} ]] && [[ "${CC}" == "amdclang" ]]; then
    let VERSIONEQ15=$(expr `${CC} -dM -E - < /dev/null 2> /dev/null|egrep 15|grep __clang_major__ |cut  -d ' ' -f 3 ` \= 15)
    if [[ ${VERSIONEQ15} == 1 ]]; then
       echo "amdclang 15 buggy. reduced optimization to O1"
       FORCETARGET+=' COMMON_OPT=-O1'
    fi
fi
if [[   -z "${FC}" ]]; then
    FC=gfortran
fi
if [[   -z "${CC}" ]]; then
    CC=cc
fi
if [[ `${CC} -dM -E - < /dev/null 2> /dev/null | grep -c GNU` > 0 ]] ; then
    let GCCVERSIONGT5=$(expr `${CC} -dumpversion | cut -f1 -d.` \> 5)
fi
# check gcc version for skylake
if [[ "$FORCETARGET" == *"SKYLAKEX"* ]]; then
    if [[ ${GCCVERSIONGT5} != 1 ]]; then
	echo
	echo you have gcc version $(${CC} -dumpversion | cut -f1 -d.)
	echo gcc version 6 and later needed for skylake
	echo
	exit 1
    fi
fi
#this fixes avx512 detection for icc
if [[ "${CC}" == "icc" ]]; then
    FORCETARGET+=HOSTCC=\"icc -xhost\"
fi

#disable threading for ppc64le since it uses OPENMP
echo arch is "$arch"
if [[ "$arch" == "ppc64le" ||  "$arch" == "riscv64" ]]; then
if [[ ${GCCVERSIONGT5} != 1 ]]; then
       echo
       echo gcc version 6 and later needed for ppc64le
       echo please specify CC=gcc
       echo where the gcc version is 6 or later
       echo
       exit 1
fi
    THREADOPT="0"
    MYNTS="1"
else
    THREADOPT="1"
    MYNTS="128"
fi
# cross compilation
GOTMINGW64=$("$CC" -dM -E - </dev/null 2> /dev/null |grep MINGW64|cut -c21)
if [[ "${GOTMINGW64}" == "1" ]]; then
    FORCETARGET+=HOSTCC=\"gcc\"
    THREADOPT="0"
    MYNTS="1"
fi

#we want openblas to use pthreads and not openmp.
#but NWChem and OpenBLAS both use USE_OPENMP
#disable USE_OPENMP if set and re-enable it later
if [[  ! -z "${USE_OPENMP}" ]]; then
    unset USE_OPENMP
    NWCHEM_USE_OPENMP=1
fi
GOTFREEBSD=$(uname -o 2>&1|awk ' /FreeBSD/ {print "1";exit}')
MYMAKE=make
MAKEJ=" -j4"
if [[  "${GOTFREEBSD}" == 1 ]]; then
MAKEJ=" "
MYMAKE=gmake
fi
echo FC is $FC
echo $MYMAKE FC=$FC $FORCETARGET LAPACK_FPFLAGS=$LAPACK_FPFLAGS_VAL  INTERFACE64=$sixty4_int BINARY=$binary NUM_THREADS=$MYNTS NO_CBLAS=1 NO_LAPACKE=1 DEBUG=0 USE_THREAD=$THREADOPT  libs netlib $MAKEJ
echo
echo OpenBLAS compilation in progress
echo output redirected to libext/openblas/OpenBLAS/openblas.log
echo
if [[ ${_FC} == xlf ]]; then
 $MYMAKE FC="xlf -qextname" $FORCETARGET  LAPACK_FPFLAGS="$LAPACK_FPFLAGS_VAL"  INTERFACE64="$sixty4_int" BINARY="$binary" NUM_THREADS=$MYNTS NO_CBLAS=1 NO_LAPACKE=1 DEBUG=0 USE_THREAD="$THREADOPT" libs netlib $MAKEJ >& openblas.log
else
 $MYMAKE FC=$FC $FORCETARGET  LAPACK_FPFLAGS="$LAPACK_FPFLAGS_VAL"  INTERFACE64="$sixty4_int" BINARY="$binary" NUM_THREADS=128 NO_CBLAS=1 NO_LAPACKE=1 DEBUG=0 USE_THREAD="$THREADOPT"  libs netlib $MAKEJ >& openblas.log
fi
if [[ "$?" != "0" ]]; then
    tail -500 openblas.log
    echo " "
    echo "OpenBLAS compilation failed"
    echo " "
    exit 1
fi

mkdir -p ../../lib
if [[ $(uname -s) == "Linux" ]]; then
    strip --strip-debug libopenblas*-*.a
fi
cp libopenblas.a ../../lib/libnwc_openblas.a
#make PREFIX=. install
if [[  ! -z "${NWCHEM_USE_OPENMP}" ]]; then
    export USE_OPENMP=1
fi
if [[ -n ${FCORG} ]]; then
    FC=${FCORG}
fi
if [[ -n ${CCORG} ]]; then
    CC=${CCORG}
fi
if [[ -n ${CRAY_ACCEL_TARGET_ORG} ]]; then
    CRAY_ACCEL_TARGET=$CRAY_ACCEL_TARGET_ORG
fi
