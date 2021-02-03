#!/usr/bin/env bash
set -v
arch=`uname -m`
VERSION=0.3.13
if [ -f  OpenBLAS-${VERSION}.tar.gz ]; then
    echo "using existing"  OpenBLAS-${VERSION}.tar.gz
else
    rm -rf OpenBLAS*
    curl -L https://github.com/xianyi/OpenBLAS/archive/v${VERSION}.tar.gz -o OpenBLAS-${VERSION}.tar.gz
fi
# patch for avx2 detection
tar xzf OpenBLAS-${VERSION}.tar.gz
ln -sf OpenBLAS-${VERSION} OpenBLAS
cd OpenBLAS-${VERSION}
# patch for apple clang -fopenmp
patch -p0 < ../clang_omp.patch
if [[  -z "${FORCETARGET}" ]]; then
FORCETARGET=" "
UNAME_S=$(uname -s)
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
if [[ "${GOTAVX2}" == "Y" ]]; then
    echo "forcing Haswell target when AVX2 is available"
    FORCETARGET=" TARGET=HASWELL "
fi
if [[ "${GOTAVX512}" == "Y" ]]; then
    echo "forcing Haswell target on SkyLake"
    FORCETARGET=" TARGET=HASWELL "
fi
fi #FORCETARGET
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
if [ -n "${USE_DYNAMIC_ARCH}" ]; then
    FORCETARGET+="DYNAMIC_ARCH=1 DYNAMIC_OLDER=1"
fi    
if [[ -n ${FC} ]] &&  [[ ${FC} == xlf ]] || [[ ${FC} == xlf_r ]] || [[ ${FC} == xlf90 ]]|| [[ ${FC} == xlf90_r ]]; then
    FORCETARGET+=' CC=gcc FC="xlf -qextname "'
    LAPACK_FPFLAGS_VAL=" -qstrict=ieeefp -O2 -g" 
elif  [[ -n ${FC} ]] && [[ "${FC}" == "flang" ]]; then
    FORCETARGET+=' F_COMPILER=FLANG '
    LAPACK_FPFLAGS_VAL=" -O1 -g -Kieee"
elif  [[ -n ${FC} ]] && [[ "${FC}" == "pgf90" ]] || [[ "${FC}" == "nvfortran" ]]; then
    FORCETARGET+=' F_COMPILER=PGI '
    LAPACK_FPFLAGS_VAL=" -O1 -g -Kieee"
elif  [[ -n ${FC} ]] && [[ "${FC}" == "ifort" ]] || [[ "${FC}" == "ifx" ]]; then
    FORCETARGET+=' F_COMPILER=INTEL '
    LAPACK_FPFLAGS_VAL=" -fp-model source -O2 -g "
else
    LAPACK_FPFLAGS_VAL=" "
fi
#disable threading for ppc64le since it uses OPENMP
echo arch is "$arch"
if [[ "$arch" == "ppc64le" ]]; then
    THREADOPT="0"
else
    THREADOPT="1"
fi

echo make $FORCETARGET  LAPACK_FPFLAGS="$LAPACK_FPFLAGS_VAL"  INTERFACE64="$sixty4_int" BINARY="$binary" NUM_THREADS=128 NO_CBLAS=1 NO_LAPACKE=1 DEBUG=0 USE_THREAD="$THREADOPT"  libs netlib -j4
 make $FORCETARGET  LAPACK_FPFLAGS="$LAPACK_FPFLAGS_VAL"  INTERFACE64="$sixty4_int" BINARY="$binary" NUM_THREADS=128 NO_CBLAS=1 NO_LAPACKE=1 DEBUG=0 USE_THREAD="$THREADOPT" libs netlib -j4

mkdir -p ../../lib
cp libopenblas.a ../../lib/libnwc_openblas.a
#make PREFIX=. install
