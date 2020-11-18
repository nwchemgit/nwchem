#!/bin/bash 
VERSION=0.3.12
if [ -f  OpenBLAS-${VERSION}.tar.gz ]; then
    echo "using existing"  OpenBLAS-${VERSION}.tar.gz
else
    rm -rf OpenBLAS*
    curl -L https://github.com/xianyi/OpenBLAS/archive/v${VERSION}.tar.gz -o OpenBLAS-${VERSION}.tar.gz
fi
tar xzf OpenBLAS-${VERSION}.tar.gz
ln -sf OpenBLAS-${VERSION} OpenBLAS
cd OpenBLAS-${VERSION}
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
if [[ -n ${FC} ]] &&  [[ ${FC} == xlf ]] || [[ ${FC} == xlf_r ]] || [[ ${FC} == xlf90 ]]|| [[ ${FC} == xlf90_r ]]; then
 make CC=gcc FC="xlf -qextname"  INTERFACE64="$sixty4_int" BINARY="$binary" USE_THREAD=0 NO_CBLAS=1 NO_LAPACKE=1 DEBUG=0 NUM_THREADS=1 LAPACK_FPFLAGS="-qstrict=ieeefp -O2 -g" libs netlib
elif  [[ -n ${FC} ]] && [[ "${FC}" == "flang" ]]; then
 make $FORCETARGET LAPACK_FPFLAGS="-O1 -g -Kieee" INTERFACE64="$sixty4_int" BINARY="$binary" USE_THREAD=0 NO_CBLAS=1 NO_LAPACKE=1 DEBUG=0 NUM_THREADS=1 libs netlib 
elif  [[ -n ${FC} ]] && [[ "${FC}" == "ifort" ]] || [[ "${FC}" == "ifx" ]]; then
 make $FORCETARGET  F_COMPILER=INTEL LAPACK_FPFLAGS="-fp-model source -O2 -g" INTERFACE64="$sixty4_int" BINARY="$binary" USE_THREAD=0 NO_CBLAS=1 NO_LAPACKE=1 DEBUG=0 NUM_THREADS=1 libs netlib 
else
 make $FORCETARGET  INTERFACE64="$sixty4_int" BINARY="$binary" USE_THREAD=0 NO_CBLAS=1 NO_LAPACKE=1 DEBUG=0 NUM_THREADS=1  libs netlib -j4
fi
mkdir -p ../../lib
cp libopenblas.a ../../lib/libnwc_openblas.a
#make PREFIX=. install
