#!/usr/bin/env bash
source ../libext_utils/get_tgz.sh
rm -rf mpich  mpich-?.?.?
#VERSION=3.4.2
#VERSION=4.0.2
VERSION=4.3.0
#curl -L http://www.mpich.org/static/downloads/${VERSION}/mpich-${VERSION}.tar.gz -o mpich.tgz
#curl -L https://github.com/pmodels/mpich/releases/download/v${VERSION}/mpich-${VERSION}.tar.gz -o mpich.tgz
get_tgz https://github.com/pmodels/mpich/releases/download/v${VERSION}/mpich-${VERSION}.tar.gz  mpich.tgz
tar xzf mpich.tgz
ln -sf mpich-${VERSION} mpich
cd mpich
if [[  -z "${FC}" ]]; then
    export FC=gfortran
fi
if [[  -z "${NWCHEM_TOP}" ]]; then
    dir3=$(dirname `pwd`)
    dir2=$(dirname "$dir3")
    NWCHEM_TOP=$(dirname "$dir2")
fi
echo FC is $FC
FC_EXTRA=$(${NWCHEM_TOP}/src/config/strip_compiler.sh ${FC})
echo FC_EXTRA is $FC_EXTRA
if  [[ ${FC_EXTRA} == nvfortran ]] ; then
    export FFLAGS+=" -fPIC "
    export FCFLAGS+=" -fPIC "
fi
#mpich crashes when F90 and F90FLAGS are set
unset F90
unset F90FLAGS
echo 'using FFLAGS=' $FFLAGS
if [[  -z "${BUILD_ELPA}" ]]; then
    SHARED_FLAGS="--disable-shared"
else
    SHARED_FLAGS="--enable-shared"
    if  [[ ${FC_EXTRA} == amdflang ]] ; then
	echo ${FC_EXTRA} not compatible with shared libraries
	exit 1
    fi
fi
#cross compilation
if [ -x "$(command -v xx-info)" ]; then
    SHARED_FLAGS+=" --host=$(xx-info triple) "
    SHARED_FLAGS+=" --build=$(uname -m)-linux-gnu "
    # 32bit or 64bit arch?
    arch=$(xx-info march)
    echo "mpich arch is " $arch
    if echo $arch |grep -q 32 ; then
	SHARED_FLAGS+=" --with-cross=../cross32.txt "
    else
	SHARED_FLAGS+=" --with-cross=../cross.txt "
    fi
fi
echo SHARED_FLAGS is $SHARED_FLAGS
if [ $(uname -s) == "Darwin" ]; then
    if pkg-config hwloc --exists; then
	HWLOC_FLAGS=" --with-hwloc=$(pkg-config hwloc  --variable=prefix) "
    fi
fi
echo HWLOC_FLAGS is $HWLOC_FLAGS
./configure --prefix=`pwd`/../.. --enable-fortran=all $SHARED_FLAGS  --disable-cxx --disable-romio  --enable-silent-rules  --enable-fortran=all \
	    --without-slurm  --with-pm=hydra \
	    --without-cuda --without-ze --without-hip --without-hcoll \
	    --with-device=ch3 \
	    $HWLOC_FLAGS
#./configure --prefix=`pwd`/../.. --enable-fortran=all $SHARED_FLAGS  --disable-cxx --enable-romio --with-pm=gforker --with-device=ch3:nemesis --disable-cuda --disable-opencl --enable-silent-rules  --enable-fortran=all
if [[ "$?" != "0" ]]; then
    cat config.log
    echo "MPICH configuration failed"
    exit 1
fi
mkdir -p ../../../lib
echo
echo mpich compilation in progress
echo output redirected to libext/mpich/mpich/mpich.log
echo
make -j3 >& mpich.log
if [[ "$?" != "0" ]]; then
    tail -500 mpich.log
    echo " "
    echo "MPICH compilation failed"
    echo " "
    exit 1
fi

make install >& mpich_install.log
${NWCHEM_TOP}/src/tools/check_libmpi.sh ${NWCHEM_TOP}
