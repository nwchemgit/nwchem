#!/usr/bin/env bash
rm -rf mpich*
#VERSION=3.4.2
VERSION=4.0.2
#curl -L http://www.mpich.org/static/downloads/${VERSION}/mpich-${VERSION}.tar.gz -o mpich.tgz
curl -L https://github.com/pmodels/mpich/releases/download/v${VERSION}/mpich-${VERSION}.tar.gz -o mpich.tgz
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
if  [[ ${FC_EXTRA} == gfortran ]] ; then
    GNUMAJOR=`$FC -dM -E - < /dev/null 2> /dev/null | grep __GNUC__ |cut -c18-`
    echo GNUMAJOR is $GNUMAJOR
    if [[ $GNUMAJOR -ge 10  ]]; then
	export FFLAGS=" -fallow-argument-mismatch "
        export FCFLAGS=" -fallow-argument-mismatch "
    fi
fi
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
./configure --prefix=`pwd`/../.. --enable-fortran=all $SHARED_FLAGS  --disable-cxx --enable-romio --with-pm=gforker --with-device=ch3:nemesis --disable-cuda --disable-opencl --enable-silent-rules  --enable-fortran=all
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
