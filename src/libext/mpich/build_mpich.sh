#!/usr/bin/env bash
rm -rf mpich*
#VERSION=3.4.2
VERSION=3.3.2
curl -L http://www.mpich.org/static/downloads/${VERSION}/mpich-${VERSION}.tar.gz -o mpich.tgz
tar xzf mpich.tgz
ln -sf mpich-${VERSION} mpich
cd mpich
GNUMAJOR=`$FC -dM -E - < /dev/null 2> /dev/null | grep __GNUC__ |cut -c18-`
if [[ $GNUMAJOR -ge 10  ]]; then
    export FFLAGS=-std=legacy
fi
#mpich crashes when F90 and F90FLAGS are set
unset F90
unset F90FLAGS
echo 'using FFLAGS=' $FFLAGS
./configure --prefix=`pwd`/../.. --enable-fortran=all --disable-shared --disable-cxx --enable-romio --with-pm=gforker --with-device=ch3:nemesis --disable-cuda --disable-opencl --enable-silent-rules  --enable-fortran=all
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
