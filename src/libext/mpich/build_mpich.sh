#!/bin/bash 
rm -rf mpich*
VERSION=3.4.2
curl -L http://www.mpich.org/static/downloads/${VERSION}/mpich-${VERSION}.tar.gz -o mpich.tgz
tar xzf mpich.tgz
#patch -p0 < mpistruct.patch
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
./configure --prefix=`pwd`/../.. --enable-fortran=all --disable-shared --disable-cxx --disable-romio --with-pm=gforker --with-device=ch3:nemesis --disable-cuda --disable-opencl
mkdir -p ../../../lib
make -j3
make install
