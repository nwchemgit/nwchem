#!/bin/bash 
rm -rf mpich*
VERSION=3.3.1
curl -L http://www.mpich.org/static/downloads/${VERSION}/mpich-${VERSION}.tar.gz -o mpich.tgz
tar xzf mpich.tgz
#patch -p0 < mpistruct.patch
ln -sf mpich-${VERSION} mpich
cd mpich
GNUMAJOR=`$FC -dM -E - < /dev/null 2> /dev/null | grep __GNUC__ |cut -c18-`
if [[ $GNUMAJOR -ge 10  ]]; then
    export FFLAGS=-std=legacy
fi
echo 'using FFLAGS=' $FFLAGS
./configure --prefix=`pwd`/../.. --enable-fast --enable-fortran=all --disable-shared
mkdir -p ../../../lib
make
make install
