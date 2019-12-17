#!/bin/bash 
rm -rf mpich*
VERSION=3.3.1
curl -L http://www.mpich.org/static/downloads/${VERSION}/mpich-${VERSION}.tar.gz -o mpich.tgz
tar xzf mpich.tgz
#patch -p0 < mpistruct.patch
ln -sf mpich-${VERSION} mpich
cd mpich
./configure --prefix=`pwd`/../.. --enable-fast --enable-fortran=all --disable-shared
mkdir -p ../../../lib
make
make install
