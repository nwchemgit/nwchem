#!/bin/bash 
if [ -z "$USE_64TO32"  ] ; then exit 0; fi
os=`uname`
#git clone https://github.com/scibuilder/scalapack.git
#svn co --non-interactive --trust-server-cert https://icl.utk.edu/svn/scalapack-dev/scalapack/trunk/ scalapack
rm -rf scalapack*
curl -LJO http://www.netlib.org/scalapack/scalapack.tgz
VERSION=2.0.2
tar xzf scalapack.tgz
ln -sf scalapack-${VERSION} scalapack
mkdir -p scalapack/build
cd scalapack/build
cmake -Wno-dev ../ -DCMAKE_BUILD_TYPE=RelWithDebInfo -DTEST_SCALAPACK=OFF  -DBUILD_TESTING=OFF -DBUILD_SHARED_LIBS=OFF  -DBLAS_openblas_LIBRARY=$BLASOPT
make V=0 -j3
mkdir -p ../../../lib
cp lib/libscalapack.a ../../../lib
