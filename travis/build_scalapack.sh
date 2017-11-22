#!/bin/bash 
cd $TRAVIS_BUILD_DIR
git clone https://github.com/scibuilder/scalapack.git
mkdir -p scalapack/build
cd scalapack/build
cmake ../ -DCMAKE_BUILD_TYPE=RelWithDebInfo -DTEST_SCALAPACK=OFF -DBUILD_TESTING=OFF -DBUILD_SHARED_LIBS=ON  -DBLAS_blas_LIBRARY="-L$TRAVIS_BUILD_DIR/lib -lopenblas" -DLAPACK_lapack_LIBRARY="-L$TRAVIS_BUILD_DIR/lib -lopenblas"  -DCMAKE_INSTALL_PREFIX=$TRAVIS_BUILD_DIR
make V=0 -j3 >& scalap.log
make install
tail -20 scalap.log
ls -lrt $TRAVIS_BUILD_DIR/lib 
