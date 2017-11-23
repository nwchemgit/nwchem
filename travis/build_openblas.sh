#!/bin/bash 
cd $TRAVIS_BUILD_DIR
git clone -b release-0.2.21 https://github.com/xianyi/OpenBLAS.git
cd OpenBLAS
make -j3 PREFIX=$TRAVIS_BUILD_DIR  USE_THREAD=0 NO_CBLAS=1 NO_LAPACKE=1 DEBUG=1 NUM_THREADS=1 all >& openblas.log
make PREFIX=$TRAVIS_BUILD_DIR install
tail -4 openblas.log
