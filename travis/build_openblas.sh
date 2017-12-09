#!/bin/bash 
cd $TRAVIS_BUILD_DIR
git clone -b release-0.2.21 https://github.com/xianyi/OpenBLAS.git
cd OpenBLAS
source $TRAVIS_BUILD_DIR/travis/nwchem.bashrc
if [[ "$USE_64TO32" == "y" ]]; then
  sixty4_int=0
else
  sixty4_int=1
fi
$TRAVIS_BUILD_DIR/travis/sleep_loop.sh make -j3 PREFIX=$TRAVIS_BUILD_DIR INTERFACE64="$sixty4_int" USE_THREAD=0 NO_CBLAS=1 NO_LAPACKE=1 DEBUG=1 NUM_THREADS=1 all
make PREFIX=$TRAVIS_BUILD_DIR install
#tail -4 openblas.log
