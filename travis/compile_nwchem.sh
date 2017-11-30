#!/bin/bash -f
# source env. variables
 source $TRAVIS_BUILD_DIR/travis/nwchem.bashrc
ls -lrt $NWCHEM_TOP|tail -3
os=`uname`
cd $NWCHEM_TOP/src
 if [[ "$os" == "Darwin" ]]; then 
     ../travis/sleep_loop.sh make  -j3 FDEBUG="-O0 -g" FOPTIMIZE="-O1 -fno-aggressive-loop-optimizations" -j3
 elif [[ "$os" == "Linux" ]]; then
     ../travis/sleep_loop.sh make  -j3 FDEBUG="-O0 -g" FOPTIMIZE="-O2 -fno-aggressive-loop-optimizations" -j3
 fi
tail -2 $NWCHEM_TOP/src/6log
head -2 $NWCHEM_TOP/src/tools/build/config.log
tail -2 $NWCHEM_TOP/src/tools/build/config.log
tail -10 $NWCHEM_TOP/src/make.log
grep HAVE_SCA $TRAVIS_BUILD_DIR/src/tools/build/config.h
#cd $NWCHEM_TOP/src && make link
