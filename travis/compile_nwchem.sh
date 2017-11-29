#!/bin/bash -f
# source env. variables
 source $TRAVIS_BUILD_DIR/travis/nwchem.bashrc
 ls -lrt $NWCHEM_TOP
cd $NWCHEM_TOP/src && ../travis/sleep_loop.sh make  -j3 FDEBUG="-O0 -g" FOPTIMIZE="-O2 -fno-aggressive-loop-optimizations" -j3
tail -2 $NWCHEM_TOP/src/6log
head -2 $NWCHEM_TOP/src/tools/build/config.log
tail -2 $NWCHEM_TOP/src/tools/build/config.log
tail -10 $NWCHEM_TOP/src/make.log
#cd $NWCHEM_TOP/src && make link
