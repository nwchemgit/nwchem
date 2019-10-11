#!/bin/bash -f
# source env. variables
 source $TRAVIS_BUILD_DIR/travis/nwchem.bashrc
ls -lrt $TRAVIS_BUILD_DIR|tail -3
os=`uname`
arch=`uname -m`
if [[ "$NWCHEM_MODULES" == "tce" ]]; then 
    export EACCSD=1
    export IPCCSD=1
fi
cd $TRAVIS_BUILD_DIR/src
if [[ "$arch" == "aarch64" ]]; then 
    FOPT2="-O1 -fno-aggressive-loop-optimizations"
else    
    FOPT2="-O2 -fno-aggressive-loop-optimizations"
fi    
 if [[ "$os" == "Darwin" ]]; then 
   if [[ "$NWCHEM_MODULES" == "tce" ]]; then
     FOPT2="-O1 -fno-aggressive-loop-optimizations"
   fi
     ../travis/sleep_loop.sh make  FDEBUG="-O0 -g" FOPTIMIZE="$FOPT2" -j3
     cd $TRAVIS_BUILD_DIR/src/64to32blas 
     make
     cd $TRAVIS_BUILD_DIR/src
     ../contrib/getmem.nwchem 500
     otool -L ../bin/MACX64/nwchem
#     printenv DYLD_LIBRARY_PATH
#     ls -lrt $DYLD_LIBRARY_PATH
#      tail -120 make.log
 elif [[ "$os" == "Linux" ]]; then
     ../travis/sleep_loop.sh make  FDEBUG="-O0 -g" FOPTIMIZE="$FOPT2" -j3
     cd $TRAVIS_BUILD_DIR/src/64to32blas 
     make
     cd $TRAVIS_BUILD_DIR/src
     $TRAVIS_BUILD_DIR/contrib/getmem.nwchem 500
 fi
#tail -2 $NWCHEM_TOP/src/6log
#head -2 $NWCHEM_TOP/src/tools/build/config.log
#tail -2 $NWCHEM_TOP/src/tools/build/config.log
#tail -10 $NWCHEM_TOP/src/make.log
#grep HAVE_SCA $TRAVIS_BUILD_DIR/src/tools/build/config.h
#cd $NWCHEM_TOP/src && make link
