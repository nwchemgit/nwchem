#!/bin/bash -f
# source env. variables
 source $TRAVIS_BUILD_DIR/travis/nwchem.bashrc
ls -lrt $NWCHEM_TOP|tail -3
os=`uname`
if [[ "$NWCHEM_MODULES" == "tce" ]]; then 
#  cd $NWCHEM_TOP/src/tce 
#  head -2 dependencies 
#  rm -f dependencies *amp
#  grep -i dot $NWCHEM_TOP/src/tce/tce_residual_t1.F
#  make include_stamp dependencies >& dep.log 
#  tail -2 dep.log
    export EACCSD=1
    export IPCCSD=1
fi
cd $NWCHEM_TOP/src
FOPT2="-O2 -fno-aggressive-loop-optimizations"
 if [[ "$os" == "Darwin" ]]; then 
   if [[ "$NWCHEM_MODULES" == "tce" ]]; then
     FOPT2="-O1 -fno-aggressive-loop-optimizations"
   fi
     ../travis/sleep_loop.sh make  FDEBUG="-O0 -g" FOPTIMIZE="$FOPT2" -j3
     ../contrib/getmem.nwchem 500
     otool -L ../bin/MACX64/nwchem
#     printenv DYLD_LIBRARY_PATH
#     ls -lrt $DYLD_LIBRARY_PATH
#      tail -120 make.log
 elif [[ "$os" == "Linux" ]]; then
     ../travis/sleep_loop.sh make  FDEBUG="-O0 -g" FOPTIMIZE="$FOPT2" -j3
     cd $NWCHEM_TOP/src/64to32blas 
     make
     cd $NWCHEM_TOP/src
     $NWCHEM_TOP/contrib/getmem.nwchem 500
 fi
#tail -2 $NWCHEM_TOP/src/6log
#head -2 $NWCHEM_TOP/src/tools/build/config.log
#tail -2 $NWCHEM_TOP/src/tools/build/config.log
#tail -10 $NWCHEM_TOP/src/make.log
#grep HAVE_SCA $TRAVIS_BUILD_DIR/src/tools/build/config.h
#cd $NWCHEM_TOP/src && make link
