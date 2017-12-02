#!/bin/bash -f
# source env. variables
 source $TRAVIS_BUILD_DIR/travis/nwchem.bashrc
 ls -lrt $NWCHEM_TOP|tail -2
 os=`uname`
 cd $NWCHEM_TOP/src
 if [[ "$os" == "Darwin" ]]; then 
 #    make NWCHEM_MODULES="nwdft driver stepper solvation"  nwchem_config
     make nwchem_config
elif [[ "$os" == "Linux" ]]; then
     make nwchem_config
#     make NWCHEM_MODULES="qmandpw"  nwchem_config
 fi
 make 64_to_32 >& 6log &

