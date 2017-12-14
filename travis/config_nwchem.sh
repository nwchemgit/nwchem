#!/bin/bash -f
# source env. variables
 source $TRAVIS_BUILD_DIR/travis/nwchem.bashrc
 ls -lrt $NWCHEM_TOP|tail -2
 cd $NWCHEM_TOP/src
     make nwchem_config
 if [[ "$USE_64TO32" == "y" ]]; then
     echo " CONVERSION 64_to_32"
     make 64_to_32 >& 6log &
 fi
 env
