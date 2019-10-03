#!/bin/bash -f
# source env. variables
 source $TRAVIS_BUILD_DIR/travis/nwchem.bashrc
# cd $TRAVIS_BUILD_DIR && pwd && cd ..; ln -sf nwchem nwchem-7.0.0 
 ls -lrt $TRAVIS_BUILD_DIR|tail -2
 cd $TRAVIS_BUILD_DIR/src
     make nwchem_config
 if [[ "$USE_64TO32" == "y" ]]; then
     echo " CONVERSION 64_to_32"
     make 64_to_32 >& 6log &
 fi
 env
