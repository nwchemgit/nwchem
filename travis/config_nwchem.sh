#!/bin/bash -f
# source env. variables
 source $TRAVIS_BUILD_DIR/travis/nwchem.bashrc
 ls -lrt $NWCHEM_TOP
 cd $NWCHEM_TOP/src && make nwchem_config &&  make 64_to_32 >& 6log &

