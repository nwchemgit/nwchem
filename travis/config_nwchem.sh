#!/bin/bash
# source env. variables
 source $TRAVIS_BUILD_DIR/travis/nwchem.bashrc
if [[ ! -z "$TARBALL" ]]; then
    cd $TRAVIS_BUILD_DIR/..
    pwd
    mv nwchem nwchem.git
    curl -LJ $TARBALL -o nwchem.tar.bz2
    tar xjf nwchem.tar.bz2
    ln -sf nwchem-7.0.0 nwchem
    mv nwchem/travis nwchem/travis.tarball
    cp -r nwchem.git/travis nwchem/.
fi
 ls -lrt $TRAVIS_BUILD_DIR|tail -2
 cd $TRAVIS_BUILD_DIR/src
     make nwchem_config
if [[ -z "$TARBALL" ]]; then
 if [[ "$USE_64TO32" == "y" ]]; then
     echo " CONVERSION 64_to_32"
     make 64_to_32 >& 6log &
 fi
fi
 env
