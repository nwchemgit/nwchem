#!/usr/bin/env bash
# source env. variables
if [[ -z "$TRAVIS_BUILD_DIR" ]] ; then
    TRAVIS_BUILD_DIR=$(pwd)
fi
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
if [[ ! -z "$USE_64TO32"  ]]; then
    echo " CONVERSION 64_to_32"
    os=`uname`
    echo output redirected to 64.log
    if [[ "$os" == "Darwin" ]]; then
        make 64_to_32  >& 64log 
    else
        make 64_to_32 CONVERT_ALL=y >& 64log
    fi
    echo @@@@ tail 64log
    tail -30 64log
    echo @@@
fi
env
