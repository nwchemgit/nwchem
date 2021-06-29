#!/usr/bin/env bash

VERSION=5.1.5

if [[ -f "libxc-${VERSION}.tar.gz" ]]; then
    echo "using existing  libxc-${VERSION}.tar.gz"
else
    echo "downloading libxc-${VERSION}.tar.gz"
    curl -L https://gitlab.com/libxc/libxc/-/archive/${VERSION}/libxc-${VERSION}.tar.gz -o libxc-${VERSION}.tar.gz
fi

tar -xzf libxc-${VERSION}.tar.gz
ln -sf libxc-${VERSION} libxc

if [[  -z "${CC}" ]]; then
    CC=cc
fi
if [[  -z "${FC}" ]]; then
#FC not defined. Look for gfortran
    if [[ ! -x "$(command -v gfortran)" ]]; then
	echo ' '
	echo 'please define FC to compile libxc'
	echo ' '
	exit 1
    else
	echo 'FC not defined, defaulting FC=gfortran'
	FC=gfortran
    fi
fi

cd libxc
mkdir build

cmake -H. -Bbuild -DCMAKE_INSTALL_PREFIX=${NWCHEM_TOP}/src/libext/libxc/install -DCMAKE_C_COMPILER=$CC -DENABLE_FORTRAN=ON -DCMAKE_Fortran_COMPILER=$FC -DDISABLE_KXC=OFF

cd build
make -j2 | tee make.log
make install

