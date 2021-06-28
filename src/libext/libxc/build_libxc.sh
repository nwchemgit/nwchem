#!/usr/bin/env bash

VERSION=5.1.5

if [[ -f "libxc-${VERSION}.tar.gz" ]]; then
    echo "using exisiting  libxc-${VERSION}.tar.gz"
else
    echo "downloading libxc-${VERSION}.tar.gz"
    curl -L https://gitlab.com/libxc/libxc/-/archive/${VERSION}/libxc-${VERSION}.tar.gz -o libxc-${VERSION}.tar.gz
fi

tar -xzf libxc-${VERSION}.tar.gz
ln -s libxc-${VERSION} libxc

cd libxc
mkdir build

cmake -H. -Bbuild -DCMAKE_INSTALL_PREFIX=${NWCHEM_TOP}/src/libext/libxc/install -DCMAKE_C_COMPILER=$CC -DENABLE_FORTRAN=ON -DCMAKE_Fortran_COMPILER=$FC -DDISABLE_KXC=OFF

cd build
make -j2 | tee make.log
make install

