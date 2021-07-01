#!/usr/bin/env bash
get_cmake38(){
	UNAME_S=$(uname -s)
	if [[ ${UNAME_S} == "Linux" ]] || [[ ${UNAME_S} == "Darwin" ]] && [[ $(uname -m) == "x86_64" ]] ; then
	    CMAKE_VER=3.16.8
	    rm -f cmake-${CMAKE_VER}-${UNAME_S}-x86_64.tar.gz
	    curl -L https://github.com/Kitware/CMake/releases/download/v${CMAKE_VER}/cmake-${CMAKE_VER}-${UNAME_S}-x86_64.tar.gz -o cmake-${CMAKE_VER}-${UNAME_S}-x86_64.tar.gz
	    tar xzf cmake-${CMAKE_VER}-${UNAME_S}-x86_64.tar.gz
	    if [[ ${UNAME_S} == "Darwin" ]] ;then
		CMAKE=`pwd`/cmake-${CMAKE_VER}-${UNAME_S}-x86_64/CMake.app/Contents/bin/cmake
	    else
		CMAKE=`pwd`/cmake-${CMAKE_VER}-${UNAME_S}-x86_64/bin/cmake
	    fi
	    return 0
	else
	    return 1
	fi

}

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
mkdir -p build

if [[ -z "${CMAKE}" ]]; then
    #look for cmake
    if [[ -z "$(command -v cmake)" ]]; then
	get_cmake38
	status=$?
	if [ $status -ne 0 ]; then
	    echo cmake required to build libxc
	    echo Please install cmake
	    echo define the CMAKE env. variable
	    exit 1
	fi
    else
	CMAKE=cmake
    fi
fi
cd build
$CMAKE -H.  -DCMAKE_INSTALL_PREFIX=${NWCHEM_TOP}/src/libext/libxc/install -DCMAKE_C_COMPILER=$CC -DENABLE_FORTRAN=ON -DCMAKE_Fortran_COMPILER=$FC -DDISABLE_KXC=OFF ..

make -j2 | tee make.log
make install

