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
CMAKE_VER_MAJ=$(${CMAKE} --version|cut -d " " -f 3|head -1|cut -d. -f1)
CMAKE_VER_MIN=$(${CMAKE} --version|cut -d " " -f 3|head -1|cut -d. -f2)
echo CMAKE_VER is ${CMAKE_VER_MAJ} ${CMAKE_VER_MIN}
if ((CMAKE_VER_MAJ < 3)) || (((CMAKE_VER_MAJ > 2) && (CMAKE_VER_MIN < 8))); then
    get_cmake38
    status=$?
    if [ $status -ne 0 ]; then
	echo cmake required to build scalapack
	echo Please install cmake
	echo define the CMAKE env. variable
	exit 1
    fi
fi

cd libxc
mkdir -p build
cd build
if [[ -z "${NWCHEM_TOP}" ]]; then
#    DIRQA=`dirname "$0"`
    MYPWD=`pwd`
#    echo DIRQ $DIRQA
    NWCHEM_TOP=`echo ${MYPWD} | sed -e 's/\/src.*//' `
fi
echo @@ NWCHEM_TOP @@ $NWCHEM_TOP
if [[ ! -z $NWCHEM_TARGET && $NWCHEM_TARGET == "LINUX" ]] ; then
    if [[ `uname -m` == x86_64 ]] ; then
	ldflags=-m32
	cflags=-m32
	fcflags=-m32
    fi
else
    ldflags=" "
    cflags=" "
    fcflags=" "
fi
CMAKE_EXE_LINKER_FLAGS=
$CMAKE  -DCMAKE_INSTALL_PREFIX=${NWCHEM_TOP}/src/libext/libxc/install -DCMAKE_C_COMPILER=$CC -DENABLE_FORTRAN=ON -DCMAKE_Fortran_COMPILER=$FC -DDISABLE_KXC=OFF \
-DCMAKE_EXE_LINKER_FLAGS=$ldflags  -DCMAKE_Fortran_FLAGS=$fcflags -DCMAKE_C_FLAGS=$cflags \
-DCMAKE_INSTALL_LIBDIR="lib" -DCMAKE_BUILD_TYPE=Release ..


make -j2 | tee make.log
make install

