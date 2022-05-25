#!/usr/bin/env bash
#set -x
source ../libext_utils/cmake.sh

check_tgz() {
    myexit=0
    [ -f $1 ] && gunzip -t $1 > /dev/null && myexit=1
    echo $myexit
}

VERSION=5.2.2
TGZ=libxc-${VERSION}.tar.gz
if [ `check_tgz $TGZ` == 1 ]; then
    echo "using existing $TGZ"
else
    echo "downloading $TGZ"
    curl -L https://gitlab.com/libxc/libxc/-/archive/${VERSION}/libxc-${VERSION}.tar.gz -o $TGZ
    if [ `check_tgz $TGZ` != 1 ]; then
	rm -f libxc-${VERSION}.tar.gz
	curl -L https://github.com/ElectronicStructureLibrary/libxc/archive/refs/tags/${VERSION}.tar.gz -o $TGZ
	if [ `check_tgz $TGZ` != 1 ]; then
	    echo
	    echo libxc download failed
	    echo
	    exit 1
	fi
    fi
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
	cmake_instdir=../libext_utils
	get_cmake_release $cmake_instdir
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
    get_cmake_release  $cmake_instdir
    status=$?
    if [ $status -ne 0 ]; then
	echo cmake required to build scalapack
	echo Please install cmake
	echo define the CMAKE env. variable
	exit 1
    fi
fi

cd libxc
# patch pk09 to avoid compiler  memory problems
patch -p0 -N < ../pk09.patch
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
rm -rf libxc/build
$CMAKE -E env CFLAGS="$cflags" LDFLAGS="$ldflags" FCFLAGS="$fcflags" FFLAGS="$fcflags" \
$CMAKE  -DCMAKE_INSTALL_PREFIX=${NWCHEM_TOP}/src/libext/libxc/install -DCMAKE_C_COMPILER=$CC -DENABLE_FORTRAN=ON -DCMAKE_Fortran_COMPILER=$FC -DDISABLE_KXC=OFF \
-DCMAKE_INSTALL_LIBDIR="lib" -DCMAKE_BUILD_TYPE=Release ..

make -j2 | tee make.log
make install

