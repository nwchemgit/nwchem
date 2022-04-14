get_cmake_release(){
    UNAME_S=$(uname -s)
    CPU=$(uname -m)
    CMAKE_VER=3.22.2
    echo pwd is `pwd`
    orgdir=`pwd`
    cmake_instdir=$1
    echo "Parameter #1 is $1"
    echo cmake_instdir is $cmake_instdir
    rm -f cmake-${CMAKE_VER}.tar.gz
    if [[ ${UNAME_S} == "Linux" ]] && [[ ${CPU} == "x86_64" ||  ${CPU} == "aarch64" || ${CPU} == "i686" ]] ; then
	cd $cmake_instdir
	if [[ ${CPU} == "i686" ]] ; then
	    CMAKE_CPU="x86_64"
	else
	    CMAKE_CPU=${CPU}
	fi
	CMAKE=`pwd`/cmake-${CMAKE_VER}-linux-${CMAKE_CPU}/bin/cmake
	CMAKE_URL=https://github.com/Kitware/CMake/releases/download/v${CMAKE_VER}/cmake-${CMAKE_VER}-linux-${CMAKE_CPU}.tar.gz
    elif [[ ${UNAME_S} == "Darwin" ]] ; then
	cd $cmake_instdir
	CMAKE=`pwd`/cmake-${CMAKE_VER}-macos-universal/CMake.app/Contents/bin/cmake
	CMAKE_URL=https://github.com/Kitware/CMake/releases/download/v${CMAKE_VER}/cmake-${CMAKE_VER}-macos-universal.tar.gz
    else
	return 1
    fi
    if [ -f ${CMAKE} ]; then
	echo using existing ${CMAKE_VER} Cmake
    else
	curl -L ${CMAKE_URL} -o cmake-${CMAKE_VER}.tar.gz
	tar xzf cmake-${CMAKE_VER}.tar.gz
	CMAKE=`pwd`/cmake-${CMAKE_VER}-linux-${CMAKE_CPU}/bin/cmake
    fi
    cd $orgdir

}

get_cmake_master(){
    CMAKE_COMMIT=09dd52c9d2684e933a3e013abc4f6848cb1befbf
    if [[ -f "cmake-$CMAKE_COMMIT.zip" ]]; then
	echo "using existing"  "cmake-$CMAKE_COMMIT.zip"
    else
	curl -L https://gitlab.kitware.com/cmake/cmake/-/archive/$CMAKE_COMMIT.zip -o cmake-$CMAKE_COMMIT.zip
    fi
    unzip -n -q cmake-$CMAKE_COMMIT.zip
    mkdir -p  cmake-$CMAKE_COMMIT/build
    cd cmake-$CMAKE_COMMIT/build
    if [[ -x "$(command -v cmake)" ]]; then
        cmake -DBUILD_CursesDialog=OFF -DBUILD_TESTING=OFF -DBUILD_QtDialog=OFF -DCMAKE_INSTALL_PREFIX=`pwd`/.. ../
    else
	../bootstrap --parallel=4 --prefix=`pwd`/..
    fi
    make -j4
    make -j4 install
    CMAKE=`pwd`/../bin/cmake
    ${CMAKE} -version
    cd ../..
    return 0
}
