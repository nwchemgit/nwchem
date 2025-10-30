export CMAKE_VERSION_REQUIRED=3.26.0
CMAKE_VER_REQ_MAJ=$(echo ${CMAKE_VERSION_REQUIRED}|cut -d . -f 1)
CMAKE_VER_REQ_MIN=$(echo ${CMAKE_VERSION_REQUIRED}|cut -d . -f 2)

get_cmake_release(){
    UNAME_S=$(uname -s)
    CPU=$(uname -m)
    CMAKE_VER=${CMAKE_VERSION_REQUIRED}
    orgdir=`pwd`
    cmake_instdir=$1
    echo "Parameter #1 is $1"
    echo cmake_instdir is $cmake_instdir
    rm -f cmake-${CMAKE_VER}.tar.gz
    if [[ ${UNAME_S} == "Linux" ]]; then
	if [[ ${CPU} == "x86_64" ||  ${CPU} == "aarch64" || ${CPU} == "i686" ]] ; then
	    cd $cmake_instdir
	    if [[ ${CPU} == "i686" ]] ; then
		CMAKE_CPU="x86_64"
	    else
		CMAKE_CPU=${CPU}
	    fi
	    export CMAKE=`pwd`/cmake-${CMAKE_VER}-linux-${CMAKE_CPU}/bin/cmake
	    export PATH=`pwd`/cmake-${CMAKE_VER}-linux-${CMAKE_CPU}/bin:$PATH
	    CMAKE_URL=https://github.com/Kitware/CMake/releases/download/v${CMAKE_VER}/cmake-${CMAKE_VER}-linux-${CMAKE_CPU}.tar.gz
	else
	        get_cmake_master
	fi
    elif [[ ${UNAME_S} == "Darwin" ]] ; then
	cd $cmake_instdir
	export CMAKE=`pwd`/cmake-${CMAKE_VER}-macos-universal/CMake.app/Contents/bin/cmake
	export PATH=`pwd`/cmake-${CMAKE_VER}-macos-universal/CMake.app/Contents/bin:$PATH
	CMAKE_URL=https://github.com/Kitware/CMake/releases/download/v${CMAKE_VER}/cmake-${CMAKE_VER}-macos-universal.tar.gz
    else
	return 1
    fi
    if [ -f ${CMAKE} ]; then
	echo using existing ${CMAKE_VER} Cmake 
    else
	curl -L ${CMAKE_URL} -o cmake-${CMAKE_VER}.tar.gz
	tar xzf cmake-${CMAKE_VER}.tar.gz
	if [[ ${UNAME_S} == "Darwin" ]] ; then
	    export CMAKE=`pwd`/cmake-${CMAKE_VER}-macos-universal/CMake.app/Contents/bin/cmake
	    export PATH=`pwd`/cmake-${CMAKE_VER}-macos-universal/CMake.app/Contents/bin:$PATH
	else
	    export CMAKE=`pwd`/cmake-${CMAKE_VER}-linux-${CMAKE_CPU}/bin/cmake
	    export PATH=`pwd`/cmake-${CMAKE_VER}-linux-${CMAKE_CPU}/bin:$PATH
	fi
    fi
    cd $orgdir

}

get_cmake_master(){
    CMAKE_COMMIT=v3.26.6
    if [[ -f "cmake-$CMAKE_COMMIT.zip" ]]; then
	echo "using existing"  "cmake-$CMAKE_COMMIT.zip" >> /tmp/cmake.log
    else
	curl -L https://github.com/Kitware/CMake/archive/refs/tags/$CMAKE_COMMIT.zip -o cmake-$CMAKE_COMMIT.zip
    fi
    rm -rf CMake*
    unzip -n -q cmake-$CMAKE_COMMIT.zip
    ln -sf CMake* cmake-$CMAKE_COMMIT
    mkdir -p  cmake-$CMAKE_COMMIT/build
    cd cmake-$CMAKE_COMMIT/build
    if [[ -x "$(command -v cmake)" ]]; then
        cmake -DBUILD_CursesDialog=OFF -DBUILD_TESTING=OFF -DBUILD_QtDialog=OFF -DCMAKE_INSTALL_PREFIX=`pwd`/.. ../
    else
	../bootstrap --parallel=4 --prefix=`pwd`/..
    fi
    make -j4
    make -j4 install
    export CMAKE=`pwd`/../bin/cmake
    export PATH=`pwd`/../bin:$PATH
    ${CMAKE} -version
    cd ../..
    return 0
}
