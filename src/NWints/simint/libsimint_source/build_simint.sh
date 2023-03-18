#!/usr/bin/env bash
# script to download simint-generator, create the simint library, compile it
# and link it in NWChem
# FC=compilername can be used to set compiler, e.g.
# SIMINT_MAXAM cane be used to set the maximum ang. momentum
# FC=ifort ./build_simint.sh 
#
#  SIMINT_MAXAM=5 ./build_simint.sh
#
mysimpwd=`pwd`
source ../../../libext/libext_utils/cmake.sh
cd $mysimpwd
if  [ -z "$(command -v python3)" ]; then
    echo python3 not installed
    echo please install python3
    exit 1
fi
if  [ -z "$(command -v curl)" ] && [ -z "$(command -v wget)" ]; then
    echo curl and wget not installed
    echo please install curl or wget
    exit 1
fi
if  [ -z "$(command -v patch)" ]; then
    echo patch not installed
    echo please install patch
    exit 1
fi
UNAME_S=$(uname -s)
if [[ ${UNAME_S} == Linux ]]; then
    CPU_FLAGS=$(cat /proc/cpuinfo | grep flags |tail -n 1)
    CPU_FLAGS_2=$(cat /proc/cpuinfo | grep flags |tail -n 1)
elif [[ ${UNAME_S} == Darwin ]]; then
    CPU_FLAGS=$(sysctl -n machdep.cpu.features)
    if [[ "$arch" == "x86_64" ]]; then
	CPU_FLAGS_2=$(sysctl -n machdep.cpu.leaf7_features)
    fi
else
    echo Operating system not supported yet
    exit 1
fi
  GOTSSE2=$(echo ${CPU_FLAGS}   | tr  'A-Z' 'a-z'| awk ' /sse2/   {print "Y"}')
   GOTAVX=$(echo ${CPU_FLAGS}   | tr  'A-Z' 'a-z'| awk ' /avx/    {print "Y"}')
  GOTAVX2=$(echo ${CPU_FLAGS_2} | tr  'A-Z' 'a-z'| awk ' /avx2/   {print "Y"}')
GOTAVX512=$(echo ${CPU_FLAGS}   | tr  'A-Z' 'a-z'| awk ' /avx512f/{print "Y"}')
   GOTSVE=$(echo ${CPU_FLAGS}   | tr  'A-Z' 'a-z'| awk ' /sve/{print "Y"}')
if [[ -n "${SIMINT_VECTOR}" ]]; then
      VEC=${SIMINT_VECTOR}
elif [[ "${GOTAVX512}" == "Y" ]]; then
    VEC=commonavx512
elif [[ "${GOTAVX2}" == "Y" ]]; then
    VEC=avx2
elif [[ "${GOTAVX}" == "Y" ]]; then
    VEC=avx
elif [[ "${GOTSSE2}" == "Y" ]]; then
    VEC=sse
elif [[ "${GOTSVE}" == "Y" ]]; then
    VEC=sve
else
    VEC=scalar
fi
echo VEC $VEC
if [[ "${VEC}" == "avx512" ]]; then
if [[   -z "${CC}" ]]; then
    CC=cc
fi
echo CC is $CC
GCC_EXTRA=$(echo $CC | cut -c 1-3)
if [ "$GCC_EXTRA" == gcc ]; then
let GCCVERSIONGT5=$(expr `${CC} -dumpversion | cut -f1 -d.` \> 5)
#echo exit code "$?"
    if [[ ${GCCVERSIONGT5} != 1 ]]; then
	echo
	echo you have gcc version $(${CC} -dumpversion | cut -f1 -d.)
	echo gcc version 6 and later needed for skylake
	echo
	exit 1
    fi
fi
fi
SRC_HOME=`pwd`
DERIV=1
if [[  -z "${SIMINT_MAXAM}" ]]; then
    SIMINT_MAXAM=3
fi
#PERMUTE_SLOW=1
PERMUTE_SLOW=${SIMINT_MAXAM}
GITHUB_USERID=edoapra
#GITHUB_USERID=simint-chem
#rm -rf simint.l${SIMINT_MAXAM}_p${PERMUTE_SLOW}_d${DERIVE}* *-chem-simint-generator-?????? simint-chem-simint-generator.tar.gz simint_lib
rm -rf simint.l${SIMINT_MAXAM}_p${PERMUTE_SLOW}_d${DERIVE}* *-chem-simint-generator-?????? simint_lib

GITHUB_URL=https://github.com/${GITHUB_USERID}/simint-generator/tarball/master
TAR_NAME=simint-chem-simint-generator.tar.gz
if [ -f  ${TAR_NAME} ]; then
    echo "using existing"  ${TAR_NAME}
else
    if  [ ! -z "$(command -v curl)" ] ; then
	curl -L "${GITHUB_URL}" -o "${TAR_NAME}"
    else
	wget -O "${TAR_NAME}" "${GITHUB_URL}"
    fi
fi
if [[ -z "${CMAKE}" ]]; then
    #look for cmake
    if [[ -z "$(command -v cmake)" ]]; then
	cmake_instdir=../../../libext/libext_utils
	get_cmake_release $cmake_instdir
	status=$?
	if [ $status -ne 0 ]; then
	    echo cmake required to build simint
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
if ((CMAKE_VER_MAJ < 3)) || (((CMAKE_VER_MAJ > 2) && (CMAKE_VER_MIN < 21))); then
    cmake_instdir=../../../libext/libext_utils/
    get_cmake_release $cmake_instdir
    status=$?
    if [ $status -ne 0 ]; then
	echo cmake 3.21 required to build simint
	echo Please install cmake
	echo define the CMAKE env. variable
	exit 1
    fi
fi
cd $mysimpwd
tar xzf simint-chem-simint-generator.tar.gz
cd *-simint-generator-???????
pwd
if [[  -z "${NWCHEM_TOP}" ]]; then
    dir4=$(dirname `pwd`)
    dir3=$(dirname "$dir4")
    dir2=$(dirname "$dir3")
    dir1=$(dirname "$dir2")
    NWCHEM_TOP=$(dirname "$dir1")
fi
mkdir -p build; cd build
if [[ -z "${SIMINT_BUILD_TYPE}" ]]; then
    SIMINT_BUILD_TYPE=Release
fi
if [[ -z "${CXX}" ]]; then
    #look for c++
    if  [ -z "$(command -v c++)" ]; then
        echo c++ not installed
        echo please install a C++ compiler and
        echo define the CXX env. variable
	exit 1
    else
	CXX=c++
    fi
fi    
if [[ -z "${CXX_FOR_BUILD}" ]]; then
    CXX_FOR_BUILD=${CXX}
fi
echo CXX_FOR_BUILD $CXX_FOR_BUILD && $CMAKE CXX=$CXX_FOR_BUILD -DCMAKE_CXX_COMPILER=$CXX_FOR_BUILD  -DCMAKE_BUILD_TYPE="${SIMINT_BUILD_TYPE}"  ../
make -j2
cd ..
#./create.py -g build/generator/ostei -l 6 -p 4 -d 1 simint.l6_p4_d1
#create.py -g build/generator/ostei -l 4 -p 4 -d 0 -ve 4 -he 4 -vg 5 -hg 5
#https://www.cc.gatech.edu/~echow/pubs/huang-chow-sc18.pdf
#workaround for PYTHONHOME crazyness
if [[ ! -z "${PYTHONHOME}" ]]; then
    export PYTHONHOMESET=${PYTHONHOME}
    unset PYTHONHOME
    echo 'PYTHONOME unset'
fi
if [[ -z "${GENERATOR_PROCESSES}" ]]; then
    GENERATOR_PROCESSES=3
    #parallel processing broken for g++-10 and later (at least on macos)
    if [[ $(expr `${CXX} -dumpversion | cut -f1 -d.` \> 9) == 1 ]]; then
	GENERATOR_PROCESSES=1
    fi
fi
echo GENERATOR_PROCESSES is ${GENERATOR_PROCESSES}
time -p ./create.py -g build/generator/ostei -l ${SIMINT_MAXAM} -p ${PERMUTE_SLOW} -d ${DERIV} ../simint.l${SIMINT_MAXAM}_p${PERMUTE_SLOW}_d${DERIV}  -ve 4 -he 4 -vg 5 -hg 5 -n ${GENERATOR_PROCESSES}
if [[ ! -z "${PYTHONHOME}" ]]; then
    export PYTHONHOME=${PYTHONHOMESET}
    unset PYTHONHOMESET
    echo 'PYTHONOME set'
fi
cd ../simint.l${SIMINT_MAXAM}_p${PERMUTE_SLOW}_d${DERIV}
mkdir -p build
cd build
if [[ -z "${FC}" ]]; then
    #look for gfortran
    if  [ -z "$(command -v gfortran)" ]; then
        echo gfortran not installed
        echo please install a Fortran compiler and
        echo define the FC env. variable
	exit 1
    else
	FC=gfortran
    fi
fi    
FC_EXTRA=$(${NWCHEM_TOP}/src/config/strip_compiler.sh ${FC})
echo FC_EXTRA $FC_EXTRA
if [[ ${FC_EXTRA} == gfortran  || ${FC_EXTRA} == flang || ${FC_EXTRA} == armflang || (${FC} == ftn && ${PE_ENV} == GNU) || (${FC} == ftn && ${PE_ENV} == AOCC) ]] ; then
    Fortran_FLAGS="-fdefault-integer-8 -cpp"
    if [[ ${FC_EXTRA} == gfortran  || (${FC} == ftn && ${PE_ENV} == GNU)]]; then
    GNUMAJOR=$(${FC} -dM -E - < /dev/null 2> /dev/null | grep __GNUC__ |cut -c18-)
    echo GNUMAJOR is $GNUMAJOR
    if [ $GNUMAJOR -ge 8 ]; then
    Fortran_FLAGS+=" -std=legacy "
    fi
fi
elif  [ ${FC} == xlf ] || [ ${FC} == xlf_r ] || [ ${FC} == xlf90 ]|| [ ${FC} == xlf90_r ]; then
    Fortran_FLAGS=" -qintsize=8 -qextname -qpreprocess"
elif  [[ ${FC} == ifort || (${FC} == ftn && ${PE_ENV} == INTEL) ]]; then
    Fortran_FLAGS="-i8 -fpp"
elif  [ ${FC} == ftn ]  && [ ${PE_ENV} == CRAY  ]; then
    Fortran_FLAGS=" -ffree -s integer64 -e F "
elif  [[ ${FC_EXTRA} == nvfortran || ${FC} == pgf90 || (${FC} == ftn && ${PE_ENV} == NVIDIA) ]]; then
    Fortran_FLAGS="-i8 -cpp"
    CC=gcc
    CXX=g++
    if  [[ ${PE_ENV} == NVIDIA ]]; then
	unset CPATH
    fi
elif  [ ${FC} == frt ] || [ ${FC} == frtpx ] ; then
    Fortran_FLAGS=" -fs -CcdLL8 -CcdII8 -cpp "
    CC=/opt/FJSVxos/devkit/aarch64/bin/aarch64-linux-gnu-gcc
    CXX=/opt/FJSVxos/devkit/aarch64/bin/aarch64-linux-gnu-g++
fi
if [[ ! -z ${FFLAGS_FORGA} ]]; then Fortran_FLAGS+=" ${FFLAGS_FORGA}" ; fi
echo Fortran_FLAGS equal "$Fortran_FLAGS"
FC="${FC}" CC="${CC}" CXX="${CXX}" $CMAKE \
 -DCMAKE_BUILD_TYPE="${SIMINT_BUILD_TYPE}" -DSIMINT_VECTOR=${VEC}  \
 -DCMAKE_INSTALL_LIBDIR=lib -DENABLE_FORTRAN=ON -DSIMINT_MAXAM=${SIMINT_MAXAM} -DSIMINT_MAXDER=${DERIV} \
 -DENABLE_TESTS=OFF     -DSIMINT_STANDALONE=OFF   \
 -DCMAKE_Fortran_FLAGS="$Fortran_FLAGS" -DCMAKE_INSTALL_PREFIX=${SRC_HOME}/simint.l${SIMINT_MAXAM}_p${PERMUTE_SLOW}_d${DERIV}.install ../
time -p make  -j4
make simint
if [[ (${FC} == ftn && ${PE_ENV} == CRAY) ]] ; then
    cp simint/SIMINTFORTRAN.mod simint/simintfortran.mod
fi
make install/fast
cd ../..
echo ln -sf  simint.l${SIMINT_MAXAM}_p${PERMUTE_SLOW}_d${DERIV}.install simint_install
ln -sf  simint.l${SIMINT_MAXAM}_p${PERMUTE_SLOW}_d${DERIV}.install simint_install
cd simint_install/lib
strip --strip-debug libsimint.a
ln -sf libsimint.a libnwc_simint.a
export SIMINT_HOME=${SRC_HOME}/simint.l${SIMINT_MAXAM}_p${PERMUTE_SLOW}_d${DERIV}.install
echo 'SIMINT library built with maximum angular momentum='${SIMINT_MAXAM}
echo SIMINT_HOME="$SIMINT_HOME"
exit 0
# remainder of script not used since flow goes back to makefile
export USE_SIMINT=1
make clean
make V=1 
pwd
cd ../api
touch `egrep -l SIM *F`
make V=1 
cd ../..
make V=1  link
echo 'NWChem built with SIMINT support. Maximum angular momentum='${SIMINT_MAXAM}
echo SIMINT_HOME="$SIMINT_HOME"
