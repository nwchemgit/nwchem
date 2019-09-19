#!/bin/bash
# script to download simint-generator, create the simint library, compile it
# and link it in NWChem
# FC=compilername can be used to set compiler, e.g.
# SIMINT_MAXAM cane be used to set the maximum ang. momentum
# FC=ifort ./build_simint.sh 
#
#  SIMINT_MAXAM=5 ./build_simint.sh
#
if  [ -z "$(command -v python3)" ]; then
    echo python3 not installed
    echo please install python3
    exit 1
fi
UNAME_S=$(uname -s)
if [[ ${UNAME_S} == Linux ]]; then
    CPU_FLAGS=$(cat /proc/cpuinfo | egrep flags)
    CPU_FLAGS_2=$(cat /proc/cpuinfo | egrep flags)
elif [[ ${UNAME_S} == Darwin ]]; then
    CPU_FLAGS=$(sysctl -n machdep.cpu.features)
    CPU_FLAGS_2=$(sysctl -n machdep.cpu.leaf7_features)
else
    echo Operating system not supported yet
    exit 1
fi
GOTSSE2=$(echo ${CPU_FLAGS} | tr  'A-Z' 'a-z'|  awk ' /sse2/  {print "Y"}')
 GOTAVX=$(echo ${CPU_FLAGS} | tr  'A-Z' 'a-z'|  awk ' /avx/  {print "Y"}')
GOTAVX2=$(echo ${CPU_FLAGS_2} | tr  'A-Z' 'a-z'| awk ' /avx2/  {print "Y"}')
GOTAVX512=$(echo ${CPU_FLAGS} | tr  'A-Z' 'a-z'| awk ' /avx512f/  {print "Y"}')
GOTAVX512="N"
if [[ "${GOTAVX512}" == "Y" ]]; then
    VEC=avx512
elif [[ "${GOTAVX2}" == "Y" ]]; then
    VEC=avx2
elif [[ "${GOTAVX}" == "Y" ]]; then
    VEC=avx
elif [[ "${GOTSSE2}" == "Y" ]]; then
    VEC=sse
else
    VEC=scalar
fi
SRC_HOME=`pwd`
DERIV=1
if [[  -z "${SIMINT_MAXAM}" ]]; then
    SIMINT_MAXAM=3
fi
#PERMUTE_SLOW=1
PERMUTE_SLOW=${SIMINT_MAXAM}
GITHUB_USERID=edoapra
rm -rf simint.l${SIMINT_MAXAM}_p${PERMUTE_SLOW}_d${DERIVE}* *-chem-simint-generator-?????? simint-chem-simint-generator.tar.gz simint_lib
curl -L https://github.com/${GITHUB_USERID}/simint-generator/tarball/master -o simint-chem-simint-generator.tar.gz
#curl -LJ https://github.com/simint-chem/simint-generator/tarball/master -o simint-chem-simint-generator.tar.gz
tar xzf simint-chem-simint-generator.tar.gz
cd *-simint-generator-???????
rm -f generator_types.patch
cat > generator_types.patch <<EOF
--- simint-chem-simint-generator-c589bd7/generator/CommandLine.hpp	2018-12-11 10:48:31.000000000 -0800
+++ modif/generator/CommandLine.hpp	2019-09-17 09:25:45.000000000 -0700
@@ -10,6 +10,7 @@
 
 #include <vector>
 #include "generator/Options.hpp"
+#include "generator/Types.hpp"
 
 
 /*! \brief Get the next argument on the command line
EOF
patch -p1 < ./generator_types.patch
pwd
mkdir -p build; cd build
if [[ -z "${CMAKE}" ]]; then
    #look for cmake
    if [[ -z "$(command -v cmake)" ]]; then
	echo cmake required to build Simint
	echo Please install cmake
	echo define the CMAKE env. variable
	exit 1
    else
	CMAKE=cmake
    fi
fi
CMAKE_VER=$(${CMAKE} --version|cut -d " " -f 3|head -1|cut -c1)
#echo CMAKE_VER is ${CMAKE_VER}
if [[ ${CMAKE_VER} -lt 3 ]]; then
    echo CMake 3.0.2 or higher is required
    echo Please install CMake 3
    echo define the CMAKE env. variable
    exit 1
fi
$CMAKE ../
make -j2
cd ..
#./create.py -g build/generator/ostei -l 6 -p 4 -d 1 simint.l6_p4_d1
#create.py -g build/generator/ostei -l 4 -p 4 -d 0 -ve 4 -he 4 -vg 5 -hg 5
#https://www.cc.gatech.edu/~echow/pubs/huang-chow-sc18.pdf
time -p ./create.py -g build/generator/ostei -l ${SIMINT_MAXAM} -p ${PERMUTE_SLOW} -d ${DERIV} ../simint.l${SIMINT_MAXAM}_p${PERMUTE_SLOW}_d${DERIV}  -ve 4 -he 4 -vg 5 -hg 5
#time -p ./create.py -g build/generator/ostei -l ${SIMINT_MAXAM} -p ${PERMUTE_SLOW} -d ${DERIV} ../simint.l${SIMINT_MAXAM}_p${PERMUTE_SLOW}_d${DERIV}  
cd ../simint.l${SIMINT_MAXAM}_p${PERMUTE_SLOW}_d${DERIV}
mkdir -p build
cd build
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
if [ ${FC} == gfortran ] || [ ${FC} == flang ] ; then
    Fortran_FLAGS="-fdefault-integer-8 -cpp"
elif  [ ${FC} == ifort ]; then
    Fortran_FLAGS="-i8 -fpp"
fi
FC="${FC}" CXX="${CXX}" $CMAKE \
 -DCMAKE_BUILD_TYPE=Release -DSIMINT_VECTOR=${VEC}  \
 -DCMAKE_INSTALL_LIBDIR=lib -DENABLE_FORTRAN=ON -DSIMINT_MAXAM=${SIMINT_MAXAM} SIMINT_MAXDER=${DERIV} \
 -DCMAKE_Fortran_FLAGS="$Fortran_FLAGS" -DCMAKE_INSTALL_PREFIX=${SRC_HOME}/simint.l${SIMINT_MAXAM}_p${PERMUTE_SLOW}_d${DERIV}.install ../
time -p make  -j2
make install
cd ../..
echo ln -sf  simint.l${SIMINT_MAXAM}_p${PERMUTE_SLOW}_d${DERIV}.install simint_install
ln -sf  simint.l${SIMINT_MAXAM}_p${PERMUTE_SLOW}_d${DERIV}.install simint_install
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
