#!/bin/bash
# script to download simint-generator, create the simint library, compile it
# and link it in NWChem
# the only argument is the maximum angular moment (if not provided, default=3)
# FC=compilername can be used to set compiler, e.g.
#
# FC=ifort ./build_simint.sh 3
#
# ./build_simint.sh
#
UNAME_S=$(uname -s)
if [[ ${UNAME_S} == Linux ]]; then
    CPU_FLAGS=$(cat /proc/cpuinfo | egrep flags)
elif [[ ${UNAME_S} == Darwin ]]; then
    CPU_FLAGS=$(sysctl -a | grep machdep.cpu.features)
else
    echo Operating system not supported yet
    exit 1
fi
GOTSSE2=$(echo ${CPU_FLAGS} | tr  'A-Z' 'a-z'| grep sse2 | tail -n 1 | awk ' /sse2/  {print "Y"}')
 GOTAVX=$(echo ${CPU_FLAGS} | tr  'A-Z' 'a-z'| grep avx  | tail -n 1 | awk ' /avx/  {print "Y"}')
GOTAVX2=$(echo ${CPU_FLAGS} | tr  'A-Z' 'a-z'| grep avx2 | tail -n 1 | awk ' /avx2/  {print "Y";exit};{print "N"}')
GOTAVX512=$(echo ${CPU_FLAGS} | tr  'A-Z' 'a-z'| grep avx512 | tail -n 1 | awk ' /avx512f/  {print "Y";exit};{print "N"}')
if [[ ! -z "${GOTAVX512}" ]]; then
    VEC=avx512
elif [[ ! -z "${GOTAVX2}" ]]; then
    VEC=avx2
elif [[ ! -z "${GOTAVX}" ]]; then
    VEC=avx
elif [[ ! -z "${GOTSSE2}" ]]; then
    VEC=sse
else
    VEC=scalar
fi
SRC_HOME=`pwd`
DERIV=1
if [ $# -eq 0 ]; then
    MAXAM=3
else
    MAXAM=$1
fi
PERMUTE_SLOW=${MAXAM}
GITHUB_USERID=edoapra
rm -rf simint.l${MAXAM}_p${PERMUTE_SLOW}_d${DERIVE}* *-chem-simint-generator-?????? simint-chem-simint-generator.tar.gz
curl -LJ https://github.com/${GITHUB_USERID}/simint-generator/tarball/master -o simint-chem-simint-generator.tar.gz
#curl -LJ https://github.com/simint-chem/simint-generator/tarball/master -o simint-chem-simint-generator.tar.gz
tar xzf simint-chem-simint-generator.tar.gz
cd *-simint-generator-???????
rm -f generator_types.patch
echo > generator_types.patch <<EOF
--- simint-chem-simint-generator-c589bd7/generator/CommandLine.hpp	2018-12-11 10:48:31.000000000 -0800
+++ modif/generator/CommandLine.hpp	2019-09-17 09:25:45.000000000 -0700
@@ -10,6 +10,7 @@
 
 #include <vector>
 #include "generator/Options.hpp"
+#include "generator/Types.hpp"
 
 
 /*! \brief Get the next argument on the command line
EOF
patch -p1 < /tmp/generator_types.patch
pwd
mkdir build; cd build
 cmake ../
make -j3
cd ..
#./create.py -g build/generator/ostei -l 6 -p 4 -d 1 simint.l6_p4_d1
#create.py -g build/generator/ostei -l 4 -p 4 -d 0 -ve 4 -he 4 -vg 5 -hg 5
#https://www.cc.gatech.edu/~echow/pubs/huang-chow-sc18.pdf
/usr/bin/time -p ./create.py -g build/generator/ostei -l ${MAXAM} -p ${PERMUTE_SLOW} -d ${DERIV} ../simint.l${MAXAM}_p${PERMUTE_SLOW}_d${DERIV}  -ve 4 -he 4 -vg 5 -hg 5
cd ../simint.l${MAXAM}_p${PERMUTE_SLOW}_d${DERIV}
mkdir build
cd build
if [[ -z "${FC}" ]]; then
    #look for gfortran
    GOTGFORTRAN=$(command -v gfortran | tail -n 1 | awk ' /gfortran/  {print "Y"}')
    if  [ ${GOTGFORTRAN} == "Y" ]; then
	FC=gfortran
    else
	exit 1
    fi
fi    
if [ ${FC} == gfortran ] || [ ${FC} == flang ] ; then
    Fortran_FLAGS="-fdefault-integer-8 -cpp"
elif  [ ${FC} == ifort ]; then
    Fortran_FLAGS="-i8 -fpp"
fi
cmake \
 -DCMAKE_BUILD_TYPE=Release -DSIMINT_VECTOR=${VEC}  \
 -DCMAKE_INSTALL_LIBDIR=lib -DENABLE_FORTRAN=ON -DSIMINT_MAXAM=${MAXAM} SIMINT_MAXDER=${DERIV} \
 -DCMAKE_Fortran_FLAGS="$Fortran_FLAGS" -DCMAKE_INSTALL_PREFIX=${SRC_HOME}/simint.l${MAXAM}_p${PERMUTE_SLOW}_d${DERIV}.install ../
/usr/bin/time -p make -j3
make install
cd ../..
export USE_SIMINT=1
export SIMINT_HOME=${SRC_HOME}/simint.l${MAXAM}_p${PERMUTE_SLOW}_d${DERIV}.install
make clean
make V=1 
pwd
cd ../api
touch `egrep -l SIM *F`
make V=1 
cd ../..
make V=1  link
echo 'NWChem built with SIMINT support. Maximum angular momentum='${MAXAM}
