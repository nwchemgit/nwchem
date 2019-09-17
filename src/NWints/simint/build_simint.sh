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
SRC_HOME=`pwd`
DERIV=1
if [ $# -eq 0 ]; then
    MAXAM=3
else
    MAXAM=$1
fi
PERMUTE_SLOW=${MAXAM}
rm -rf simint.l${MAXAM}_p${PERMUTE_SLOW}_d${DERIVE}* simint-chem-simint-generator*
curl -LJ https://github.com/simint-chem/simint-generator/tarball/master -o simint-chem-simint-generator.tar.gz
tar xzf simint-chem-simint-generator.tar.gz
cd simint-chem-simint-generator-*
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
GOTSSE2=$(cat /proc/cpuinfo | egrep sse2 | tail -n 1 | awk ' /sse2/  {print "Y"}')
GOTAVX=$(cat /proc/cpuinfo | egrep avx | tail -n 1 | awk ' /avx/  {print "Y"}')
GOTAVX2=$(cat /proc/cpuinfo | egrep avx2 | tail -n 1 | awk ' /avx2/  {print "Y"}')
if [ ${GOTAVX2} == "Y" ]; then
    VEC=avx2
elif [ ${GOTAVX} == "Y" ]; then
    VEC=avx
elif [ ${GOTSSE2} == "Y" ]; then
    VEC=sse
else
    VEC=scalar
fi
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
