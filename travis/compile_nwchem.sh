#!/usr/bin/env bash
echo "start compile"
echo "BLAS_SIZE is " $BLAS_SIZE
set -e
# source env. variables
if [[ -z "$TRAVIS_BUILD_DIR" ]] ; then
    TRAVIS_BUILD_DIR=$(pwd)
fi
echo TRAVIS_BUILD_DIR is $TRAVIS_BUILD_DIR
source $TRAVIS_BUILD_DIR/travis/nwchem.bashrc
echo ============================================================
env|egrep BLAS     || true
env|egrep USE_6   || true
ls -lrt $TRAVIS_BUILD_DIR|tail -3
echo ============================================================
os=`uname`
arch=`uname -m`
if [[ "$NWCHEM_MODULES" == "tce" ]]; then 
    export EACCSD=1
    export IPCCSD=1
fi
cd $TRAVIS_BUILD_DIR/src
#FDOPT="-O0 -g"
export MPICH_FC=$FC
if [[ "$arch" == "aarch64" ]]; then 
    if [[ "$FC" == "flang" ]]  ; then
	export BUILD_MPICH=1
        FOPT="-O2  -ffast-math"
    elif [[ "$(basename -- $FC | cut -d \- -f 1)" == "nvfortran" ]] ; then
	export USE_FPICF=1
#	export MPICH_FC=nvfortran
	env|egrep FC
    else
#should be gfortran	
	if [[ "$NWCHEM_MODULES" == "tce" ]]; then 
	    FOPT="-O0 -fno-aggressive-loop-optimizations"
	else
	    FOPT="-O1 -fno-aggressive-loop-optimizations"
	fi
    fi
else
    if [[ "$FC" == "ifort" ]] || [[ "$FC" == "ifx" ]] ; then
	FOPT=-O2
    if [[ -z "$BUILD_OPENBLAS" ]] ; then
	if [[ "$os" == "Darwin" ]]; then
	    export BUILD_MPICH=1
 	    export BLASOPT="-L$MKLROOT  -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core  -lpthread -lm -ldl"
	else
	    export USE_FPICF=Y
            export BLASOPT="-L$MKLROOT/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core  -lpthread -lm -ldl"
            export SCALAPACK_LIB="-L$MKLROOT/lib/intel64 -lmkl_scalapack_ilp64 -lmkl_blacs_intelmpi_ilp64 -lpthread -lm -ldl"
	    export SCALAPACK_SIZE=8
	    unset BUILD_SCALAPACK
	fi
        unset BUILD_OPENBLAS
	export BLAS_SIZE=8
	export LAPACK_LIB="$BLASOPT"
    fi
	export I_MPI_F90="$FC"
    elif [[ "$FC" == "flang" ]] || [[ "$(basename -- $FC | cut -d \- -f 1)" == "nvfortran" ]] ; then
	export BUILD_MPICH=1
        if [[ "$FC" == "flang" ]]; then
	    FOPT="-O2  -ffast-math"
	fi
        if [[ "$FC" == "nvfortran" ]]; then
	    export USE_FPICF=1
#	    FOPT="-O2 -tp haswell"
	fi
    fi
fi
#check linear algebra
if [[ -z "$BLASOPT" ]] && [[ -z "$BUILD_OPENBLAS" ]] && [[ -z "$USE_INTERNALBLAS" ]] ; then
    echo " no existing BLAS settings, defaulting to BUILD_OPENBLAS=y "
    export BUILD_OPENBLAS=1
fi
# install armci-mpi if needed
if [[ "$ARMCI_NETWORK" == "ARMCI" ]]; then
    cd tools
    ./install-armci-mpi
    export EXTERNAL_ARMCI_PATH=$NWCHEM_TOP/external-armci
    cd ..
fi    

#compilation
 if [[ "$os" == "Darwin" ]]; then 
   if [[ "$NWCHEM_MODULES" == "tce" ]]; then
     FOPT="-O1 -fno-aggressive-loop-optimizations"
   fi
   if [[ ! -z "$USE_SIMINT" ]] ; then 
       FOPT="-O0 -fno-aggressive-loop-optimizations"
       SIMINT_BUILD_TYPE=Debug
       export PATH="/usr/local/bin:$PATH"
#       export LDFLAGS="-L/usr/local/opt/python@3.7/lib:$LDFLAGS"
   fi
   if [[ -z "$TRAVIS_HOME" ]]; then
       env
       mkdir -p ../bin/MACX64
       gcc -o ../bin/MACX64/depend.x config/depend.c
       make nwchem_config
       cd libext   && make V=-1  && cd ..
       cd tools    && make V=-1  && cd ..
       nohup make USE_INTERNALBLAS=y deps_stamp  >& deps.log &
       sleep 120s
       echo tail deps.log '@@@'
       tail -10  deps.log
       echo done tail deps.log '@@@'
       export QUICK_BUILD=1
       if [[ -z "$FOPT" ]]; then
	   make V=0   -j3
       else
	   make V=0 FOPTIMIZE="$FOPT"   -j3
       fi
   else
       ../travis/sleep_loop.sh make V=1 FOPTIMIZE="$FOPT"   -j3
   fi
     unset QUICK_BUILD
     cd $TRAVIS_BUILD_DIR/src/64to32blas 
     make
     cd $TRAVIS_BUILD_DIR/src
     ../contrib/getmem.nwchem 1000
     otool -L ../bin/MACX64/nwchem
#     printenv DYLD_LIBRARY_PATH
#     ls -lrt $DYLD_LIBRARY_PATH
#      tail -120 make.log
 elif [[ "$os" == "Linux" ]]; then
     export MAKEFLAGS=-j3
     echo    "$FOPT$FDOPT"
if [[ -z "$TRAVIS_HOME" ]]; then
    mkdir -p ../bin/LINUX64
    gcc -o ../bin/LINUX64/depend.x config/depend.c
    make nwchem_config
    cd libext   && make V=-1  && cd ..
    cd tools    && make V=-1  && cd ..
    nohup make USE_INTERNALBLAS=y deps_stamp  >& deps.log &
    cd hessian
    nohup make USE_INTERNALBLAS=y dependencies include_stamp >& ../deps2.log &
    cd ../nwdft/xc
    nohup make USE_INTERNALBLAS=y dependencies include_stamp >& ../../deps3.log &
    cd ../..
    sleep 360s
    echo tail deps.log '11@@@'
    tail -10  deps.log
    echo done tail deps.log '11@@@'
    echo tail deps2.log '11@@@'
    tail deps2.log || true
    echo tail deps3.log '11@@@'
    tail deps3.log || true
    export QUICK_BUILD=1
    if [[ -z "$FOPT" ]]; then
	make V=0   -j3
    else
	make V=0 FOPTIMIZE="$FOPT"   -j3
    fi
else
    ../travis/sleep_loop.sh make V=1 FOPTIMIZE="$FOPT"  -j3
fi
     unset QUICK_BUILD
     cd $TRAVIS_BUILD_DIR/src/64to32blas 
     make
     cd $TRAVIS_BUILD_DIR/src
     $TRAVIS_BUILD_DIR/contrib/getmem.nwchem 1000
 fi
 #caching
 mkdir -p $TRAVIS_BUILD_DIR/.cachedir/binaries/$NWCHEM_TARGET $TRAVIS_BUILD_DIR/.cachedir/files
 cp $TRAVIS_BUILD_DIR/bin/$NWCHEM_TARGET/nwchem  $NWCHEM_EXECUTABLE
 echo === ls binaries cache ===
 ls -lrt $TRAVIS_BUILD_DIR/.cachedir/binaries/$NWCHEM_TARGET/ 
 echo =========================
 rsync -av $TRAVIS_BUILD_DIR/src/basis/libraries.bse  $TRAVIS_BUILD_DIR/.cachedir/files/.
 rsync -av $TRAVIS_BUILD_DIR/src/basis/libraries  $TRAVIS_BUILD_DIR/.cachedir/files/.
 rsync -av $TRAVIS_BUILD_DIR/src/nwpw/libraryps  $TRAVIS_BUILD_DIR/.cachedir/files/.
