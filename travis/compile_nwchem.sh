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
        FOPT="-O2  -ffast-math"
    elif [[ "$(basename -- $FC | cut -d \- -f 1)" == "nvfortran" ]] ; then
	export USE_FPICF=1
	export MPICH_FC=nvfortran
	env|egrep FC
	nvfortran -V
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
#	FOPT=-O2
    if [[ -z "$BUILD_OPENBLAS" ]] ; then
	if [[ "$os" == "Darwin" ]]; then
 	    export BLASOPT="-L$MKLROOT  -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core  -lpthread -lm -ldl  -L`gfortran -print-file-name=libquadmath.0.dylib|sed -e s'/libquadmath.0.dylib//'` "
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
        if [[ "$FC" == "flang" ]]; then
	    FOPT="-O2  -ffast-math"
	fi
        if [[ "$FC" == "nvfortran" ]]; then
	    export USE_FPICF=1
            export MPICH_FC=nvfortran
	    nvfortran -V
	fi
    fi
fi
if [[ "$BLAS_ENV" == lib*openblas* ]] || [[ "$BLAS_ENV" == "brew_openblas" ]]; then
    if [[ "$BLAS_ENV" == *openblas64* ]]; then
        myob="openblas64"
    else
        myob="openblas"
    fi
    if [[ "$BLAS_ENV" == "brew_openblas" ]]; then
	if [ -z "$HOMEBREW_PREFIX" ] ; then
	    HOMEBREW_PREFIX=/usr/local
	fi
	export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$HOMEBREW_PREFIX/opt/openblas/lib/pkgconfig
    fi
    export BLASOPT=$(pkg-config --libs $myob)
    export LAPACK_LIB=$BLASOPT
    echo "BLASOPT and LAPACK_LIB are" $BLASOPT $LAPACK_LIB
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
# try to use ubuntu flaky GA pkg 
if [[ "$ARMCI_NETWORK" == "GA_DEBIAN" ]]; then
    export EXTERNAL_GA_PATH=/usr
    export EXTERNAL_ARMCI_PATH=/usr
    unset ARMCI_NETWORK
fi    

if [[ "$FC" == "gfortran" ]]; then
   if [[ "$($FC -dM -E - < /dev/null 2> /dev/null | grep __GNUC__ |cut -c 18-)" -lt 9 ]]; then
#disable xtb  if gfortran version < 9
     unset USE_TBLITE
     export NWCHEM_MODULES=$(echo $NWCHEM_MODULES |sed  's/xtb//')
   fi
fi


#compilation
export MAKEFLAGS=-j3
echo    "$FOPT$FDOPT"
    python3 -V || true
    python3-config  --ldflags || true
    if [[ -z "$FOPT" ]]; then
	make V=-1   -j3
    else
	make V=-1 FOPTIMIZE="$FOPT"   -j3
    fi
     cd $TRAVIS_BUILD_DIR/src/64to32blas 
     make
     cd $TRAVIS_BUILD_DIR/src
     $TRAVIS_BUILD_DIR/contrib/getmem.nwchem 1000
 #caching
 mkdir -p $TRAVIS_BUILD_DIR/.cachedir/binaries/$NWCHEM_TARGET $TRAVIS_BUILD_DIR/.cachedir/files
 cp $TRAVIS_BUILD_DIR/bin/$NWCHEM_TARGET/nwchem  $NWCHEM_EXECUTABLE
 echo === ls binaries cache ===
 ls -lrt $TRAVIS_BUILD_DIR/.cachedir/binaries/$NWCHEM_TARGET/ 
 echo =========================
 rsync -av $TRAVIS_BUILD_DIR/src/basis/libraries.bse  $TRAVIS_BUILD_DIR/.cachedir/files/.
 rsync -av $TRAVIS_BUILD_DIR/src/basis/libraries  $TRAVIS_BUILD_DIR/.cachedir/files/.
 rsync -av $TRAVIS_BUILD_DIR/src/nwpw/libraryps  $TRAVIS_BUILD_DIR/.cachedir/files/.
