#- env: == default ==
os=`uname`
arch=`uname -m`
if [[ -z "$TRAVIS_BUILD_DIR" ]] ; then
    TRAVIS_BUILD_DIR=$(pwd)
fi
export NWCHEM_TOP=$TRAVIS_BUILD_DIR
#TARBALL=https://github.com/nwchemgit/nwchem/releases/download/v7.0.0-beta1/nwchem-7.0.0-release.revision-5bcf0416-src.2019-11-01.tar.bz2
export USE_MPI=y
if [[ "$os" == "Darwin" ]]; then 
  export NWCHEM_TARGET=MACX64 
  export DYLD_LIBRARY_PATH=$TRAVIS_BUILD_DIR/lib:$DYLD_LIBRARY_PATH
  if [[ "$MPI_IMPL" == "openmpi" ]]; then
    export PATH=/usr/local/opt/open-mpi/bin/:$PATH 
  fi
  if [[ "$MPI_IMPL" == "mpich" ]]; then 
    export PATH=/usr/local/opt/mpich/bin/:$PATH 
  fi
  #python3 on brew
  export PATH=/usr/local/bin:$PATH

fi
if [[ "$os" == "Linux" ]]; then 
   export NWCHEM_TARGET=LINUX64 
fi
export OMP_NUM_THREADS=1
export USE_NOIO=1
if [[ "$BLAS_SIZE" == "4" ]]; then
  export USE_64TO32=y
fi
if [[ -z "$USE_INTERNALBLAS" ]]; then
  export BUILD_OPENBLAS="y"
  export BUILD_SCALAPACK="y"
  export BLAS_SIZE=8
  export SCALAPACK_SIZE=8
fi
export NWCHEM_EXECUTABLE=$TRAVIS_BUILD_DIR/.cachedir/binaries/$NWCHEM_TARGET/nwchem_"$arch"_`echo $NWCHEM_MODULES|sed 's/ /-/g'`_"$MPI_IMPL"
