#- env: == default ==
os=`uname`
export NWCHEM_TOP=$TRAVIS_BUILD_DIR
export USE_MPI=y
export USE_64TO32=y
 if [[ "$os" == "Darwin" ]]; then 
  export BLASOPT="-L/usr/local/opt/openblas/lib -lopenblas"
  if [[ "$MPI_IMPL" == "openmpi" ]]; then
    export SCALAPACK="-L/usr/local/lib -lscalapack -lopenblas"
  fi
  export NWCHEM_TARGET=MACX64 
  export DYLD_LIBRARY_PATH=$TRAVIS_BUILD_DIR/lib:$DYLD_LIBRARY_PATH
  if [[ "$MPI_IMPL" == "mpich" ]]; then 
    export PATH=/usr/local/opt/mpich/bin/:$PATH 
  fi
  if [[ "$MPI_IMPL" == "openmpi" ]]; then 
    export PATH=/usr/local/opt/open-mpi/bin/:$PATH 
  fi
fi
if [[ "$os" == "Linux" ]]; then 
   export BLASOPT="-L$TRAVIS_BUILD_DIR/lib -lopenblas"
   export SCALAPACK="-L$TRAVIS_BUILD_DIR/lib  -lscalapack -lopenblas"
   export NWCHEM_TARGET=LINUX64 
   export LD_LIBRARY_PATH=$TRAVIS_BUILD_DIR/lib:$LD_LIBRARY_PATH
fi
export OMP_NUM_THREADS=1
export USE_NOIO=1
export BLAS_SIZE=4
export SCALAPACK_SIZE=4
export USE_PYTHONCONFIG=y
export PYTHONVERSION=2.7
export PYTHONHOME=/usr
