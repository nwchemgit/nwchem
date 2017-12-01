#- env: == default ==
os=`uname`
export NWCHEM_TOP=$TRAVIS_BUILD_DIR
export USE_MPI=y
export USE_64TO32=y
 if [[ "$os" == "Darwin" ]]; then 
export BLASOPT="-lopenblas"
export SCALAPACK=" -lscalapack -lopenblas"
fi
if [[ "$os" == "Linux" ]]; then 
   export BLASOPT="-L$TRAVIS_BUILD_DIR/lib -lopenblas"
   export SCALAPACK="-L$TRAVIS_BUILD_DIR/lib  -lscalapack -lopenblas"
fi
export BLAS_SIZE=4
export SCALAPACK_SIZE=4
export USE_PYTHONCONFIG=y
export PYTHONVERSION=2.7
export PYTHONHOME=/usr
