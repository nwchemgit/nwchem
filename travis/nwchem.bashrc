#- env: == default ==
export NWCHEM_TOP=$TRAVIS_BUILD_DIR
export NWCHEM_MODULES="qmandpw"
export USE_MPI=y
export USE_64TO32=y
export BLASOPT="-L$TRAVIS_BUILD_DIR/lib -lopenblas -lpthread -lrt"
export BLAS_SIZE=4
export SCALAPACK="-L$TRAVIS_BUILD_DIR/lib  -lscalapack -lnwclapack -lopenblas -lpthread -lrt"
export SCALAPACK_SIZE=4
export USE_PYTHONCONFIG=y
export PYTHONVERSION=2.7
export PYTHONHOME=/usr
