#- env: == default ==
os=`uname`
arch=`uname -m`
#export NWCHEM_TOP=$TRAVIS_BUILD_DIR
#TARBALL=https://github.com/nwchemgit/nwchem/releases/download/v7.0.0-beta1/nwchem-7.0.0-release.revision-5bcf0416-src.2019-11-01.tar.bz2
export USE_MPI=y
if [[ "$os" == "Darwin" ]]; then 
   if [[ "$NWCHEM_MODULES" = "tce" && -z "$TARBALL" ]]; then
     export USE_INTERNALBLAS=y
   else
     export USE_64TO32="y"
     export BLASOPT="-L/usr/local/opt/openblas/lib -lopenblas"
     export LAPACK_LIB="-L/usr/local/opt/openblas/lib -lopenblas"
     if [[ "$MPI_IMPL" == "openmpi" ]]; then
       export SCALAPACK="-L/usr/local/lib -lscalapack -lopenblas"
     fi
   fi
  export NWCHEM_TARGET=MACX64 
  export DYLD_LIBRARY_PATH=$TRAVIS_BUILD_DIR/lib:$DYLD_LIBRARY_PATH
  if [[ "$MPI_IMPL" == "openmpi" ]]; then
    export PATH=/usr/local/opt/open-mpi/bin/:$PATH 
  fi
  if [[ "$MPI_IMPL" == "mpich" ]]; then 
    export PATH=/usr/local/opt/mpich/bin/:$PATH 
  fi
  export PATH=/usr/local/opt/python@3.8/bin:$PATH
  export PYTHONVERSION=3.8

fi
if [[ "$os" == "Linux" ]]; then 
   export NWCHEM_TARGET=LINUX64 
   if [[ -z "$USE_SIMINT" ]] && [[ "$arch" != "aarch64" ]] && [[ "$NWCHEM_MODULES" != "tce" ]] ; then 
     export BUILD_OPENBLAS="y"
     export BUILD_SCALAPACK="y"
   else
     export BLASOPT="-lopenblas"
     export LAPACK_LIB="-lopenblas"
     if [[ "$MPI_IMPL" == "mpich" ]]; then 
       export SCALAPACK="-lscalapack-mpich -lopenblas"
     elif [[ "$MPI_IMPL" == "openmpi" ]]; then
       export SCALAPACK="-lscalapack-openmpi -lopenblas"
     fi
   fi
     export USE_64TO32="y"
     if [[ "$arch" == "aarch64" ]]; then
         if [[ "$NWCHEM_MODULES" == "tce" && -z "$TARBALL" ]]; then
	     unset BLASOPT
	     unset LAPACK_LIB
	     unset SCALAPACK
	     export USE_INTERNALBLAS=y
	     unset USE_64TO32
	     unset BLAS_SIZE
	     unset SCALAPACK_SIZE
	 fi
     fi
     
#   fi
fi
export OMP_NUM_THREADS=1
export USE_NOIO=1
if [[ "$USE_64TO32" == "y" ]]; then
  export BLAS_SIZE=4
  export SCALAPACK_SIZE=4
fi
export NWCHEM_EXECUTABLE=$TRAVIS_BUILD_DIR/.cachedir/binaries/$NWCHEM_TARGET/nwchem_"$arch"_`echo $NWCHEM_MODULES|sed 's/ /-/g'`_"$MPI_IMPL"
