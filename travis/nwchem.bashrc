#- env: == default ==
os=`uname`
arch=`uname -m`
if test -f "/usr/lib/os-release"; then
    dist=$(grep ID= /etc/os-release |head -1 |cut -c4-| sed 's/\"//g')
fi
if [ -z "$DISTR" ] ; then
    DISTR=$dist
fi
echo DISTR is "$DISTR"
if [[ -z "$TRAVIS_BUILD_DIR" ]] ; then
    TRAVIS_BUILD_DIR=$(pwd)
else
    export NWCHEM_TOP=$TRAVIS_BUILD_DIR
fi
echo NWCHEM_TOP is $NWCHEM_TOP
#TARBALL=https://github.com/nwchemgit/nwchem/releases/download/v7.0.0-beta1/nwchem-7.0.0-release.revision-5bcf0416-src.2019-11-01.tar.bz2
export USE_MPI=y
if [[ "$FC" == "flang" ]]; then
    export PATH=/usr/lib/aomp_11.12-0/bin/:$PATH
#    export PATH=/opt/rocm-4.0.0/llvm/bin:$PATH
fi
if [[ "$FC" == "ifort" ]]; then
    source /opt/intel/oneapi/compiler/latest/env/vars.sh
    ifort -V
    if [ -f /opt/intel/oneapi/mkl/latest/env/vars.sh ] ; then
	source /opt/intel/oneapi/mkl/latest/env/vars.sh
    fi
fi
if [[ "$MPI_IMPL" == "intel" ]]; then
    source /opt/intel/oneapi/mpi/latest/env/vars.sh
    mpif90 -v
    mpif90 -show
    if [ -f /opt/intel/oneapi/mkl/latest/env/vars.sh ] ; then
	source /opt/intel/oneapi/mkl/latest/env/vars.sh
    fi
fi
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

if [[ "$DISTR" == "fedora" ]] || [[ "$DISTR" == "centos" ]]; then
    export PATH=/usr/lib64/"$MPI_IMPL"/bin:$PATH
    export LD_LIBRARY_PATH=/usr/lib64/"$MPI_IMPL"/lib:$LD_LIBRARY_PATH
fi
if [[ "$BLAS_ENV" == "internal" ]]; then
    export USE_INTERNALBLAS=1
    export BLAS_SIZE=8
elif [[ "$BLAS_ENV" == "build_openblas" ]]; then
    export BUILD_OPENBLAS="y"
    export BLAS_SIZE=8
fi
if [[ -z "$USE_INTERNALBLAS" ]]; then
    if [[ -z "$BLASOPT" ]] ; then
	export BUILD_OPENBLAS="y"
	export BLAS_SIZE=8
    else
	unset BUILD_OPENBLAS
    fi
    if [[ "$SCALAPACK_ENV" == "off" ]]; then
	unset BUILD_SCALAPACK
	unset SCALAPACK
	unset SCALAPACK_SIZE
    else
	if [[ -z "$SCALAPACK" ]] ; then
	    export BUILD_SCALAPACK="y"
	    export SCALAPACK_SIZE=8
	else
	    unset BUILD_SCALAPACK
	fi
    fi
fi
export NWCHEM_EXECUTABLE=$TRAVIS_BUILD_DIR/.cachedir/binaries/$NWCHEM_TARGET/nwchem_"$arch"_`echo $NWCHEM_MODULES|sed 's/ /-/g'`_"$MPI_IMPL"
