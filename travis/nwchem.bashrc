#- env: == default ==
echo "start nwchem.bashrc"
echo "BLAS_SIZE is " $BLAS_SIZE
echo "BLASOPT is " $BLASOPT
echo "BUILD_OPENBLAS is " $BUILD_OPENBLAS
os=`uname`
arch=`uname -m`
if test -f "/usr/lib/os-release"; then
    dist=$(grep ID= /etc/os-release |head -1 |cut -c4-| sed 's/\"//g')
fi
if [ -z "$CC" ] ; then
    CC=cc
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
    if [[ "USE_AOMP" == "Y" ]]; then
	aomp_major=13
	aomp_minor=0-2
	export PATH=/usr/lib/aomp_"$aomp_major"."$aomp_minor"/bin/:$PATH
	export LD_LIBRARY_PATH=/usr/lib/aomp_"$aomp_major"."$aomp_minor"/lib:$LD_LIBRARY_PATH
    else
	source ${NWCHEM_TOP}/setenv_AOCC.sh
    fi
     export BUILD_MPICH=1
fi
if [[ "$FC" == "amdflang" ]]; then
    rocm_version=5.1.1
    export PATH=/opt/rocm-"$rocm_version"/bin:$PATH
    export LD_LIBRARY_PATH=/opt/rocm-"$rocm_version"/lib:/opt/rocm-"$rocm_version"/llvm/lib:$LD_LIBRARY_PATH
fi

if [[ "$FC" == "nvfortran" ]]; then
#    source /etc/profile.d/lmod.sh
#    module use /opt/nvidia/hpc_sdk/modulefiles
#    module load nvhpc
#     export BUILD_MPICH=1
     nv_major=22
     nv_minor=3
     nverdot="$nv_major"."$nv_minor"
     export PATH=/opt/nvidia/hpc_sdk/Linux_"$arch"/"$nverdot"/compilers/bin:$PATH
     export LD_LIBRARY_PATH=/opt/nvidia/hpc_sdk/Linux_"$arch"/"$nverdot"/compilers/lib:$LD_LIBRARY_PATH
     sudo /opt/nvidia/hpc_sdk/Linux_"$arch"/"$nverdot"/compilers/bin/makelocalrc -x
     export FC=nvfortran
     export MPICH_FC=nvfortran
#	if [ -z "$BUILD_MPICH" ] ; then
##use bundled openmpi
#	export PATH=/opt/nvidia/hpc_sdk/Linux_"$arch"/"$nverdot"/comm_libs/mpi/bin:$PATH
#	export LD_LIBRARY_PATH=/opt/nvidia/hpc_sdk/Linux_"$arch"/"$nverdot"/comm_libs/mpi/lib:$LD_LIBRARY_PATH
#	fi
#    export CC=gcc
fi
if [[ "$FC" == "ifort" ]] || [[ "$FC" == "ifx" ]] ; then
    case "$os" in
	Darwin)
	    IONEAPI_ROOT=~/apps/oneapi
	    ;;
	Linux)	
	    IONEAPI_ROOT=/opt/intel/oneapi
	    ;;		
    esac			
    source "$IONEAPI_ROOT"/compiler/latest/env/vars.sh
    export I_MPI_F90="$FC"
    "$FC" -V
    if [ -f "$IONEAPI_ROOT"/mkl/latest/env/vars.sh ] ; then
	source "$IONEAPI_ROOT"/mkl/latest/env/vars.sh
    fi

fi
if [[ "$MPI_IMPL" == "intel" ]]; then
    source /opt/intel/oneapi/mpi/latest/env/vars.sh
    export I_MPI_F90="$FC"
    mpif90 -v
    mpif90 -show
    if [ -f /opt/intel/oneapi/mkl/latest/env/vars.sh ] ; then
	source /opt/intel/oneapi/mkl/latest/env/vars.sh
    fi
fi
if [[ "$os" == "Darwin" ]]; then 
  export NWCHEM_TARGET=MACX64 
  export DYLD_FALLBACK_LIBRARY_PATH=$TRAVIS_BUILD_DIR/lib:$DYLD_FALLBACK_LIBRARY_PATH 
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
  if [[ "$MPI_IMPL" == "mpich" ]]; then
      export MPICH_FC=$FC
      if [[ "$arch" != "aarch64" ]]; then
	  export BUILD_MPICH=1
      fi
  fi
fi
export OMP_NUM_THREADS=1
export USE_NOIO=1

if [[ "$DISTR" == "fedora" ]] || [[ "$DISTR" == "centos" ]]; then
    export PATH=/usr/lib64/"$MPI_IMPL"/bin:$PATH
    export LD_LIBRARY_PATH=/usr/lib64/"$MPI_IMPL"/lib:$LD_LIBRARY_PATH
fi
if [[ -z "$BLAS_SIZE" ]] ; then
    echo "BLAS_SIZE not set, setting = 8"
    export BLAS_SIZE=8
fi

if [[ "$BLAS_ENV" == "internal" ]]; then
    export USE_INTERNALBLAS=1
    export SCALAPACK_ENV="off"
elif [[ "$BLAS_ENV" == "build_openblas" ]]; then
    export BUILD_OPENBLAS="y"
elif [[ "$BLAS_ENV" == "accelerate" ]]; then
    export BLASOPT="-framework Accelerate"
    export BLAS_LIB=${BLASOPT}
    export LAPACK_LIB=${BLASOPT}
    export BLAS_SIZE=4
fi
if [[ "$BLAS_SIZE" == "4" ]]; then
  export SCALAPACK_SIZE=4
  export USE_64TO32=y
fi
if [[ -z "$USE_INTERNALBLAS" ]]; then
    if [[ "$SCALAPACK_ENV" == "off" ]]; then
	unset BUILD_SCALAPACK
	unset SCALAPACK
	unset SCALAPACK_SIZE
    else
	if [[ -z "$SCALAPACK" ]] ; then
	    export BUILD_SCALAPACK="y"
	    if [[ -z "$SCALAPACK_SIZE" ]] ; then
		export SCALAPACK_SIZE=8
	    fi
#elpa
	    GFORTRAN_EXTRA=$(echo $FC | cut -c 1-8)
	    if  [[ ${FC} == gfortran ]] || [[ ${GFORTRAN_EXTRA} == gfortran ]] ; then
#	    if  [[ ${FC} == gfortran ]]  ; then
		if [[ `${CC} -dM -E - < /dev/null 2> /dev/null | grep -c clang` == 0 ]] && [[ `${CC} -dM -E - < /dev/null 2> /dev/null | grep -c GNU` > 0 ]] && [[ "$(expr `${CC} -dumpversion | cut -f1 -d.` \> 7)" == 1 ]]; then
		    if [[ "$os" == "Linux" ]] && [[ "$arch" == "x86_64" ]]; then
			if [[ ! -z "$BUILD_OPENBLAS" ]]; then
#			    if [[ "$SCALAPACK_SIZE" == "8" ]]; then
				export BUILD_ELPA=1
#			    fi
			fi
		    fi
		fi
	    fi
	else
	    unset BUILD_SCALAPACK
	fi
    fi
fi
#summary
echo "from nwchem.bashrc"
echo "BLAS_SIZE = " "$BLAS_SIZE"
echo "SCALAPACK_SIZE = " "$SCALAPACK_SIZE"
if [[ ! -z "$USE_64TO32" ]]; then
    echo "USE_64TO32 = " "$USE_64TO32"
fi
if [[ ! -z "$BUILD_ELPA" ]]; then
echo "BUILD_ELPA = " "$BUILD_ELPA"
fi
export NWCHEM_EXECUTABLE=$TRAVIS_BUILD_DIR/.cachedir/binaries/$NWCHEM_TARGET/nwchem_"$arch"_`echo $NWCHEM_MODULES|sed 's/ /-/g'`_"$MPI_IMPL"
