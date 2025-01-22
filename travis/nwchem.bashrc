#- env: == default ==
echo "start nwchem.bashrc"
echo "BLAS_SIZE is " $BLAS_SIZE
echo "BLASOPT is " $BLASOPT
echo "BUILD_OPENBLAS is " $BUILD_OPENBLAS
os=`uname`
arch=`uname -m`
if [[ -z "$APPTAINER_NAME" ]] || [[ -z "$SINGULARITY_NAME" ]] ; then
    MYSUDO=sudo
else
    MYSUDO=" "
fi
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
	aomp_major=18
	aomp_minor=0-0
	export PATH=/usr/lib/aomp_"$aomp_major"."$aomp_minor"/bin/:$PATH
	export LD_LIBRARY_PATH=/usr/lib/aomp_"$aomp_major"."$aomp_minor"/lib:$LD_LIBRARY_PATH
    else
	source ${NWCHEM_TOP}/setenv_AOCC.sh
    fi
     export BUILD_MPICH=1
fi
if [[ "$FC" == "amdflang" ]]; then
    export PATH=/opt/rocm/bin:$PATH
    export LD_LIBRARY_PATH=/opt/rocm-"$rocm_version"/lib:/opt/rocm/llvm/lib:$LD_LIBRARY_PATH
     export BUILD_MPICH=1
fi
if [[ "$FC" == 'flang-new-'* ]]; then
    export BUILD_MPICH=1
fi

if [[ "$FC" == "nvfortran" ]]; then
     export PATH=/opt/nvidia/hpc_sdk/Linux_"$arch"/latest/compilers/bin:$PATH
     echo "**** before hpc_sdk makelocalrc"
     $MYSUDO /opt/nvidia/hpc_sdk/Linux_"$arch"/latest/compilers/bin/makelocalrc -x
     echo "**** done hpc_sdk makelocalrc"
     export LD_LIBRARY_PATH=/opt/nvidia/hpc_sdk/Linux_"$arch"/latest/compilers/lib:$LD_LIBRARY_PATH
     export FC=nvfortran
     export MPICH_FC=nvfortran
fi
if [[ "$FC" == "ifort" ]] || [[ "$FC" == "ifx" ]] ; then
    if [[ "$os" == "Darwin" ]]; then
	IONEAPI_ROOT=~/apps/oneapi
    else
	IONEAPI_ROOT=/opt/intel/oneapi
# fix runtime mpi_init error
#	export FI_LOG_LEVEL=TRACE
	export FI_PROVIDER=shm
	echo "*** output of fi_info ***"
	echo $(fi_info -l) || true
    fi
#    source "$IONEAPI_ROOT"/compiler/latest/env/vars.sh
    source "$IONEAPI_ROOT"/setvars.sh --force
    export I_MPI_F90="$FC"
#force icc on macos to cross-compile x86 on arm64
# icx not available on macos
    if [[ "$os" == "Darwin" ]]; then
     if [[ "$FC" != "gfortran" ]] && [[ "$FC" == "gfortran*" ]]; then
	 mygccver=$(echo "$FC"|cut -d - -f 2)
	 export PATH=$HOMEBREW_CELLAR/../opt/gcc@"$mygccver"/bin:$PATH
	 echo gfortran is $(gfortran -v)
	 echo gfortran-"$mygccver" is $(gfortran-"$mygccver" -v)
     fi
     if [[ "$CC" != gcc ]] && [[ "$CC" == gcc* ]]; then
	 mygccver=$(echo "$CC"|cut -d - -f 2)
	 export PATH=$HOMEBREW_CELLAR/../opt/gcc@"$mygccver"/bin:$PATH
	 echo gcc is $(gcc -v)
	 echo gcc-"$mygccver" is $(gcc-"$mygccver" -v)
     fi
#	CC=icc
#	CXX=icc
# Intel MPI not available on macos       
#       export BUILD_MPICH=1
       unset BUILD_PLUMED
       
    fi
    "$FC" -v
    "$CC" -v
fi
if [[ -f "$IONEAPI_ROOT"/mpi/latest/env/vars.sh ]]; then
    source "$IONEAPI_ROOT"/mpi/latest/env/vars.sh
    export I_MPI_F90="$FC"
    if [[ "$MPI_IMPL" == "intel" ]]; then
	mpif90 -v
	mpif90 -show
    fi
    if [ -f /opt/intel/oneapi/mkl/latest/env/vars.sh ] ; then
	source /opt/intel/oneapi/mkl/latest/env/vars.sh
    fi
fi
if [[ "$MPI_IMPL" == "build_mpich" ]]; then 
  export BUILD_MPICH=1
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

fi
if [[ "$os" == "Linux" ]]; then 
   export NWCHEM_TARGET=LINUX64 
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
elif [[ "$BLAS_ENV" == lib*openblas* ]] || [[ "$BLAS_ENV" == "brew_openblas" ]]; then
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
		if [[ "$MPI_IMPL" == "mpich" ]]; then
		    export MPICH_FC=$FC
		    if [[ "$arch" != "aarch64" ]]; then
			export BUILD_MPICH=1
		    fi
		fi
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
export NWCHEM_EXECUTABLE=$TRAVIS_BUILD_DIR/.cachedir/binaries/$NWCHEM_TARGET/nwchem_"$arch"_`echo $NWCHEM_MODULES|sed 's/ /-/g'`_"$MPI_IMPL"_"$FC"

if [[ "$FC" == "gfortran" ]]; then
   if [[ "$($FC -dM -E - < /dev/null 2> /dev/null | grep __GNUC__ |cut -c 18-)" -lt 9 ]]; then
#disable xtb  if gfortran version < 9
     unset USE_TBLITE
     export NWCHEM_MODULES=$(echo $NWCHEM_MODULES |sed  's/xtb//')
   fi
fi

if [[ "$USE_LIBXC" == "-1" ]]; then
    unset USE_LIBXC
    export LIBXC_LIB=/usr/lib/x86_64-linux-gnu
    export LIBXC_INCLUDE=/usr/include
fi

#export PYTHONVERSION=$(python -V| cut -d ' ' -f 2 |cut -d . -f 1-2)
