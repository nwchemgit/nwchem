#!/usr/bin/env bash
declare -a LIST_MPIFWRAP
LIST_MPIFWRAP=( 'mpif90' 'mpifort' 'mpiifort' 'mpifc' 'mpixlf_r' 'mpixlf' 'mpif77' 'mpifrt' 'mpifrtpx' 'mpif90-mpich-mp' 'mpifort-mpich-mp' 'mpif90-openmpi-mp' 'mpifort-openmpi-mp')
export LIST_MPIFWRAP
declare -a LIST_MPICWRAP
LIST_MPICWRAP=( 'mpicc' 'mpicc'    'mpiicc'  'mpiicc' 'mpixlc_r' 'mpixlc' 'mpicc' 'mpifcc' 'mpifccpx' 'mpicc-mpich-mp' 'mpicc-mpich-mp' 'mpicc-openmpi-mp' 'mpicc-openmpi-mp')
export LIST_MPICWRAP
myscalapwd=`pwd`
source ../libext_utils/cmake.sh
cd $myscalapwd

if [[ -z "${MPIF90}" ]]; then
if [[ "$FC" = "ftn"  ]] ; then
    MPIF90="ftn"
    MPICC="cc"
else
    if ! [ -x "$(command -v mpif90)" ]; then
	length=${#LIST_MPIFWRAP[*]}
	indx=0
	while [ "${indx}" -lt "${length}" ] ; do
	    candidate="${LIST_MPIFWRAP[${indx}]}"
	    CMDGUESS="$(command -v $candidate)"
	    if [ $? -eq 0 ] ; then
		MPIF90="${LIST_MPIFWRAP[${indx}]}"
		MPICC="${LIST_MPICWRAP[${indx}]}"
		echo found mpi wrappers $MPIF90 $MPICC
		break
	    else   ((indx++)); fi
	done
	if [ -z ${MPIF90} ] ; then
	    echo
	    echo mpif90 not installed
	    echo mpif90 is required for building Scalapack
	    echo
	    exit 1
	fi
    else
	MPIF90="mpif90"
	if [[ -z "${MPICC}" ]]; then
            MPICC=mpicc
	fi
    fi
fi
fi
if [[  -z "${FC}" ]]; then
    FC=$($MPIF90 -show|cut -d " " -f 1)
fi
if [[  -z "${CC}" ]]; then
    CC=$($MPICC -show|cut -d " " -f 1)
fi
if [[  -z "${NWCHEM_TOP}" ]]; then
    dir3=$(dirname `pwd`)
    dir2=$(dirname "$dir3")
    NWCHEM_TOP=$(dirname "$dir2")
fi
# take care of xcode 15 quirks
source ${NWCHEM_TOP}/src/config/fix_xcode15.sh

if [[ ! -z "${BUILD_MPICH}" ]]; then
    export PATH=${NWCHEM_TOP}/src/libext/bin:$PATH
    if [ -x "$(command -v pkg-config1)" ]; then
        export LDFLAGS=`pkg-config --libs-only-L hwloc`
    else
	if [ -x "$(command -v brew)" ]; then
	    export LDFLAGS=-L`brew --prefix`/lib/
	else
	    echo 'WARNING: cannot guess the location of the hwloc library'
#	    exit 1
	fi
    fi
    echo LDFLAGS for hwloc is $LDFLAGS
fi
if [[ "$FC" = "ftn"  ]] || [[ ! -z "$USE_CMAKE_MASTER" ]] ; then
    get_cmake_master
else
if [[ -z "${CMAKE}" ]]; then
    #look for cmake
    if [[ -z "$(command -v cmake)" ]]; then
	cmake_instdir=../libext_utils
	get_cmake_release $cmake_instdir
	status=$?
	if [ $status -ne 0 ]; then
	    echo cmake required to build scalapack
	    echo Please install cmake
	    echo define the CMAKE env. variable
	    exit 1
	fi
    else
	CMAKE=cmake
    fi
fi
fi
CMAKE_VER_MAJ=$(${CMAKE} --version|cut -d " " -f 3|head -1|cut -d. -f1)
CMAKE_VER_MIN=$(${CMAKE} --version|cut -d " " -f 3|head -1|cut -d. -f2)
echo CMAKE_VER is ${CMAKE_VER_MAJ} ${CMAKE_VER_MIN}
if ((CMAKE_VER_MAJ < 3)) || (((CMAKE_VER_MAJ > 2) && (CMAKE_VER_MIN < 99))); then
    cmake_instdir=../libext_utils
    get_cmake_release $cmake_instdir
    status=$?
    if [ $status -ne 0 ]; then
	echo cmake required to build scalapack
	echo Please install cmake
	echo define the CMAKE env. variable
	exit 1
    fi
fi
cd $myscalapwd
pwd

#if [[ "$SCALAPACK_SIZE" != "4"  ]] ; then
#    echo SCALAPACK_SIZE must be equal to 4
#    exit 1
#fi
#if [[ "$BLAS_SIZE" != "4"  ]] ; then
#    echo BLAS_SIZE must be equal to 4 for SCALAPACK
#    exit 1
#fi
if [[ "$BLAS_SIZE" != "$SCALAPACK_SIZE"  ]] ; then
    echo "BLAS_SIZE must be the same as SCALAPACK_SIZE"
    echo "BLAS_SIZE = " "$BLAS_SIZE"
    echo "SCALAPACK_SIZE = " "$SCALAPACK_SIZE"
    exit 1
fi

if [[  -z "${SCALAPACK_SIZE}" ]]; then
   SCALAPACK_SIZE=8
fi
if [[ "$BLAS_SIZE" == 4 ]] && [[ -z "$USE_64TO32"   ]] ; then
    if [[ "$NWCHEM_TARGET" != "LINUX" ]] && [[ "$NWCHEM_TARGET" != "MACX" ]] ; then
    echo USE_64TO32 must be set when BLAS_SIZE=4 on 64-bit architectures
    exit 1
    fi
fi
if [[ ! -z "$BUILD_OPENBLAS"   ]] ; then
    BLASOPT="-L`pwd`/../lib -lnwc_openblas"
fi
#git clone https://github.com/scibuilder/scalapack.git
#svn co --non-interactive --trust-server-cert https://icl.utk.edu/svn/scalapack-dev/scalapack/trunk/ scalapack
VERSION=2.1.0
#curl -L https://github.com/Reference-ScaLAPACK/scalapack/archive/v${VERSION}.tar.gz -o scalapack.tgz
#COMMIT=bc6cad585362aa58e05186bb85d4b619080c45a9
#COMMIT=ea5d20668a6b8bbee645b7ffe44623c623969d33
COMMIT=5bad7487f496c811192334640ce4d3fc5f88144b
COMMIT=782e739f8eb0e7f4d51ad7dd23fc1d03dc99d240
COMMIT=a23c2cdc6586c427686f6097ae66bb54ef693571
rm -rf scalapack 
if [[ -f "scalapack-$COMMIT.tar.gz" ]]; then
    echo "using existing"  "scalapack-$COMMIT.tar.gz"
else
    echo "downloading"  "scalapack-$COMMIT.tar.gz"
    rm -f scalapack-$COMMIT.tar.gz
    tries=1 ; until [ "$tries" -ge 6 ] ; do
		  if [ "$tries" -gt 1 ]; then sleep 9; echo attempt no.  $tries ; fi
		  curl -L https://github.com/Reference-ScaLAPACK/scalapack/archive/$COMMIT.tar.gz -o scalapack-$COMMIT.tar.gz
		  # check tar.gz integrity
		  gzip -t scalapack-$COMMIT.tar.gz >&  /dev/null
		  if [ $? -eq 0 ]; then break ;  fi
		  tries=$((tries+1)) ;  done
fi
tar xzf scalapack-$COMMIT.tar.gz
ln -sf scalapack-*$COMMIT scalapack
#ln -sf scalapack-${VERSION} scalapack
#curl -L http://www.netlib.org/scalapack/scalapack-${VERSION}.tgz -o scalapack.tgz
#tar xzf scalapack.tgz
cd scalapack
# macos accelerate does not contain dcombossq
if [[  ! -z "$USE_INTERNALBLAS" ]]; then
    export USE_DCOMBSSQ=1
    BLASOPT="-L${NWCHEM_TOP}/lib/${NWCHEM_TARGET} -lnwclapack -lnwcblas"
fi
if [[ $(echo "$LAPACK_LIB" |awk '/Accelerate/ {print "Y"; exit}' ) == "Y" ]]; then
    export USE_DCOMBSSQ=1
fi
if [[ $(echo "$LAPACK_LIB" |awk '/lapack/ {print "Y"; exit}' ) == "Y" ]]; then
    export USE_DCOMBSSQ=1
fi
if [[ $(echo "$LAPACK_LIB" |awk '/lfjlapack/ {print "Y"; exit}'  ) == "Y" ]]; then
    export USE_DCOMBSSQ=1
fi
# ubuntu 24.04 openblas pkg does not contain dcombossq
if [[ $(echo "$LAPACK_LIB" |awk '/openblas/ {print "Y"; exit}'  ) == "Y" ]]; then
    export USE_DCOMBSSQ=1
fi
if [[  -z "$USE_DCOMBSSQ" ]]; then
    patch -p0 -s -N < ../dcombssq.patch
fi
patch -p0 -s -N < ../cmake.patch
#curl -LJO https://github.com/Reference-ScaLAPACK/scalapack/commit/189c84001bcd564296a475c5c757afc9f337e828.patch
#patch -p1 < 189c84001bcd564296a475c5c757afc9f337e828.patch
rm -rf build
mkdir -p build
cd build
if  [[ -n ${FC} ]] &&   [[ ${FC} == xlf ]] || [[ ${FC} == xlf_r ]] || [[ ${FC} == xlf90 ]]|| [[ ${FC} == xlf90_r ]]; then
    Fortran_FLAGS=" -qintsize=4 -qextname "
elif [[ -n ${FC} ]] &&   [[ ${FC} == ftn ]]; then
    if [[ ${PE_ENV} == INTEL ]]; then
	Fortran_FLAGS="-O2 -g -axCORE-AVX2"
    fi
#elif [[ -n ${FC} ]] &&   [[ ${FC} == flang ]]; then
# unset FC=flang since cmake gets lost
#       unset FC
fi
#if [[ ! -z "$BUILD_SCALAPACK"   ]] ; then
#    Fortran_FLAGS+=-I"$NWCHEM_TOP"/src/libext/include
#fi
echo "SCALAPACK_SIZE" is $SCALAPACK_SIZE
if [[ ${FC} == ftn ]]; then
    if [[ ${PE_ENV} == PGI ]]; then
          FC=pgf90
    fi
    if [[ ${PE_ENV} == INTEL ]]; then
	FC=ifort
    fi
    if [[ ${PE_ENV} == GNU ]]; then
	FC=gfortran
    fi
    if [[ ${PE_ENV} == AOCC ]]; then
	FC=flang
    fi
    if [[ ${PE_ENV} == NVIDIA ]]; then
	FC=nvfortran
    fi
    if [[ ${PE_ENV} == CRAY ]]; then
	FC=crayftn
	CC=clang
	#fix for libunwind.so link problem
        export LD_LIBRARY_PATH=/opt/cray/pe/cce/$CRAY_FTN_VERSION/cce-clang/x86_64/lib:/opt/cray/pe/lib64/cce/:$LD_LIBRARY_PATH
    fi
fi
FC_EXTRA=$(${NWCHEM_TOP}/src/config/strip_compiler.sh ${FC})

if [[  -z "$MPICH_FC"   ]] ; then
    export MPICH_FC="$FC"
fi
echo MPICH_FC is "$MPICH_FC"
if [[  -z "$MPICH_CC"   ]] ; then
    export MPICH_CC="$CC"
fi
echo MPICH_CC is "$MPICH_CC"
echo $(${MPICC} -show)
#Intel MPI
if [[  -z "$I_MPI_F90"   ]] ; then
    export I_MPI_F90="$FC"
fi
if [[  -z "$I_MPI_CC"   ]] ; then
    export I_MPI_CC="$CC"
fi
echo I_MPI_F90 is "$I_MPI_F90"
echo I_MPI_CC is "$I_MPI_CC"
if [[  -z "$PE_ENV"   ]] ; then
    #check if mpif90 and FC are consistent
    MPIF90_EXTRA=$(${NWCHEM_TOP}/src/config/strip_compiler.sh `${MPIF90} -show`)
    if [[ $MPIF90_EXTRA != $FC_EXTRA ]]; then
        echo which mpif90 is `which mpif90`
        echo mpif90show `${MPIF90} -show`
	echo FC and MPIF90 are not consistent
	echo FC is $FC_EXTRA
	echo MPIF90 is $MPIF90_EXTRA
	exit 1
    fi
fi
#fix for clang 12 error in implicit-function-declaration
GOTCLANG=$( "$MPICC" -dM -E - </dev/null 2> /dev/null |grep __clang__|head -1|cut -c19)
if [[ `${MPICC} -dM -E - < /dev/null 2> /dev/null | grep -c GNU` > 0 ]] ; then
    let GCCVERSIONGT12=$(expr `${MPICC} -dumpversion | cut -f1 -d.` \> 12)
    let GCCVERSIONGT13=$(expr `${MPICC} -dumpversion | cut -f1 -d.` \> 13)
    let GCCVERSIONGT14=$(expr `${MPICC} -dumpversion | cut -f1 -d.` \> 14)
fi

if [[ ${GOTCLANG} == "1" ]] || [[ ${GCCVERSIONGT12} == "1" ]]  ; then
    C_FLAGS=" -Wno-error=implicit-function-declaration "
fi
if [[ ${GCCVERSIONGT13} == "1" ]]  ; then
    C_FLAGS+=" -std=gnu17 "
fi
if [[ ${GCCVERSIONGT14} == "1" ]]  ; then
    C_FLAGS+=" -fPIE "
    Fortran_FLAGS+=" -fPIE "
fi
if [[  "$SCALAPACK_SIZE" == 8 ]] ; then
    if  [[ ${FC} == f95 ]] || [[ ${FC_EXTRA} == gfortran ]] ; then
    Fortran_FLAGS+=" -fdefault-integer-8 -w "
    elif  [[ ${FC_EXTRA} == flang ]] ; then
    Fortran_FLAGS+=" -fdefault-integer-8 "
    elif  [[ ${FC} == xlf ]] || [[ ${FC} == xlf_r ]] || [[ ${FC} == xlf90 ]]|| [[ ${FC} == xlf90_r ]]; then
    Fortran_FLAGS=" -qintsize=8 -qextname "
    elif  [[ ${FC} == crayftn ]]; then
    Fortran_FLAGS=" -s integer64 -h nopattern"
    elif  [[ ${FC} == frtpx ]] || [[ ${FC} == frt ]]; then
    Fortran_FLAGS=" -fs -CcdLL8 -CcdII8 "
    else
    Fortran_FLAGS+=" -i8 "
    fi
    C_FLAGS+=" -DInt=long"
fi
#cross-compilation: we set CDEFS
#https://github.com/Reference-ScaLAPACK/scalapack/commit/1bdf63ec17bf8e827b8c5abd292f0e41bdc2f56e
CMAKE_EXTRA=" "
if  [[ ${FC} == frtpx ]] ||  [ -x "$(command -v xx-info)" ]; then
    CMAKE_EXTRA="-DCDEFS=Add_"
fi    
#skip argument check for gfortran
arch=`uname -m`
echo arch is $arch
if  [[ ${FC_EXTRA} == nvfortran ]]; then
echo 'nvfortran -V is ' `nvfortran -V`
    if  [[ ${USE_HWOPT} == n ]]; then
      if [[ "$arch" == "x86_64" ]]; then
	Fortran_FLAGS+=" -tp px "
      fi
    fi
fi
if  [[ ${FC_EXTRA} == flang ]]; then
    Fortran_FLAGS+=" -fPIC "
fi
if  [[ ${FC_EXTRA} == gfortran ]] || [[ ${FC} == f95 ]]; then
    Fortran_FLAGS+=" -fPIC "
    if [[ "$(expr `${FC} -dumpversion | cut -f1 -d.` \> 7)" == 1 ]]; then
	Fortran_FLAGS+=" -std=legacy "
    fi
    LDFLAGS+=" -fno-lto "
fi
if [[ ${PE_ENV} == NVIDIA ]] || [[ ${FC} == nvfortran ]] ; then
  Fortran_FLAGS+=" -fPIC "
fi
if [[ "$CRAY_CPU_TARGET" == "mic-knl" ]]; then
    module swap craype-mic-knl craype-haswell
    KNL_SWAP=1
fi

# force -m32 flag on 32-bit x86 linux to avoid -mx32
Fortran_FLAGS_RELWITHDEB=" -O2 -g -DNDEBUG "
if [[ "$arch" == "i686" ]] || [[ "$arch" == "x86_64" ]]; then
    if [[ ${NWCHEM_TARGET} == LINUX ]] && [[ ${FC_EXTRA} == gfortran ]] ; then
       Fortran_FLAGS+=" -m32 "
       C_FLAGS+=" -m32 "
    fi
fi
echo " $MPIF90 -show is " `$MPIF90 -show`
echo LDFLAGS is $LDFLAGS
if [[ ${FC} == nvfortran ]] ; then
    echo compiling with CC="$MPICC"  FC=$FC MPIF90=$MPIF90 CFLAGS="$C_FLAGS" FFLAGS="$Fortran_FLAGS" $CMAKE -Wno-dev ../ -DCMAKE_BUILD_TYPE=MinSizeRel -DCMAKE_C_FLAGS="$C_FLAGS"  -DCMAKE_Fortran_FLAGS="$Fortran_FLAGS" -DTEST_SCALAPACK=OFF  -DBUILD_TESTING=OFF -DBUILD_SHARED_LIBS=OFF  -DBLAS_openblas_LIBRARY="$BLASOPT"  -DBLAS_LIBRARIES="$BLASOPT"  -DLAPACK_openblas_LIBRARY="$BLASOPT"  -DLAPACK_LIBRARIES="$BLASOPT" -DCMAKE_Fortran_FLAGS_RELWITHDEBINFO="-O2 -g -DNDEBUG  $Fortran_FLAGS"  $CMAKE_EXTRA  -DMPI_Fortran_COMPILE_OPTIONS="$Fortran_FLAGS"
    CC="$MPICC"  FC=$FC MPIF90=$MPIF90 CFLAGS="$C_FLAGS" FFLAGS="$Fortran_FLAGS" $CMAKE -Wdev ../ -DCMAKE_BUILD_TYPE=MinSizeRel -DCMAKE_C_FLAGS="$C_FLAGS"  -DCMAKE_Fortran_FLAGS="$Fortran_FLAGS" -DTEST_SCALAPACK=OFF  -DBUILD_TESTING=OFF -DBUILD_SHARED_LIBS=OFF  -DBLAS_openblas_LIBRARY="$BLASOPT"  -DBLAS_LIBRARIES="$BLASOPT"  -DLAPACK_openblas_LIBRARY="$BLASOPT"  -DLAPACK_LIBRARIES="$BLASOPT" -DCMAKE_Fortran_FLAGS_RELWITHDEBINFO="-O2 -g -DNDEBUG  $Fortran_FLAGS" $CMAKE_EXTRA -DMPI_Fortran_COMPILE_OPTIONS="$Fortran_FLAGS"
else
    echo compiling with CC="$MPICC"  CFLAGS="$C_FLAGS" FC=$MPIF90 FFLAGS="$Fortran_FLAGS" $CMAKE -Wno-dev ../ -DCMAKE_BUILD_TYPE=MinSizeRel -DCMAKE_C_FLAGS="$C_FLAGS"  -DCMAKE_Fortran_FLAGS="$Fortran_FLAGS" -DTEST_SCALAPACK=OFF  -DBUILD_TESTING=OFF -DBUILD_SHARED_LIBS=OFF  -DBLAS_openblas_LIBRARY="$BLASOPT"  -DBLAS_LIBRARIES="$BLASOPT"  -DLAPACK_openblas_LIBRARY="$BLASOPT"  -DLAPACK_LIBRARIES="$BLASOPT" -DCMAKE_Fortran_FLAGS_RELWITHDEBINFO="-O2 -g -DNDEBUG  $Fortran_FLAGS"  $CMAKE_EXTRA
    CC="$MPICC"  FC=$MPIF90 CFLAGS="$C_FLAGS" FFLAGS="$Fortran_FLAGS" $CMAKE -Wdev ../ -DCMAKE_BUILD_TYPE=MinSizeRel -DCMAKE_C_FLAGS="$C_FLAGS"  -DCMAKE_Fortran_FLAGS="$Fortran_FLAGS" -DTEST_SCALAPACK=OFF  -DBUILD_TESTING=OFF -DBUILD_SHARED_LIBS=OFF  -DBLAS_openblas_LIBRARY="$BLASOPT"  -DBLAS_LIBRARIES="$BLASOPT"  -DLAPACK_openblas_LIBRARY="$BLASOPT"  -DLAPACK_LIBRARIES="$BLASOPT" -DCMAKE_Fortran_FLAGS_RELWITHDEBINFO="-O2 -g -DNDEBUG  $Fortran_FLAGS" $CMAKE_EXTRA
fi
if [[ "$?" != "0" ]]; then
    echo " "
    echo "cmake failed"
    echo " "
    cat $(find . -name *log)
    exit 1
fi
make V=0 -j3 scalapack/fast
if [[ "$?" != "0" ]]; then
    echo " "
    echo "compilation failed"
    echo " "
    exit 1
fi
mkdir -p ../../../lib
if [[ $(uname -s) == "Linux" ]]; then
    if [ -x "$(command -v xx-info)" ]; then
	MYSTRIP=$(xx-info)-strip
    else
	MYSTRIP=strip
    fi
    echo MYSTRIP is $MYSTRIP
    $MYSTRIP --strip-debug lib/libscalapack.a
fi
cp lib/libscalapack.a ../../../lib/libnwc_scalapack.a
if [[ "$KNL_SWAP" == "1" ]]; then
    module swap  craype-haswell craype-mic-knl
fi
