#!/usr/bin/env bash
os=`uname`
dist="ubuntu"
arch=`uname -m`
env | grep FC || true
env | grep CC || true
if test -f "/usr/lib/os-release"; then
    dist=$(grep ID= /etc/os-release |head -1 |cut -c4-| sed 's/\"//g')
fi
if test -f "/usr/lib/fedora-release"; then
    dist="fedora"
fi
if test -f "/usr/lib/centos-release"; then
    dist="centos"
fi
echo dist is "$dist"
if [ -z "$DISTR" ] ; then
    DISTR=$dist
fi
echo DISTR is "$DISTR"
case "$os" in
    Darwin)
	IONEAPI_ROOT=~/apps/oneapi
	;;
    Linux)
	IONEAPI_ROOT=/opt/intel/oneapi
	;;
esac
 if [[ "$os" == "Darwin" ]]; then 
#  HOMEBREW_NO_AUTO_UPDATE=1 brew cask uninstall oclint || true  
#  HOMEBREW_NO_INSTALL_CLEANUP=1  HOMEBREW_NO_AUTO_UPDATE=1 brew install gcc "$MPI_IMPL" openblas python3 ||true
     HOMEBREW_NO_INSTALL_CLEANUP=1  HOMEBREW_NO_AUTO_UPDATE=1 brew install gcc "$MPI_IMPL" python3 gsed grep automake autoconf ||true
     if [[ "$FC" == "ifort" ]] || [[ "$FC" == "ifx" ]] ; then
         if [[ -f ~/apps/oneapi/setvars.sh ]]; then 
	     echo ' using intel cache installation '
	 else
	mkdir -p ~/mntdmg ~/apps/oneapi || true
	cd ~/Downloads
	dir_base="18342"
	dir_hpc="18341"
	base="m_BaseKit_p_2022.1.0.92_offline"
	hpc="m_HPCKit_p_2022.1.0.86_offline"
	curl -LJO https://registrationcenter-download.intel.com/akdlm/irc_nas/"$dir_base"/"$base".dmg
	curl -LJO https://registrationcenter-download.intel.com/akdlm/irc_nas/"$dir_hpc"/"$hpc".dmg
	echo "installing BaseKit"
	hdiutil attach "$base".dmg  -mountpoint ~/mntdmg -nobrowse
	sudo ~/mntdmg/bootstrapper.app/Contents/MacOS/install.sh --cli  --eula accept \
	     --action install --components default  --install-dir ~/apps/oneapi
	hdiutil detach ~/mntdmg
	#
	echo "installing HPCKit"
	hdiutil attach "$hpc".dmg  -mountpoint ~/mntdmg -nobrowse
	sudo ~/mntdmg/bootstrapper.app/Contents/MacOS/install.sh --cli  --eula accept \
	     --action install --components default --install-dir ~/apps/oneapi
	hdiutil detach ~/mntdmg
	ls -lrta ~/apps/oneapi ||true
	sudo rm -rf "$IONEAPI_ROOT"/intelpython "$IONEAPI_ROOT"/dal "$IONEAPI_ROOT"/advisor \
	     "$IONEAPI_ROOT"/ipp "$IONEAPI_ROOT"/conda_channel 	"$IONEAPI_ROOT"/dnnl \
	     "$IONEAPI_ROOT"/installer "$IONEAPI_ROOT"/vtune_profiler "$IONEAPI_ROOT"/tbb || true
	fi
	 source "$IONEAPI_ROOT"/setvars.sh || true
	 export I_MPI_F90="$FC"
	ls -lrta ~/apps/oneapi ||true
	df -h 
	rm -f *dmg || true
	df -h
	"$FC" -V
	icc -V
     else
	 #hack to fix Github actions mpif90
	 gccver=`brew list --versions gcc| head -1 |cut -c 5-`
	 echo brew gccver is $gccver
	 ln -sf /usr/local/Cellar/gcc/$gccver/bin/gfortran-* /usr/local/Cellar/gcc/$gccver/bin/gfortran || true
	 ln -sf /usr/local/Cellar/gcc/$gccver/bin/gfortran-* /usr/local/bin/gfortran || true
#	 ln -sf /usr/local/bin/$FC /usr/local/bin/gfortran
	 $FC --version
	 gfortran --version
     fi
     #hack to get 3.10 as default
     brew install python@3.10
     brew link --force --overwrite python@3.10
#  if [[ "$MPI_IMPL" == "openmpi" ]]; then
#      HOMEBREW_NO_INSTALL_CLEANUP=1 HOMEBREW_NO_AUTO_UPDATE=1 brew install scalapack
#  fi
fi
 if [[ "$os" == "Linux" ]]; then
     if [[ "$DISTR" == "fedora" ]] || [[ "$DISTR" == "centos" ]] ; then
	 rpminst=dnf
	 if [[ "$DISTR" == "centos" ]] ; then
	     rpminst=yum
	 fi
	 sudo $rpminst udate;  sudo $rpminst -y install perl perl python3-devel time patch openblas-serial64 openmpi-devel cmake gcc-gfortran unzip which make tar bzip2 openssh-clients rsync
	 #	 module load mpi
	 if [[ "$MPI_IMPL" == "openmpi" ]]; then
	     sudo $rpminst -y install  openmpi-devel
	 else
	     echo ready only for openmpi
	     exit 1
	 fi
	 export PATH=/usr/lib64/"$MPI_IMPL"/bin:$PATH
	 export LD_LIBRARY_PATH=/usr/lib64/"$MPI_IMPL"/lib:$LD_LIBRARY_PATH
	 which mpif90
	 mpif90 -show
     else
    if [[ "$MPI_IMPL" == "openmpi" ]]; then
	mpi_bin="openmpi-bin" ; mpi_libdev="libopenmpi-dev" scalapack_libdev="libscalapack-openmpi-dev"
    fi
    if [[ "$MPI_IMPL" == "mpich" ]]; then
        mpi_bin="mpich" ; mpi_libdev="libmpich-dev" scalapack_libdev="libscalapack-mpich-dev"
    fi
    if [[ "$MPI_IMPL" == "intel" || "$FC" == "ifort" || "$FC" == "ifx" ]]; then
	export APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=1
        tries=0 ; until [ "$tries" -ge 10 ] ; do \
	wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
            && sudo -E apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB  \
            && rm -f GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB || true \
            && echo "deb https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list \
            && sudo add-apt-repository "deb https://apt.repos.intel.com/oneapi all main"  \
	    && sudo apt-get update && break ;\
            tries=$((tries+1)) ; echo attempt no.  $tries    ; sleep 30 ;  done

        mpi_bin="  " ; mpi_libdev="intel-oneapi-mpi-devel" scalapack_libdev="intel-oneapi-mkl"
    fi
    sudo apt-get update
    sudo apt-get -y install software-properties-common
    sudo add-apt-repository universe && sudo apt-get update
#    sudo apt-get -y install gfortran python3-dev python-dev cmake "$mpi_libdev" "$mpi_bin" "$scalapack_libdev"  make perl  libopenblas-dev python3 rsync
    sudo apt-get -y install gfortran python3-dev python-dev cmake "$mpi_libdev" "$mpi_bin"  make perl  python3 rsync
    if [[ "$FC" == "gfortran-11" ]] || [[ "$CC" == "gcc-11" ]]; then
	sudo  add-apt-repository -y ppa:ubuntu-toolchain-r/test 
        sudo  apt-get -y install gcc-11 gfortran-11 g++-11
    fi
    if [[ "$FC" == "ifort" ]] || [[ "$FC" == "ifx" ]]; then
	sudo apt-get -y install intel-oneapi-ifort intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic  intel-oneapi-mkl
	if [[ "$?" != 0 ]]; then
	    echo "apt-get install failed: exit code " "${?}"
	    exit 1
	fi
	sudo apt-get -y install intel-oneapi-mpi-devel
    fi
    if [[ "$FC" == "flang" ]]; then
	if [[ "USE_AOMP" == "Y" ]]; then
	    aomp_major=13
	    aomp_minor=0-6
	    wget https://github.com/ROCm-Developer-Tools/aomp/releases/download/rel_"$aomp_major"."$aomp_minor"/aomp_Ubuntu2004_"$aomp_major"."$aomp_minor"_amd64.deb
	    sudo dpkg -i aomp_Ubuntu2004_"$aomp_major"."$aomp_minor"_amd64.deb
	    export PATH=/usr/lib/aomp_"$aomp_major"."$aomp_minor"/bin/:$PATH
	    export LD_LIBRARY_PATH=/usr/lib/aomp_"$aomp_major"."$aomp_minor"/lib:$LD_LIBRARY_PATH
	    ls -lrt /usr/lib | grep aomp ||true
	else
	    aocc_version=3.2.0
	    aocc_dir=aocc-compiler-${aocc_version}
	    curl -LJO https://developer.amd.com/wordpress/media/files/${aocc_dir}.tar
	    tar xf ${aocc_dir}.tar
	    ./${aocc_dir}/install.sh
	    source setenv_AOCC.sh
	    pwd
	fi
	flang -v
	which flang
    fi
    if [[ "$FC" == "amdflang" ]]; then
	sudo apt-get install -y wget gnupg2 coreutils dialog tzdata
	rocm_version=4.5.2
	wget -q -O - https://repo.radeon.com/rocm/rocm.gpg.key |  sudo apt-key add -
	echo 'deb [arch=amd64] https://repo.radeon.com/rocm/apt/'$rocm_version'/ ubuntu main' | sudo tee /etc/apt/sources.list.d/rocm.list
	sudo apt-get  update -y && sudo apt-get -y install rocm-llvm openmp-extras
	export PATH=/opt/rocm-"$rocm_version"/bin:$PATH
	export LD_LIBRARY_PATH=/opt/rocm-"$rocm_version"/lib:/opt/rocm-"$rocm_version"/llvm/lib:$LD_LIBRARY_PATH
	amdflang -v
	amdclang -v
    fi
    if [[ "$FC" == "nvfortran" ]]; then
	sudo apt-get -y install lmod g++ libtinfo5 libncursesw5 lua-posix lua-filesystem lua-lpeg lua-luaossl
	nv_major=22
	nv_minor=3
	nverdot="$nv_major"."$nv_minor"
	nverdash="$nv_major"-"$nv_minor"
	arch_dpkg=`dpkg --print-architecture`
	echo 'deb [trusted=yes] https://developer.download.nvidia.com/hpc-sdk/ubuntu/'$arch_dpkg' /' | sudo tee /etc/apt/sources.list.d/nvhpc.list
	echo '*** added hpc-sdk source to /etc/aps ***'
	ls -lrt /etc/apt/sources.list.d/ || true
	ls -lrt	/etc/apt/sources.list.d/nvhpc.list || true
	sudo cat /etc/apt/sources.list.d/nvhpc.list || true
	sudo apt-get update -y
	apt-cache search nvhpc
        sudo apt-get install -y nvhpc-"$nverdash"
	export PATH=/opt/nvidia/hpc_sdk/Linux_"$arch"/"$nverdot"/compilers/bin:$PATH
	export LD_LIBRARY_PATH=/opt/nvidia/hpc_sdk/Linux_"$arch"/"$nverdot"/compilers/lib:$LD_LIBRARY_PATH
	sudo /opt/nvidia/hpc_sdk/Linux_"$arch"/"$nverdot"/compilers/bin/makelocalrc -x

#	source /etc/profile.d/lmod.sh
#        module use /opt/nvidia/hpc_sdk/modulefiles
#	module load nvhpc
	export FC=nvfortran
#	if [ -z "$BUILD_MPICH" ] ; then
##use bundled openmpi
#	export PATH=/opt/nvidia/hpc_sdk/Linux_"$arch"/"$nverdot"/comm_libs/mpi/bin:$PATH
#	export LD_LIBRARY_PATH=/opt/nvidia/hpc_sdk/Linux_"$arch"/"$nverdot"/comm_libs/mpi/lib:$LD_LIBRARY_PATH
#	fi
	export CC=gcc
	nvfortran -v
	nvfortran
	which nvfortran
    fi
    fi
fi
