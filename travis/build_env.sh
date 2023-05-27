#!/usr/bin/env bash
if [[ -z "$TRAVIS_BUILD_DIR" ]] ; then
    TRAVIS_BUILD_DIR=$(pwd)
fi
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
if [[ $(grep -c fedora /etc/os-release) > 0 ]]; then
    dist="fedora"
fi
echo dist is "$dist"
if [ -z "$DISTR" ] ; then
    DISTR=$dist
fi
echo DISTR is "$DISTR"
	IONEAPI_ROOT=~/apps/oneapi
 if [[ "$os" == "Darwin" ]]; then 
#  HOMEBREW_NO_AUTO_UPDATE=1 brew cask uninstall oclint || true  
#  HOMEBREW_NO_INSTALL_CLEANUP=1  HOMEBREW_NO_AUTO_UPDATE=1 brew install gcc "$MPI_IMPL" openblas python3 ||true
     HOMEBREW_NO_INSTALL_CLEANUP=1  HOMEBREW_NO_AUTO_UPDATE=1 brew install gcc "$MPI_IMPL" python3 gsed grep automake autoconf ||true
     if [[ "$FC" != "gfortran" ]]; then
	 #install non default gfortran, ie gfortran-9
	 #get version
	 mygccver=$(echo "$FC"|cut -d - -f 2)
	 echo mygccver is "$mygccver"
	 HOMEBREW_NO_INSTALL_CLEANUP=1  HOMEBREW_NO_AUTO_UPDATE=1 brew reinstall gcc@"$mygccver" || true
     fi
     #hack to fix Github actions mpif90
     gccver=`brew list --versions gcc| head -1 |cut -c 5-`
     echo brew gccver is $gccver
     if [ -z "$HOMEBREW_CELLAR" ] ; then
	 HOMEBREW_CELLAR=/usr/local/Cellar
     fi
     ln -sf $HOMEBREW_CELLAR/gcc/$gccver/bin/gfortran-* $HOMEBREW_CELLAR/gcc/$gccver/bin/gfortran || true
     ln -sf $HOMEBREW_CELLAR/gcc/$gccver/bin/gfortran-* /usr/local/bin/gfortran || true
     #	 ln -sf /usr/local/bin/$FC /usr/local/bin/gfortran
     $FC --version
     gfortran --version
     if [[ "$FC" == "ifort" ]] || [[ "$FC" == "ifx" ]] ; then
         if [[ -f "$IONEAPI_ROOT"/setvars.sh ]]; then 
	     echo ' using intel cache installation '
	 else
	mkdir -p ~/mntdmg $IONEAPI_ROOT || true
	cd ~/Downloads
	dir_base="2516a0a0-de4d-4f3d-9e83-545b32127dbb"
	dir_hpc="a99cb1c5-5af6-4824-9811-ae172d24e594"
	base="m_BaseKit_p_2023.1.0.45568"
	hpc="m_HPCKit_p_2023.1.0.44543"
	curl -sS -LJO https://registrationcenter-download.intel.com/akdlm/IRC_NAS/"$dir_base"/"$base".dmg
	curl -sS -LJO https://registrationcenter-download.intel.com/akdlm/IRC_NAS/"$dir_hpc"/"$hpc".dmg
	echo "installing BaseKit"
	hdiutil attach "$base".dmg  -mountpoint ~/mntdmg -nobrowse
	sudo  ~/mntdmg/bootstrapper.app/Contents/MacOS/install.sh  -c -s --action install  \
        --components intel.oneapi.mac.mkl.devel  --install-dir $IONEAPI_ROOT --eula accept
	hdiutil detach ~/mntdmg
        #fix slow ifort https://community.intel.com/t5/Intel-oneAPI-HPC-Toolkit/slow-execution-of-ifort-icpc-on-MacOSX-catalina/m-p/1203190
	#
	echo "installing HPCKit"
	hdiutil attach "$hpc".dmg  -mountpoint ~/mntdmg -nobrowse
	sudo  ~/mntdmg/bootstrapper.app/Contents/MacOS/install.sh -c -s  --eula accept \
	     --action install --components default --install-dir $IONEAPI_ROOT
	hdiutil detach ~/mntdmg
	$TRAVIS_BUILD_DIR/travis/fix_xcodebuild.sh
	sudo cp xcodebuild "$IONEAPI_ROOT"/compiler/latest/mac/bin/intel64/.
	ls -lrta $IONEAPI_ROOT ||true
	sudo rm -rf "$IONEAPI_ROOT"/intelpython "$IONEAPI_ROOT"/dal "$IONEAPI_ROOT"/advisor \
	     "$IONEAPI_ROOT"/ipp "$IONEAPI_ROOT"/conda_channel 	"$IONEAPI_ROOT"/dnnl \
	     "$IONEAPI_ROOT"/installer "$IONEAPI_ROOT"/vtune_profiler "$IONEAPI_ROOT"/tbb || true
	fi
	 source "$IONEAPI_ROOT"/setvars.sh --force || true
	 export I_MPI_F90="$FC"
	ls -lrta $IONEAPI_ROOT ||true
	rm -f *dmg || true
	"$FC" -V
	icc -V
     fi
     #hack to get 3.10 as default
     brew install python@3.10
     brew link --force --overwrite python@3.10
     if [[ "$MPI_IMPL" == "mpich" ]]; then
	 #         brew install mpich && brew upgrade mpich && brew unlink openmpi && brew unlink mpich && brew link --overwrite  mpich ||true
	 brew update || true
	 brew unlink open-mpi && brew install mpich && brew upgrade mpich  && brew link --overwrite  mpich || true
     fi
#  if [[ "$MPI_IMPL" == "openmpi" ]]; then
#      HOMEBREW_NO_INSTALL_CLEANUP=1 HOMEBREW_NO_AUTO_UPDATE=1 brew install scalapack
#  fi
fi
if [[ "$os" == "Linux" ]]; then
    if [[ "$DISTR" == "fedora" ]] || [[ "$DISTR" == "centos" ]] ; then
	env
	rpminst=dnf
	if [[ "$DISTR" == "centos" ]] ; then
	    rpminst=yum
	fi
	if [[ "$HOSTNAME" != "fedoraqemuwe40672" ]]; then
	    sudo $rpminst update;  sudo $rpminst -y install perl perl python3-devel time patch cmake gcc-gfortran unzip which make tar bzip2 openssh-clients rsync
	    sudo $rpminst -y install openblas-serial64 || true
	    #	 module load mpi
	    if [[ "$MPI_IMPL" == "openmpi" ]]; then
		sudo $rpminst -y install  openmpi-devel
            elif [[ "$MPI_IMPL" == "mpich" ]]; then
		sudo $rpminst -y install mpich  mpich-devel
	    else
		echo ready only for openmpi
		exit 1
	    fi
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
	    export TERM=dumb
            rm -f l_Base*sh l_HP*sh
	    tries=0 ; until [ "$tries" -ge 10 ] ; do \
			  dir_base="7deeaac4-f605-4bcf-a81b-ea7531577c61"
			  dir_hpc="1ff1b38a-8218-4c53-9956-f0b264de35a4"
			  base="l_BaseKit_p_2023.1.0.46401_offline"
			  hpc="l_HPCKit_p_2023.1.0.46346_offline"
			  wget -nv https://registrationcenter-download.intel.com/akdlm/IRC_NAS/"$dir_hpc"/"$hpc".sh \
			      && wget -nv  https://registrationcenter-download.intel.com/akdlm/IRC_NAS/"$dir_base"/"$base".sh \
			      && break ;\
			      tries=$((tries+1)) ; echo attempt no.  $tries    ; sleep 30 ;  done

	    if [[ "$MPI_IMPL" == "intel" ]]; then
		mpi_bin="  " ; mpi_libdev=" " scalapack_libdev=" "
	    fi
	fi
	if [[ "$GITHUB_WORKFLOW" != "NWChem_CI_selfhosted" ]]; then
	    sudo apt-get update
	    sudo apt-get -y install software-properties-common
	    sudo add-apt-repository universe && sudo apt-get update
	    if [[ "$FC" == "gfortran-11" ]] || [[ "$CC" == "gcc-11" ]]; then
		sudo  add-apt-repository -y ppa:ubuntu-toolchain-r/test 
		pkg_extra+="gcc-11 gfortran-11 g++-11"
	    fi
	    if [[ "$USE_LIBXC" == "-1" ]]; then
		pkg_extra+=" libxc-dev"
	    fi
	    echo pkg to install: gfortran python3-dev  make perl  python3 rsync "$mpi_libdev" "$mpi_bin" $pkg_extra
            tries=0 ; until [ "$tries" -ge 10 ] ; do \
			  sudo apt-get -y install gfortran python3-dev  make perl  python3 rsync "$mpi_libdev" "$mpi_bin" $pkg_extra \
			      && break ;\
			  tries=$((tries+1)) ; echo attempt no.  $tries    ; sleep 30 ;  done

	fi
	if [[ "$FC" == "ifort" ]] || [[ "$FC" == "ifx" ]]; then
            sh ./"$base".sh -a -c -s --action remove --install-dir $IONEAPI_ROOT   --eula accept
            sh ./"$hpc".sh -a -c -s --action remove --install-dir  $IONEAPI_ROOT  --eula accept
	    
            sh ./"$base".sh -a -c -s --action install --components intel.oneapi.lin.mkl.devel --install-dir $IONEAPI_ROOT  --eula accept
	    intel_components="intel.oneapi.lin.ifort-compiler:intel.oneapi.lin.dpcpp-cpp-compiler-pro"
	    if [[ "$MPI_IMPL" == "intel" ]]; then
		intel_components+=":intel.oneapi.lin.mpi.devel"
	    fi
            sh ./"$hpc".sh -a -c -s --action install \
               --components  "$intel_components"  \
               --install-dir $IONEAPI_ROOT     --eula accept
	    rm -f ./"$hpc".sh ./"$base".sh
	    if [[ "$?" != 0 ]]; then
		echo "apt-get install failed: exit code " "${?}"
		exit 1
	    fi
            source "$IONEAPI_ROOT"/setvars.sh || true
	    export I_MPI_F90="$FC"
	    "$FC" -V
	    icc -V

	fi
	if [[ "$FC" == "flang" ]]; then
	    if [[ "USE_AOMP" == "Y" ]]; then
		aomp_major=16
		aomp_minor=0-3
		wget -nv https://github.com/ROCm-Developer-Tools/aomp/releases/download/rel_"$aomp_major"."$aomp_minor"/aomp_Ubuntu2004_"$aomp_major"."$aomp_minor"_amd64.deb
		sudo dpkg -i aomp_Ubuntu2004_"$aomp_major"."$aomp_minor"_amd64.deb
		export PATH=/usr/lib/aomp_"$aomp_major"."$aomp_minor"/bin/:$PATH
		export LD_LIBRARY_PATH=/usr/lib/aomp_"$aomp_major"."$aomp_minor"/lib:$LD_LIBRARY_PATH
		ls -lrt /usr/lib | grep aomp ||true
	    else
		aocc_version=4.0.0
		aocc_dir=aocc-compiler-${aocc_version}
#		curl -sS -LJO https://developer.amd.com/wordpress/media/files/${aocc_dir}.tar
		tries=0 ; until [ "$tries" -ge 10 ] ; do \
                curl -sS -LJO https://download.amd.com/developer/eula/aocc-compiler/${aocc_dir}.tar \
                && break ; \
                tries=$((tries+1)) ; echo attempt no.  $tries    ; sleep 30 ;  done
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
	    rocm_version=5.4.3
	    wget -q -O - https://repo.radeon.com/rocm/rocm.gpg.key |  sudo apt-key add -
	    echo 'deb [arch=amd64] https://repo.radeon.com/rocm/apt/'$rocm_version'/ ubuntu main' | sudo tee /etc/apt/sources.list.d/rocm.list
	    sudo apt-get  update -y && sudo apt-get -y install rocm-llvm openmp-extras
	    export PATH=/opt/rocm/bin:$PATH
	    export LD_LIBRARY_PATH=/opt/rocm/lib:/opt/rocm/llvm/lib:$LD_LIBRARY_PATH
	    amdflang -v
	    amdclang -v
	fi
	if [[ "$FC" == "nvfortran" ]]; then
	    sudo apt-get -y install lmod g++ libtinfo5 libncursesw5 lua-posix lua-filesystem lua-lpeg lua-luaossl
	    nv_major=23
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
	    
	    export FC=nvfortran
	    export CC=gcc
	    nvfortran -V
	    which nvfortran
	fi
    fi
    # check for mpif90 command and exit if not present
    if [[ ! $(command -v mpif90) ]]; then echo "mpif90 not present"; exit 1; fi
    echo "mpif90 -show output is " `mpif90 -show` || true
    echo "which mpif90 output is " `which mpif90` ||  true
fi
