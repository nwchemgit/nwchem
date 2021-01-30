#!/bin/bash
os=`uname`
arch=`uname -m`
dist="ubuntu"
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
 if [[ "$os" == "Darwin" ]]; then 
#  HOMEBREW_NO_AUTO_UPDATE=1 brew cask uninstall oclint || true  
#  HOMEBREW_NO_INSTALL_CLEANUP=1  HOMEBREW_NO_AUTO_UPDATE=1 brew install gcc "$MPI_IMPL" openblas python3 ||true
     HOMEBREW_NO_INSTALL_CLEANUP=1  HOMEBREW_NO_AUTO_UPDATE=1 brew install gcc "$MPI_IMPL" python3 ||true
     #hack to fix Github actions mpif90
     ln -sf /usr/local/bin/$FC /usr/local/bin/gfortran
     $FC --version
     gfortran --version
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
    if [[ "$MPI_IMPL" == "intel" ]]; then
	export APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=1
        tries=0 ; until [ "$tries" -ge 5 ] ; do \
	wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
            && sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB  \
            && rm -f GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB || true \
            && echo "deb https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list \
            && sudo add-apt-repository "deb https://apt.repos.intel.com/oneapi all main"  \
	    && sudo apt-get update && break ;\
            tries=$((tries+1)) ; echo attempt no.  $tries    ; sleep 15 ;  done

	sudo apt-cache search intel-oneapi-mpi
        mpi_bin="  " ; mpi_libdev="intel-oneapi-mpi-devel" scalapack_libdev="intel-oneapi-mkl"
    fi
    sudo add-apt-repository universe && sudo apt-get update
#    sudo apt-get -y install gfortran python3-dev python-dev cmake "$mpi_libdev" "$mpi_bin" "$scalapack_libdev"  make perl  libopenblas-dev python3 rsync
    sudo apt-get -y install gfortran python3-dev python-dev cmake "$mpi_libdev" "$mpi_bin"  make perl  python3 rsync
    if [[ "$FC" == "ifort" ]]; then
	sudo apt-get -y install intel-oneapi-ifort intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic  intel-oneapi-mkl
	sudo apt-get -y install intel-oneapi-mpi-devel
    fi
    if [[ "$FC" == "flang" ]]; then
	wget https://github.com/ROCm-Developer-Tools/aomp/releases/download/rel_11.12-0/aomp_Ubuntu2004_11.12-0_amd64.deb
	sudo dpkg -i aomp_Ubuntu2004_11.12-0_amd64.deb
	export PATH=/usr/lib/aomp_11.12-0/bin/:$PATH
    fi
    fi
fi
